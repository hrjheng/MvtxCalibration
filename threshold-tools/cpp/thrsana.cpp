#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <limits>

#include <immintrin.h>
#include <mm_malloc.h>
#include <signal.h>
#include <unistd.h>

#include "progressbar.h"
#include "utils.h"

using namespace std;

int runnumber = 0;

constexpr uint8_t FLX_WORD_SZ = 32;

struct decoder_t
{
    std::ifstream file;

    static constexpr uint8_t nlanes = 3;

    uint32_t filebuf_size;
    uint32_t lanebuf_size;
    uint32_t nbytesleft;
    size_t nbytesread;
    uint8_t *filebuffer;
    uint8_t *lanebuffer;

    std::unordered_map<uint64_t, size_t> m_strobeHits = {};
    std::vector<uint32_t> m_hits = {};

    std::array<uint8_t *, nlanes> lane_ends;

    uint8_t *ptr;
    uint8_t *prev_packet_ptr;

    std::vector<uint32_t> feeids;
    std::array<uint8_t, nlanes> chipIds = {0xFF, 0xFF, 0xFF};

    trg_t trigger;

    int thscan_nInj;
    int thscan_nChg;
    int feeid = -1;

    bool isThr = false;
    bool isTun = false;

    uint32_t evt_cnts[Trg::BitMap::nBitMap] = {};

    decoder_t(std::string &flName)
        : filebuf_size(8192 * 1001)
        , lanebuf_size(1024 * 1024)
        , nbytesleft(0)
        , nbytesread(0)
        , thscan_nInj(0)
        , thscan_nChg(0)
    {
        // TODO: check meaningfullness of params...
        file.open(flName, std::ios::binary);
        if (!file.is_open())
        {
            std::cerr << "Error while trying to open file: " << strerror(errno);
            std::cerr << ". Exiting" << std::endl;
            exit(-1);
        }

        assert(!(filebuf_size % 8192));
        // memory allocations
        filebuffer = reinterpret_cast<uint8_t *>(_mm_malloc(filebuf_size * sizeof(uint8_t), 4096));
        lanebuffer = reinterpret_cast<uint8_t *>(_mm_malloc(nlanes * lanebuf_size * sizeof(uint8_t), 4096));
        m_hits.resize(100000);

        if (!filebuffer || !lanebuffer)
        {
            std::cerr << "Error while allocating memory: " << strerror(errno);
            std::cerr << ". Exiting." << std::endl;
            exit(-1);
        }
        ptr = filebuffer;
        trigger = {};
    }

    ~decoder_t()
    {
        _mm_free(filebuffer);
        _mm_free(lanebuffer);
    }

    inline void packet_init(bool reset_hits = false)
    {
        for (uint8_t i = 0; i < nlanes; ++i)
        {
            lane_ends[i] = lanebuffer + i * lanebuf_size;
        }
        if (reset_hits)
        {
            clearVec(m_hits);
            m_strobeHits.clear();
        }
    }

    inline bool has_lane_data()
    {
        uint8_t i = 0;
        for (const auto &lane_end : lane_ends)
        {
            if (lane_end != (lanebuffer + i * lanebuf_size))
            {
                return true;
            }
            ++i;
        }
        return false;
    }

    size_t ptr_pos(uint8_t *_ptr = nullptr) { return (((!_ptr) ? ptr : _ptr) - filebuffer); }
};

uint32_t nHB = 0, nHB_with_data = 0, nStrobe = 0, nPhysTrg_in_HB = 0;
std::vector<uint64_t> physTrg_bco;

void reset_stat()
{
    nHB = 0;
    nHB_with_data = 0;
    nStrobe = 0;
}

void printStat(const size_t n_events, const size_t n_evt_with_payload, const size_t nTrg)
{
    std::cout << "Read " << n_events << " events. " << n_evt_with_payload;
    std::cout << " with ALPIDE payload and " << nTrg << " strobes" << std::endl;
    std::cout << std::endl;
}

void printTrgCnts(decoder_t *decoder)
{
    std::cout << "Trigger counts:" << endl;
    for (const auto &trg : Trg::allBitMap)
    {
        std::cout << Trg::BitMapName[trg].c_str() << ": " << decoder->evt_cnts[trg] << std::endl;
    }
    std::cout << std::endl;
}

void save_block(decoder_t *decoder)
{
    ostringstream ss;
    ss << "last_chunk_" << decoder->feeid << ".err";
    ofstream fdump(ss.str().c_str(), ios::trunc);
    fdump.write(reinterpret_cast<char *>(decoder->filebuffer), decoder->ptr - decoder->filebuffer + decoder->nbytesleft);
    fdump.close();
}

void meanrms(float *m, float *s, float *data, size_t n)
{
    float s1 = 0.;
    float s2 = 0.;
    int nn0 = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (data[i] > 0)
        {
            float x = data[i];
            s1 += x;
            s2 += x * x;
            ++nn0;
        }
    }
    s1 /= nn0;
    s2 /= nn0;
    *m = s1;
    *s = sqrtf(s2 - s1 * s1);
}

static inline void threshold_next_row(float *thrs, float *rmss, float *sumd, float *sumd2, int nch, int ninj)
{
    for (int x = 0; x < 3 * 1024; ++x)
    {
        float den = sumd[x];
        float m = sumd2[x];
        if (den > 0)
        {
            m /= den;
        }
        float u = den;
        float s = sqrtf(m - u * u);
        thrs[x] = m;
        rmss[x] = s;
    }
}

static inline void threshold_next_charge(float *sumd, float *sumd2, int ch, int *lasthist, int *hist, int ninj)
{
    int ch1 = ch - 1;
    int ch2 = ch;

    float ddV = 1.0 * (ch2 - ch1);
    float V1 = ch1;
    float V2 = ch2;
    float meandV = 0.5 * (V2 + V1);

    for (int x = 0; x < 3 * 1024; ++x)
    {
        float f = 1.f / (1.f * ninj);
        float n2 = hist[x] * f;
        float n1 = lasthist[x] * f;

        float dn = n2 - n1;
        float den = dn / ddV;
        float m = meandV * dn / ddV;
        sumd[x] += den;
        sumd2[x] += m;
    }
}

// temporal-noise
static inline void tnoise_accumulate(float *S0,           //
                                     float *S1,           //
                                     float *S2,           //
                                     int ch,              //
                                     const int *lasthist, //
                                     const int *hist,     //
                                     int ninj,            //
                                     int npix /* = nChipsPerLane*1024 */)
{
    const float invN = 1.0f / static_cast<float>(ninj);
    const float c_mid = static_cast<float>(ch) - 0.5f;

    for (int x = 0; x < npix; ++x)
    {
        float dn = static_cast<float>(hist[x] - lasthist[x]) * invN;
        if (dn <= 0.0f)
            continue;                // suppress small negative fluctuations
        S0[x] += dn;                 // total weight (area of derivative)
        S1[x] += dn * c_mid;         // first moment
        S2[x] += dn * c_mid * c_mid; // second moment
    }
}

static inline void tnoise_finalize_row(float *thr_mu_map, //
                                       float *tnoise_map, //
                                       const float *S0,   //
                                       const float *S1,   //
                                       const float *S2,   //
                                       int row,           //
                                       int nChipsPerLane  //
)
{
    const int stride = nChipsPerLane * 1024;
    float *mu_row = thr_mu_map + row * stride;
    float *sig_row = tnoise_map + row * stride;

    for (int x = 0; x < stride; ++x)
    {
        const float s0 = S0[x];
        if (s0 > 0.0f)
        {
            const float mu = S1[x] / s0;
            float var = (S2[x] / s0) - mu * mu;
            if (var < 0.0f)
                var = 0.0f; // numerical safety
            mu_row[x] = mu;
            sig_row[x] = std::sqrt(var);
        }
        else
        {
            mu_row[x] = std::numeric_limits<float>::quiet_NaN();
            sig_row[x] = std::numeric_limits<float>::quiet_NaN();
        }
    }
}

static inline void fillrowhist(int *hist, uint32_t *hits, int n, int row)
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverflow"
    const __m128i masky = _mm_set1_epi32(0x000003FE);
    const __m128i mask0 = _mm_set1_epi32(0x00000001);
    const __m128i maska = _mm_set1_epi32(0x0000FFFE);
#pragma GCC diagnostic pop

    const __m128i bad = _mm_set1_epi32(3 * 1024); // bad pixel ID
    __m128i y = _mm_set1_epi32(row << 1);
    while (n > 0)
    {
        __m128i a = _mm_load_si128((__m128i *)hits);
        __m128i yp = _mm_and_si128(a, masky);
        __m128i ok = _mm_cmpeq_epi32(yp, y);

        __m128i t1 = _mm_srli_epi32(a, 1);
        __m128i t2 = _mm_srli_epi32(a, 9);
        __m128i t3 = _mm_and_si128(t2, maska);
        __m128i t4 = _mm_xor_si128(t1, a);
        __m128i t5 = _mm_and_si128(t4, mask0);
        __m128i x = _mm_or_si128(t3, t5);
        __m128i b = _mm_blendv_epi8(bad, x, ok); // epi8 is OK, due to ok

        if (n > 3)
            ++hist[_mm_extract_epi32(b, 3)];
        if (n > 2)
            ++hist[_mm_extract_epi32(b, 2)];
        if (n > 1)
            ++hist[_mm_extract_epi32(b, 1)];

        ++hist[_mm_extract_epi32(b, 0)];

        hits += 4;
        n -= 4;
    }
}

static inline void transformhits(uint32_t *hits, int n)
{
    const __m128i masky = _mm_set1_epi32(0x000001FF);
    const __m128i mask0 = _mm_set1_epi32(0x00000001);
    const __m128i maska = _mm_set1_epi32(0x0000FFFE);
    while (n > 0)
    {
        __m128i a = _mm_load_si128((__m128i *)hits);

        __m128i t1 = _mm_srli_epi32(a, 1);    // >> 1 each 4 32'b elements
        __m128i y = _mm_and_si128(t1, masky); // & each 4 32'b element with 0x000001FF

        __m128i t2 = _mm_srli_epi32(a, 9);
        __m128i t3 = _mm_and_si128(t2, maska);
        __m128i t4 = _mm_xor_si128(t1, a);
        __m128i t5 = _mm_and_si128(t4, mask0);
        __m128i x = _mm_or_si128(t3, t5);

        __m128i t6 = _mm_slli_epi32(y, 16);
        __m128i yx = _mm_or_si128(t6, x);

        _mm_store_si128((__m128i *)hits, yx);

        hits += 4;
        n -= 4;
    }
}

static inline void fillhitmap(uint32_t *map, uint32_t *hits, int n)
{
    const __m128i masky __attribute__((unused)) = _mm_set1_epi32(0x01FF0000);
    const __m128i maskx = _mm_set1_epi32(0x00003FFF);
    const __m128i stride = _mm_set1_epi32(3 * 1024);
    while (n > 0)
    {
        __m128i hs = _mm_stream_load_si128((__m128i *)hits);
        __m128i t1 = _mm_srli_epi32(hs, 16);
        __m128i t2 = _mm_mullo_epi32(t1, stride);
        __m128i t3 = _mm_and_si128(hs, maskx);
        __m128i a = _mm_add_epi32(t2, t3);
        if (n > 3)
            ++map[_mm_extract_epi32(a, 3)];
        if (n > 2)
            ++map[_mm_extract_epi32(a, 2)];
        if (n > 1)
            ++map[_mm_extract_epi32(a, 1)];

        ++map[_mm_extract_epi32(a, 0)];

        hits += 4;
        n -= 4;
    }
}

void check_APE(const uint8_t &laneptr)
{
    switch (laneptr)
    {
    case 0xF2:
        std::cerr << " APE_STRIP_START" << std::endl;
        break;
    case 0xF4:
        std::cerr << " APE_DET_TIMEOUT" << std::endl;
        break;
    case 0xF5:
        std::cerr << " APE_OOT" << std::endl;
        break;
    case 0xF6:
        std::cerr << " APE_PROTOCOL_ERROR" << std::endl;
        break;
    case 0xF7:
        std::cerr << " APE_LANE_FIFO_OVERFLOW_ERROR" << std::endl;
        break;
    case 0xF8:
        std::cerr << " APE_FSM_ERROR" << std::endl;
        break;
    case 0xF9:
        std::cerr << " APE_PENDING_DETECTOR_EVENT_LIMIT" << std::endl;
        break;
    case 0xFA:
        std::cerr << " APE_PENDING_LANE_EVENT_LIMIT" << std::endl;
        break;
    case 0xFB:
        std::cerr << " APE_O2N_ERROR" << std::endl;
        break;
    case 0xFC:
        std::cerr << " APE_RATE_MISSING_TRG_ERROR" << std::endl;
        break;
    case 0xFD:
        std::cerr << " APE_PE_DATA_MISSING" << std::endl;
        break;
    case 0xFE:
        std::cerr << " APE_OOT_DATA_MISSING" << std::endl;
        break;
    default:
        std::cerr << " Unknown APE code" << std::endl;
    }
    return;
}

static inline void decode(uint8_t *laneptr, uint8_t *laneend, decoder_t *decoder, uint8_t &chipId)
{
    if (laneptr == laneend)
    {
        return;
    }

    uint8_t busy_on = 0;
    uint8_t busy_off = 0;

    uint8_t laneId = 0;
    uint8_t reg = 0;
    uint8_t chip_header_found = 0;
    uint8_t chip_trailer_found = 0;

    while (laneptr < laneend)
    {
        // TODO: check out of bounds problem (better: ensure that the 2 bytes following
        // laneend are readable)
        if (*laneptr == 0xF1)
        { // BUSY ON
            ++busy_on;
            ++laneptr;
        }
        else if (*laneptr == 0xF0)
        { // BUSY OFF
            ++busy_off;
            ++laneptr;
        }
        else if ((*laneptr & 0xF0) == 0xE0)
        { // EMPTY
            chip_header_found = 0;
            chip_trailer_found = 1;
            chipId = laneptr[0] & 0xF;
            laneId = chipId % 3;
            // abc = laneptr[1];
            busy_on = busy_off = 0;
            laneptr += 2;
        }
        else
        {
            if (chip_header_found)
            {
                if ((laneptr[0] & 0xE0) == 0xC0)
                { // REGION HEADER
                    // TODO: move first region header out of loop, asserting its existence
                    reg = laneptr[0] & 0x1F;
                    ++laneptr;
                }
                if ((laneptr[0] & 0xC0) == 0x40)
                { // DATA SHORT
                    int addr = (laneptr[0] & 0x3F) << 8 | laneptr[1];
                    addr |= (laneId << 19) | (reg << 14);
                    decoder->m_hits.push_back(addr);
                    laneptr += 2;
                }
                else if ((laneptr[0] & 0xC0) == 0x00)
                { // DATA LONG
                    int addr = (laneptr[0] & 0x3F) << 8 | laneptr[1];
                    addr |= (laneId << 19) | (reg << 14);
                    decoder->m_hits.push_back(addr);

                    uint8_t hitmap = laneptr[2]; // TODO: assert that bit 8 is 0?
                    if (hitmap & 0x01)
                        decoder->m_hits.push_back(addr + 1);
                    if (hitmap & 0x7E)
                    { // provide early out (mostly 2-pixel clusters...)
                        if (hitmap & 0x02)
                            decoder->m_hits.push_back(addr + 2);
                        if (hitmap & 0x04)
                            decoder->m_hits.push_back(addr + 3);
                        if (hitmap & 0x08)
                            decoder->m_hits.push_back(addr + 4);
                        if (hitmap & 0x10)
                            decoder->m_hits.push_back(addr + 5);
                        if (hitmap & 0x20)
                            decoder->m_hits.push_back(addr + 6);
                        if (hitmap & 0x40)
                            decoder->m_hits.push_back(addr + 7);
                    }
                    laneptr += 3;
                }
                else if ((laneptr[0] & 0xF0) == 0xB0)
                { // CHIP TRAILER
                    chip_trailer_found = 1;
                    busy_on = busy_off = chip_header_found = 0;
                    ++laneptr;
                }
                else if (laneptr[0] == 0xF0)
                { // BUSY_OFF
                    ++laneptr;
                    ++busy_off;
                }
                else if (laneptr[0] == 0xF1)
                { // BUSY_ON
                    ++laneptr;
                    ++busy_on;
                }
                else if ((laneptr[0] & 0xF0) == 0xF0)
                { // APE
                    std::cerr << " Chip " << (int)chipId << ":";
                    check_APE(laneptr[0]);
                    chip_trailer_found = 1;
                    busy_on = busy_off = chip_header_found = 0;
                    ++laneptr;
                }
                else
                {
                    // ERROR
                    std::cerr << "ERROR: invalid byte 0x" << std::hex << (int)(laneptr[0]) << std::endl;
                    while (laneptr != laneend)
                    {
                        std::cerr << " " << std::hex << (int)(*(uint8_t *)laneptr) << " ";
                        // printf(" %02X ", *(uint8_t*)laneptr);
                        ++laneptr;
                    }
                }
            }
            else
            { // chip_header
                if ((laneptr[0] & 0xF0) == 0xA0)
                {
                    chip_header_found = 1;
                    chip_trailer_found = 0;
                    chipId = laneptr[0] & 0xF;
                    laneId = chipId % 3;
                    // abc = laneptr[1];
                    reg = 0;
                    laneptr += 2;
                }
                else if (laneptr[0] == 0x00)
                { // padding
                    ++laneptr;
                }
                else if (laneptr[0] == 0xf0)
                { // BUSY_OFF
                    ++laneptr;
                    ++busy_off;
                }
                else if (laneptr[0] == 0xf1)
                { // BUSY_ON
                    ++laneptr;
                    ++busy_on;
                }
                else if ((laneptr[0] & 0xF0) == 0xF0)
                { // APE
                    std::cerr << " Chip " << (int)chipId << ":";
                    check_APE(laneptr[0]);
                    chip_trailer_found = 1;
                    busy_on = busy_off = chip_header_found = 0;
                    ++laneptr;
                }
                else
                { // ERROR
                    std::cerr << "ERROR: invalid byte 0x" << std::hex << (int)(laneptr[0]) << std::endl;
                    while (laneptr != laneend)
                    {
                        std::cerr << " " << std::hex << (int)(*(uint8_t *)laneptr) << " ";
                        ++laneptr;
                    }
                }
            } // data
        } // busy_on, busy_off, chip_empty, other
    } // while
    if (!chip_trailer_found)
    {
        std::cerr << "\nERROR: ALPIDE data end without data trailer.";
        std::cerr << " HB: " << nHB << ", Strobe: " << decoder->m_strobeHits.size() << std::endl;
    }
    return;
}

static inline void decoder_decode_lanes_into_hits(decoder_t *decoder)
{
    uint8_t chipId = 0xFF;
    size_t strobe_nhits = decoder->m_hits.size();
    for (int i = 0; i < decoder->nlanes; ++i)
    {
        if (decoder->lanebuffer + decoder->lanebuf_size * i != decoder->lane_ends[i])
        {
            decode(decoder->lanebuffer + decoder->lanebuf_size * i, decoder->lane_ends[i], decoder, chipId);
            assert(!(decoder->trigger.strobe_id >> 52));
            //     assert(chipId != 0xFF);
            decoder->chipIds[i] = chipId;
        }
    }
    decoder->m_strobeHits[decoder->trigger.strobe_id] = decoder->m_hits.size() - strobe_nhits;
}

size_t pull_data(decoder_t *decoder)
{
    ssize_t nread = 0;
    uint8_t *buf_ptr = decoder->filebuffer; // set buf_ptr to filebuffer init
    decoder->nbytesread += decoder->ptr - buf_ptr;
    memcpy(buf_ptr, decoder->ptr, decoder->nbytesleft); // move byte left to filebuffer start
    decoder->ptr = buf_ptr;
    buf_ptr += decoder->nbytesleft;
    size_t len = (decoder->filebuf_size - decoder->nbytesleft);
    len = (len > 0x1000) ? len & ~0xFFF : len; // in chunks of 256 bytes multiples
    uint64_t n;
    do
    {

        decoder->file.read((char *)buf_ptr, len);
        n = decoder->file.gcount();
        len -= n;
        buf_ptr += n;
        nread += n;
    } while (len > 0 && n > 0);

    return nread;
}

enum EXIT_CODE
{
    NO_FLX_HEADER = -1,
    BAD_READ = -2,
    BAD_END_FILE = -3,
    BAD_GBT_CNT = -4,
    DONE = 0,
    HB_DATA_DONE = 1,
    HB_NO_DATA_DONE = 2,
    N_EXIT_CODE
};

static inline EXIT_CODE decoder_read_event_into_lanes(decoder_t *decoder)
{
    rdh_t rdh;
    uint32_t nStopBit = 0;
    uint32_t prev_pck_cnt = 0;

    bool header_found = false;
    bool prev_evt_complete = false;
    bool haspayload = false;

    size_t n_no_continuation = 0, n_packet_done = 0;

    // loop over pages in files
    while (true)
    {
        if (decoder->nbytesleft > 0)
        {
            bool padding_found = false;
            while ((*(reinterpret_cast<uint16_t *>(decoder->ptr + 30)) == 0xFFFF) && decoder->nbytesleft)
            {
                padding_found = true;
                decoder->ptr += FLX_WORD_SZ;
                decoder->nbytesleft -= FLX_WORD_SZ;
                decoder->nbytesleft += (!decoder->nbytesleft) ? pull_data(decoder) : 0;
            }
            if (padding_found)
            {
                padding_found = false;
                ASSERT(!((decoder->nbytesread + decoder->ptr_pos()) & 0xFF), decoder, "FLX header is not properly aligned in byte %lu of current chunk, previous packet %ld", decoder->ptr_pos(), decoder->ptr_pos(decoder->prev_packet_ptr));
            }
        }

        uint32_t pagesize = 0;
        if (decoder->nbytesleft >= (2 * FLX_WORD_SZ))
        { // at least FLX header and RDH
            if (*(reinterpret_cast<uint16_t *>(decoder->ptr + 30)) == 0xAB01)
            {
                // rdh.decode(decoder->ptr);

                if (!rdh.decode(decoder->ptr))
                {
                    std::cerr << "Failed to decode RDH - invalid rdhGBTcounter" << std::endl;
                    return BAD_GBT_CNT;
                }
                pagesize = (rdh.pageSize + 1) * FLX_WORD_SZ;
            }
            else
            {

                return NO_FLX_HEADER;
            }
        }

        if (!pagesize || decoder->nbytesleft < pagesize)
        { // pagesize = 0 read rdh
            // at least the RDH needs to be there...
            if (decoder->nbytesleft < 0)
            {
                printf("ERROR: d_nbytesleft: %d, less than zero \n", decoder->nbytesleft);
                return BAD_READ;
            }

            size_t nread = pull_data(decoder);

            decoder->nbytesleft += nread;

            if (!nread)
                return (decoder->nbytesleft < pagesize) ? BAD_END_FILE : DONE;

            continue;
        }

        decoder->prev_packet_ptr = decoder->ptr;
        uint8_t *flx_ptr = decoder->ptr; // payload: TODO: check header
        decoder->ptr += pagesize;
        decoder->nbytesleft -= pagesize;

        if (!pagesize)
            continue; // TODO...

        if (decoder->feeid < 0)
        {
            if (rdh.stopBit)
            {
                nStopBit++;
            }

            if (std::find(decoder->feeids.cbegin(), decoder->feeids.cend(), rdh.feeId) == decoder->feeids.cend())
            {
                decoder->feeids.push_back(rdh.feeId);
            }
            else if (nStopBit > 10 * decoder->feeids.size())
            {
                return DONE;
            }
            else
            {
                continue;
            }
        }

        if (rdh.feeId != decoder->feeid)
            continue; // TODO: ...

        flx_ptr += 2 * FLX_WORD_SZ; // skip RDH
        // printf("flx_ptr: %d \n", *(uint8_t*)flx_ptr);

        const size_t nFlxWords = (pagesize - (2 * FLX_WORD_SZ)) / FLX_WORD_SZ;
        // TODO assert pagesize > 2
        ASSERT(((!rdh.packetCounter) || (rdh.packetCounter == prev_pck_cnt + 1)), decoder, "Incorrect pages count %d in byte %ld of current chunk", rdh.packetCounter, decoder->ptr_pos());
        prev_pck_cnt = rdh.packetCounter;

        if (!rdh.packetCounter) // begin of HB
        {
            nHB++;
            nPhysTrg_in_HB = 0;

            clearVec(physTrg_bco);

            // TODO check right protocol for SOC/EOC or SOT/EOT
            decoder->packet_init(true);
            for (const auto &trg : Trg::allBitMap)
            {
                if ((rdh.trgType >> trg) & 1)
                {
                    decoder->evt_cnts[trg]++;
                }
            }
        }
        else
        {
            if (!rdh.stopBit)
            {
                ASSERT(!prev_evt_complete, decoder, "Previous event was already complete, byte %ld of current chuck", decoder->ptr_pos());
            }
        }

        struct tdh_t tdh;
        struct tdt_t tdt;
        struct cdw_t cdw;
        int prev_gbt_cnt = 3;
        for (size_t iflx = 0; iflx < nFlxWords; ++iflx)
        {
            __m256i data = _mm256_stream_load_si256((__m256i *)flx_ptr);
            const uint16_t gbt_cnt = _mm256_extract_epi16(data, 15) & 0x3FF;
            ASSERT((gbt_cnt - prev_gbt_cnt) <= 3, decoder, "Error. Bad gbt counter in the flx packet at byte %ld", decoder->ptr_pos(flx_ptr));
            prev_gbt_cnt = gbt_cnt;
            const uint16_t n_gbt_word = ((gbt_cnt - 1) % 3) + 1;

            uint8_t *gbt_word;
            for (size_t igbt = 0; igbt < n_gbt_word; ++igbt)
            {
                gbt_word = flx_ptr + (igbt * 10);
                uint8_t lane = *(reinterpret_cast<uint8_t *>(gbt_word + 9));

                if (lane == 0xE0)
                {
                    // lane heder: needs to be present: TODO: assert this
                    // TODO assert first word after RDH and active lanes
                    haspayload = false;
                }
                else if (lane == 0xE8)
                { // TRIGGER DATA HEADER (TDH)
                    tdh.decode(gbt_word);
                    header_found = true;
                    uint64_t strobe_id = ((tdh.orbit << 12) | tdh.bc);
                    if (tdh.bc) // statistic trigger for first bc already filled on RDH
                    {
                        for (const auto &trg : Trg::allBitMap)
                        {
                            if (trg == Trg::BitMap::FE_RST) //  TDH save first 12 bits only
                                break;
                            if (((tdh.trigger_type >> trg) & 1))
                            {
                                decoder->evt_cnts[trg]++;
                            }
                        }
                    }

                    if ((tdh.trigger_type >> Trg::PHYSICS) & 0x1) // LL1 in MVTX
                    {
                        ++nPhysTrg_in_HB;
                        physTrg_bco.push_back(tdh.orbit);
                    }

                    if (!tdh.continuation)
                    {
                        n_no_continuation++;
                        decoder->trigger.strobe_id = strobe_id;
                    } // end if not cont
                } // end TDH
                else if (lane == 0xF8)
                { // CALIBRATION DATA WORD (CDW)
                    cdw.decode(gbt_word);
                    if (decoder->isThr)
                    {
                        uint16_t new_row = cdw.user_field & 0xFFFF;
                        uint16_t new_charge = (cdw.user_field >> 16) & 0xFFFF;
                        if ((!decoder->trigger.thscan_inj) || (decoder->trigger.thscan_inj == decoder->thscan_nInj))
                        {
                            if (!decoder->isTun)
                            {
                                ASSERT(new_row >= decoder->trigger.thscan_row, decoder, "Row not increasing after thscan_injections: previous %d new %d", decoder->trigger.thscan_row, new_row);
                            }
                            if (new_row == decoder->trigger.thscan_row)
                            { // rolling charge at change of row
                                ASSERT(new_charge > decoder->trigger.thscan_chg, decoder, "Charge not increasing after max thscan_injections: previous %d, new %d [previous row %d current row %d]", decoder->trigger.thscan_chg, new_charge, decoder->trigger.thscan_row, new_row);
                            }
                            //       else
                            //       {
                            //         if(thscan_current_charge != -1 and rdh['triggers'] != ['EOT'])
                            //         {
                            //          assert new_charge < thscan_current_charge,
                            //             "Charge not decreasing after row change: previous {thscan_current_charge}
                            //            new {new_charge}\nRDH: {rdh}"
                            //         }
                            //       }
                            decoder->trigger.thscan_inj = 1;
                            decoder->trigger.thscan_row = new_row;
                            decoder->trigger.thscan_chg = new_charge;
                        }
                        else
                        {
                            ASSERT(new_row == decoder->trigger.thscan_row, decoder, "Row not correct before reaching max thscan_injections: expected %d got %d", decoder->trigger.thscan_row, new_row);
                            ASSERT(new_charge == decoder->trigger.thscan_chg, decoder, "Charge not correct before reaching max thscan_injections: expected %d got %d, [previous row %d, current row %d observed injections %d]", decoder->trigger.thscan_chg, new_charge, decoder->trigger.thscan_row, new_row, decoder->trigger.thscan_inj);
                            decoder->trigger.thscan_inj += 1;
                        }
                    }
                }
                else if (lane == 0xF0)
                { // lane trailer
                    tdt.decode(gbt_word);
                    prev_evt_complete = tdt.packet_done;
                    if (tdt.packet_done)
                    {
                        n_packet_done++;
                        // assert(n_no_continuation == n_packet_done);
                    }
                }
                else if (lane == 0xE4)
                { // DIAGNOSTIC DATA WORD (DDW)
                    // TODO add diagnostic dataword decoder
                    assert(rdh.stopBit);
                }
                else if (((lane >> 5) & 0x7) == 0x5)
                { // IB DIAGNOSTIC DATA
                    // decode IB diagnostic word
                    cerr << "WARNING!!! IB diagnostic data word received and skipped." << endl;
                }
                else
                { // lane payload
                    ASSERT(((lane >> 5) & 0x7) == 0x1, decoder, "Wrong GBT Word %x in byte %ld, it is not an IB word", lane, decoder->ptr_pos(flx_ptr));
                    haspayload = true;
                    ASSERT(header_found, decoder, "Trigger header not found before chip data, in byte %ld", decoder->ptr_pos(flx_ptr));
                    lane &= 0x1F; // TODO: assert range + map IDs
                    lane %= 3;
                    memcpy(decoder->lane_ends[lane], gbt_word, 9);
                    decoder->lane_ends[lane] += 9;
                }

                if (prev_evt_complete)
                {
                    if (decoder->has_lane_data())
                    {
                        decoder_decode_lanes_into_hits(decoder);
                        decoder->packet_init();
                    }
                    prev_evt_complete = false;
                    header_found = false;
                }
            } // for igbt
            flx_ptr += FLX_WORD_SZ;
        }

        if (rdh.stopBit)
        {
            return (haspayload) ? HB_DATA_DONE : HB_NO_DATA_DONE;
        }
    } // while(true)
}

static inline EXIT_CODE check_next_event(decoder_t *decoder)
{
    EXIT_CODE ret = decoder_read_event_into_lanes(decoder);

    switch (ret)
    {
    case NO_FLX_HEADER:
        ASSERT(false, decoder, "Error reading file, wrong felix header position in byte %ld, previous flx header in byte %ld", decoder->ptr_pos(), decoder->ptr_pos(decoder->prev_packet_ptr));
        break;
    case BAD_READ:
        std::cerr << "Error while reading file BAD_READ: " << strerror(errno);
        std::cerr << " Exiting." << std::endl;
        break;
    case BAD_END_FILE:
        std::cerr << "Error while reading file: " << strerror(errno);
        std::cerr << ". Last read was incomplete. Exiting (some events might be ignored).";
        std::cerr << std::endl;
        break;
    default:;
    };
    return ret;
}

void save_file(std::string fname, uint32_t n_row, uint32_t nChipsPerLane, float *&data)
{
    ofstream fileMap(fname.data(), ios_base::trunc);
    fileMap.write(reinterpret_cast<const char *>(data), n_row * nChipsPerLane * 1024 * sizeof(float));
    fileMap.close();
}

void run_thrs_tuning(struct decoder_t *decoder, string &prefix, const int &n_vcasn_ithr = 1)
{

    if (!decoder->isTun)
    {
        std::cout << "Runing THR analysis for feeid " << decoder->feeid << "..." << std::endl;
    }
    else
    {
        std::cout << "Runing THR tunning analysis for feeid " << decoder->feeid;
        if (n_vcasn_ithr < 0)
        {
            std::cout << ", expected EOX event!!!" << std::endl;
        }
        else
        {
            std::cout << " and n " << n_vcasn_ithr << "..." << std::endl;
        }
    }

    constexpr uint8_t nChipsPerLane = 3;
    int16_t n_row = decoder->isTun ? 6 : 512;
    float *thrs = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * n_row * 1024 * sizeof(float), 4096));
    float *rmss = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * n_row * 1024 * sizeof(float), 4096));
    float *sumd = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * 1 * 1024 * sizeof(float), 4096));
    float *sumd2 = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * 1 * 1024 * sizeof(float), 4096));
    int *rowhist = reinterpret_cast<int *>(_mm_malloc((nChipsPerLane * 1024 + 1) * sizeof(int), 4096));
    int *lastrowhist = reinterpret_cast<int *>(_mm_malloc((nChipsPerLane * 1024 + 1) * sizeof(int), 4096));

    float *thr_mu  = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * n_row * 1024 * sizeof(float), 4096));
    float *tnoise  = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * n_row * 1024 * sizeof(float), 4096));
    float *S0      = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * 1 * 1024 * sizeof(float), 4096));
    float *S1      = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * 1 * 1024 * sizeof(float), 4096));
    float *S2      = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * 1 * 1024 * sizeof(float), 4096));
    bzero(thr_mu,  nChipsPerLane * n_row * 1024 * sizeof(float));
    bzero(tnoise,  nChipsPerLane * n_row * 1024 * sizeof(float));
    bzero(S0,      nChipsPerLane * 1 * 1024 * sizeof(float));
    bzero(S1,      nChipsPerLane * 1 * 1024 * sizeof(float));
    bzero(S2,      nChipsPerLane * 1 * 1024 * sizeof(float));

    int ngood = 0;
    int nbad = 0;
    int iInj = 0;
    int prev_row = -1;

    if (decoder->isTun)
    {
        reset_stat();
    }
    bzero(rowhist, (nChipsPerLane * 1024 + 1) * sizeof(int));
    std::array<int, 6> chipRow = {{1, 2, 254, 255, 509, 510}};
    while (true)
    {
        EXIT_CODE ret = check_next_event(decoder);
        if (ret == DONE)
        {
            break;
        }
        else if (ret == HB_DATA_DONE)
        {
            nHB_with_data++;
            int iRow = decoder->trigger.thscan_row;
            int iChg = decoder->trigger.thscan_chg;
            if (prev_row != iRow)
            {
                printf("Row %4d : ", iRow);
                fflush(stdout);
                ngood = 0;
                nbad = 0;
                bzero(sumd, nChipsPerLane * 1024 * sizeof(float));
                bzero(sumd2, nChipsPerLane * 1024 * sizeof(float));
                bzero(S0, nChipsPerLane * 1024 * sizeof(float));
                bzero(S1, nChipsPerLane * 1024 * sizeof(float));
                bzero(S2, nChipsPerLane * 1024 * sizeof(float));
                prev_row = iRow;
            }
            // only one strobe per HB for Threshold
            assert(decoder->m_strobeHits.size() == 1);
            nStrobe += decoder->m_strobeHits.size();
            int nhits = decoder->m_hits.size();
            fillrowhist(rowhist, decoder->m_hits.data(), nhits, (decoder->isTun) ? chipRow[iRow] : iRow);
            iInj++;
            if (iInj == decoder->thscan_nInj)
            {
                iInj = 0;
                int nhit = 0;
                for (int i = 0; i < nChipsPerLane * 1024; ++i)
                {
                    nhit += rowhist[i];
                }

                ngood += nhit;
                nbad += rowhist[nChipsPerLane * 1024];
                if (iChg)
                {
                    threshold_next_charge(sumd, sumd2, iChg, lastrowhist, rowhist, decoder->thscan_nInj);

                    tnoise_accumulate(S0, S1, S2, iChg, lastrowhist, rowhist, decoder->thscan_nInj, nChipsPerLane * 1024);
                }
                int *tmp = lastrowhist;
                lastrowhist = rowhist;
                rowhist = tmp;
                bzero(rowhist, (nChipsPerLane * 1024 + 1) * sizeof(int));
                iChg++;
            }
            if (iChg == decoder->thscan_nChg)
            {
                printf("thscan_row %4d ", iRow);
                iChg = 0;
                threshold_next_row(thrs + iRow * nChipsPerLane * 1024, rmss + iRow * nChipsPerLane * 1024, sumd, sumd2, decoder->thscan_nChg, decoder->thscan_nInj);
                tnoise_finalize_row(thr_mu, tnoise, S0, S1, S2, iRow, nChipsPerLane);
                float m, merr, s, serr;
                meanrms(&m, &merr, thrs + decoder->trigger.thscan_row * nChipsPerLane * 1024, nChipsPerLane * 1024);
                meanrms(&s, &serr, rmss + decoder->trigger.thscan_row * nChipsPerLane * 1024, nChipsPerLane * 1024);
                printf(" (mean: %5.2f +/- %4.2f ; RMS: %5.2f +/- %4.2f ; good/bad hits: %d / %d)\n", m, merr, s, serr, ngood, nbad);
                if (decoder->isTun && iRow == 5)
                {
                    ostringstream fname;
                    fname << "./output/" << runnumber << "/" << prefix << ((prefix != "") ? "_" : "") << "thr_map_" << decoder->feeid;
                    fname << "-" << n_vcasn_ithr << ".dat";
                    save_file(fname.str().data(), n_row, nChipsPerLane, thrs);
                    break;

                    {
                        ostringstream f2;
                        f2 << "./output/" << runnumber << "/" << prefix << ((prefix != "") ? "_" : "") << "tnoise_map_" << decoder->feeid << "-" << n_vcasn_ithr << ".dat";
                        save_file(f2.str(), n_row, nChipsPerLane, tnoise);
                    }
                    {
                        ostringstream f3;
                        f3 << "./output/" << runnumber << "/" << prefix << ((prefix != "") ? "_" : "") << "thrmu_map_" << decoder->feeid << "-" << n_vcasn_ithr << ".dat";
                        save_file(f3.str(), n_row, nChipsPerLane, thr_mu);
                    }
                }
            }
        }
        else if (ret == HB_NO_DATA_DONE)
        {
            cout << "Event HB " << nHB << " has no data." << endl;
            continue;
        }
        else
            exit(-1);
    }
    if (!decoder->isTun)
    {
        ostringstream fname;
        fname << "./output/" << runnumber << "/" << prefix << ((prefix != "") ? "_" : "") << "thr_map_" << decoder->feeid << ".dat";
        save_file(fname.str(), n_row, nChipsPerLane, thrs);

        fname.str("");
        fname.clear();
        fname << "./output/" << runnumber << "/" << prefix << ((prefix != "") ? "_" : "") << "rtn_map_" << decoder->feeid << ".dat";
        save_file(fname.str(), n_row, nChipsPerLane, rmss);
    }
    printStat(nHB, nHB_with_data, nStrobe);
    _mm_free(lastrowhist);
    _mm_free(rowhist);
    _mm_free(sumd2);
    _mm_free(sumd);
    _mm_free(rmss);
    _mm_free(thrs);

    _mm_free(S2);
    _mm_free(S1);
    _mm_free(S0);
    _mm_free(tnoise);
    _mm_free(thr_mu);
}

void run_thrs_scan(struct decoder_t *decoder, string &prefix)
{
    std::cout << "Runing THR analysis for feeid " << decoder->feeid << "..." << std::endl;

    constexpr uint8_t nChipsPerLane = 3;
    int16_t n_row = 512;
    float *thrs = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * n_row * 1024 * sizeof(float), 4096));
    float *rmss = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * n_row * 1024 * sizeof(float), 4096));
    float *sumd = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * 1 * 1024 * sizeof(float), 4096));
    float *sumd2 = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * 1 * 1024 * sizeof(float), 4096));
    int *rowhist = reinterpret_cast<int *>(_mm_malloc((nChipsPerLane * 1024 + 1) * sizeof(int), 4096));
    int *lastrowhist = reinterpret_cast<int *>(_mm_malloc((nChipsPerLane * 1024 + 1) * sizeof(int), 4096));

    float *thr_mu  = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * n_row * 1024 * sizeof(float), 4096));
    float *tnoise  = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * n_row * 1024 * sizeof(float), 4096));
    float *S0      = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * 1 * 1024 * sizeof(float), 4096));
    float *S1      = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * 1 * 1024 * sizeof(float), 4096));
    float *S2      = reinterpret_cast<float *>(_mm_malloc(nChipsPerLane * 1 * 1024 * sizeof(float), 4096));
    bzero(thr_mu,  nChipsPerLane * n_row * 1024 * sizeof(float));
    bzero(tnoise,  nChipsPerLane * n_row * 1024 * sizeof(float));
    bzero(S0,      nChipsPerLane * 1 * 1024 * sizeof(float));
    bzero(S1,      nChipsPerLane * 1 * 1024 * sizeof(float));
    bzero(S2,      nChipsPerLane * 1 * 1024 * sizeof(float));

    int ngood = 0;
    int nbad = 0;
    int iInj = 0;
    int prev_row = -1;

    bzero(rowhist, (nChipsPerLane * 1024 + 1) * sizeof(int));
    std::array<int, 6> chipRow = {{1, 2, 254, 255, 509, 510}};
    while (true)
    {
        EXIT_CODE ret = check_next_event(decoder);
        if (ret == DONE)
        {
            break;
        }
        else if (ret == HB_DATA_DONE)
        {
            nHB_with_data++;
            int iRow = decoder->trigger.thscan_row;
            int iChg = decoder->trigger.thscan_chg;
            if (prev_row != iRow)
            {
                fflush(stdout);
                ngood = 0;
                nbad = 0;
                bzero(sumd, nChipsPerLane * 1024 * sizeof(float));
                bzero(sumd2, nChipsPerLane * 1024 * sizeof(float));

                bzero(S0, nChipsPerLane * 1024 * sizeof(float));
                bzero(S1, nChipsPerLane * 1024 * sizeof(float));
                bzero(S2, nChipsPerLane * 1024 * sizeof(float));
                prev_row = iRow;
            }

            // only one strobe per HB for Threshold
            assert(decoder->m_strobeHits.size() == 1);
            nStrobe += decoder->m_strobeHits.size();
            int nhits = decoder->m_hits.size();
            fillrowhist(rowhist, decoder->m_hits.data(), nhits, (decoder->isTun) ? chipRow[iRow] : iRow);
            iInj++;
            if (iInj == decoder->thscan_nInj)
            {
                iInj = 0;
                int nhit = 0;
                for (int i = 0; i < nChipsPerLane * 1024; ++i)
                {
                    nhit += rowhist[i];
                }

                ngood += nhit;
                nbad += rowhist[nChipsPerLane * 1024];
                if (iChg)
                {
                    threshold_next_charge(sumd, sumd2, iChg, lastrowhist, rowhist, decoder->thscan_nInj);
                    tnoise_accumulate(S0, S1, S2, iChg, lastrowhist, rowhist, decoder->thscan_nInj, nChipsPerLane * 1024);
                }
                int *tmp = lastrowhist;
                lastrowhist = rowhist;
                rowhist = tmp;
                bzero(rowhist, (nChipsPerLane * 1024 + 1) * sizeof(int));
                iChg++;
            }
            if (iChg == decoder->thscan_nChg)
            {
                iChg = 0;
                threshold_next_row(thrs + iRow * nChipsPerLane * 1024, rmss + iRow * nChipsPerLane * 1024, sumd, sumd2, decoder->thscan_nChg, decoder->thscan_nInj);
                tnoise_finalize_row(thr_mu, tnoise, S0, S1, S2, iRow, nChipsPerLane);
                float m, merr, s, serr;
                meanrms(&m, &merr, thrs + decoder->trigger.thscan_row * nChipsPerLane * 1024, nChipsPerLane * 1024);
                meanrms(&s, &serr, rmss + decoder->trigger.thscan_row * nChipsPerLane * 1024, nChipsPerLane * 1024);
            }
        }
        else if (ret == HB_NO_DATA_DONE)
        {
            cout << "Event HB " << nHB << " has no data." << endl;
            continue;
        }
        else
            exit(-1);
    }

    ostringstream fname;
    fname << "./output/" << runnumber << "/" << prefix << "_thrs.dat";
    save_file(fname.str(), n_row, nChipsPerLane, thrs);

    ostringstream fnamerms;
    fnamerms << "./output/" << runnumber << "/" << prefix << "_rmss.dat";
    save_file(fnamerms.str(), n_row, nChipsPerLane, rmss);

    {
        ostringstream f2;
        f2 << "./output/" << runnumber << "/" << prefix << "_tnoise.dat";
        save_file(f2.str(), n_row, nChipsPerLane, tnoise);
    }
    {
        ostringstream f3;
        f3 << "./output/" << runnumber << "/" << prefix << "_thrmu.dat";
        save_file(f3.str(), n_row, nChipsPerLane, thr_mu);
    }

    printStat(nHB, nHB_with_data, nStrobe);
    _mm_free(lastrowhist);
    _mm_free(rowhist);
    _mm_free(sumd2);
    _mm_free(sumd);
    _mm_free(rmss);
    _mm_free(thrs);

    _mm_free(S2);
    _mm_free(S1);
    _mm_free(S0);
    _mm_free(tnoise);
    _mm_free(thr_mu);
}

void run_badpix_scan(struct decoder_t *decoder, string &prefix)
{
    std::cout << "Dead pixel scan for feeid " << decoder->feeid << "..." << std::endl;

    uint32_t print_hb_cnt = 10000;
    progressbar bar(100);
    ostringstream ss;

    bar.set_tail_s(ss.str());

    uint32_t *hitmap = reinterpret_cast<uint32_t *>(_mm_malloc(3 * 1024 * 512 * sizeof(uint32_t), 4096));

    bzero(hitmap, 3 * 1024 * 512 * sizeof(uint32_t));

    while (true)
    {
        EXIT_CODE ret = check_next_event(decoder);

        if (!(nHB % print_hb_cnt))
        {
            bar.update();
        }
        if (nHB && (!(nHB % (100 * print_hb_cnt))))
        {
            std::cout << std::endl;
            bar.reset();
        }

        if (ret == DONE)
        {
            break;
        }
        else if (ret == HB_DATA_DONE)
        {
            nHB_with_data++;
            nStrobe += decoder->m_strobeHits.size();
            size_t nhits = decoder->m_hits.size();
            transformhits(decoder->m_hits.data(), nhits);
            fillhitmap(hitmap, decoder->m_hits.data(), nhits);
        }
        else if (ret == HB_NO_DATA_DONE)
        {
            continue;
        }
        else
            exit(-1);
    }

    std::cout << std::endl;

    ostringstream fname;
    fname << "./output/" << runnumber << "/" << prefix << "_deadpix.dat";
    ofstream fhitmap(fname.str().data(), ios_base::trunc);
    fhitmap.write(reinterpret_cast<char *>(hitmap), decoder->nlanes * 1024 * 512 * sizeof(uint32_t));
    fhitmap.close();

    size_t ntot[decoder->nlanes] = {0, 0, 0};
    size_t ntotal = 0;
    for (int lane = 0; lane < decoder->nlanes; ++lane)
    {
        for (int y = 0; y < 512; ++y)
            for (int x = (lane * 1024); x < ((lane + 1) * 1024); ++x)
            {
                ntot[lane] += hitmap[y * decoder->nlanes * 1024 + x];
            }
        ntotal += ntot[lane];
    }
    printf("Total number of hits: %lu. \n", ntotal);
    printStat(nHB, nHB_with_data, nStrobe);
    printTrgCnts(decoder);
}

void get_all_feeids(decoder_t *decoder)
{
    while (true)
    {
        if (check_next_event(decoder) == DONE)
            break;
        else
            exit(-1);
    }

    std::sort(decoder->feeids.begin(), decoder->feeids.end());
    std::cout << decoder->feeids.size() << " feeids founds." << std::endl;
    ostringstream ss;
    ss << "[ ";
    for (auto &_feeid : decoder->feeids)
    {
        ss << _feeid << ((_feeid == decoder->feeids.back()) ? "" : ", ");
    }
    ss << " ]";

    std::cout << ss.str() << std::endl;
}

void signal_callback_handler(int signum)
{
    std::cout << "Signal " << signum << " catched. Exit" << std::endl;
    exit(signum);
}

static void display_help()
{
    printf("Usage: mvtx-thrscan [OPTIONS] <file_name>\n");
    printf("Decode MVTX raw data.\n\n");
    printf("If not file_name given using /dev/stdin by default\n");
    printf("General options:\n");
    printf("  -f <feeid>         Default: 0.\n");
    printf("  -h                 Display help.\n");
    printf("  -i <thr_inj>       number of injection for threshold.\n");
    printf("                     default: 25.\n");
    printf("  -n <n_thr>         number of DAC settings used for threshold tuning.\n");
    printf("  -p <prefix>        prefix text add to the output data file.\n");
    printf("  -t <run_test>      test analysis to run: \n");
    printf("                     0 : run dead pixel scan.\n");
    printf("                     1 : run threshold scan analysis.\n");
    printf("                     2 : run threshold tuning.\n");
}

int main(int argc, char **argv)
{
    enum RUN_TEST
    {
        RUN_DEADPIX,
        RUN_THR,
        RUN_TUNING,
        NO_RUN_TEST
    };
    int opt = -1;
    int _feeid = -1;
    int test = -1;

    int thr_inj = 25;
    int thr_chg = 50;

    int n_thr = 1;
    std::string prefix("");
    std::string filename("/dev/stdin");

    while ((opt = getopt(argc, argv, ":f:hi:n:p:t:")) != -1)
    {
        switch (opt)
        {
        case 'f':
            if (sscanf(optarg, "%d", &_feeid) != 1)
            {
                display_help();
                exit(-1);
            }
            break;

        case 'h':
            display_help();
            exit(0);
            break;

        case 'i':
            if (sscanf(optarg, "%d", &thr_inj) != 1)
            {
                exit(-1);
            }
            break;

        case 'n':
            if (sscanf(optarg, "%d", &n_thr) != 1)
            {
                display_help();
                exit(-1);
            }
            break;

        case 'p':
            prefix = optarg;
            break;

        case 't':
            if (sscanf(optarg, "%d", &test) != 1)
            {
                display_help();
                exit(-1);
            }
            break;

        default:
            display_help();
            exit(-1);
        }
    }

    if (optind != argc) // No file given
    {
        filename = std::string(argv[optind++]);
    }

    if (optind != argc)
    {
        std::cout << "###ERROR!!! Extra arguments given" << std::endl;
        display_help();
        exit(-1);
    }

    // get the run number from prefix, prefix is something like "{Run type}_00067472_0000_FEEID3", the number after the first underscore is the run number
    // split the prefix by underscore and get the second element
    size_t firstUnderscore = prefix.find('_');
    size_t secondUnderscore = prefix.find('_', firstUnderscore + 1);
    std::string numberStr = prefix.substr(firstUnderscore + 1, secondUnderscore - firstUnderscore - 1);
    runnumber = std::stoi(numberStr);
    std::cout << "Processing run number: " << runnumber << std::endl;
    // create output directory if not exists
    std::ostringstream mkdir_cmd;
    mkdir_cmd << "mkdir -p ./output/" << runnumber;
    system(mkdir_cmd.str().c_str());

    signal(SIGINT, signal_callback_handler);

    std::unique_ptr<decoder_t> decoder = std::make_unique<decoder_t>(filename);
    decoder->thscan_nInj = thr_inj;
    decoder->thscan_nChg = thr_chg;
    decoder->feeid = _feeid;
    if (_feeid < 0)
    {
        get_all_feeids(decoder.get());
        exit(0);
    }

    RUN_TEST run_test = NO_RUN_TEST;
    switch (test)
    {
    case 0:
        run_test = RUN_DEADPIX;
        break;

    case 1:
        run_test = RUN_THR;
        break;

    case 2:
        run_test = RUN_TUNING;
        break;

    default:
        run_test = NO_RUN_TEST;
    }

    if (run_test == NO_RUN_TEST)
    {
        std::cout << "### ERROR: no run test provided, doing nothing" << std::endl;
        display_help();
        exit(-1);
    }

    switch (run_test)
    {
    case RUN_DEADPIX:
        decoder.get()->isThr = false;
        run_badpix_scan(decoder.get(), prefix);
        break;
    case RUN_THR:
        decoder.get()->isThr = true;
        decoder.get()->isTun = false;
        run_thrs_scan(decoder.get(), prefix);
        break;
    case RUN_TUNING:
        decoder.get()->isThr = true;
        decoder.get()->isTun = true;
        for (int i = 0; i < n_thr; ++i)
        {
            run_thrs_tuning(decoder.get(), prefix, i);
        }
        run_thrs_tuning(decoder.get(), prefix, -1);
        break;
    default:
        std::cout << "### ERROR: no run test provided, doing nothing" << std::endl;
        exit(-1);
    }

    return 0;
}
