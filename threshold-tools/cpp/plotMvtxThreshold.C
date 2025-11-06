#include "MVTX_OnlMonPoly.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <fstream>
#include <sstream>

int layer_starting_index[3] = {0, 12, 28};

void Pal_blue()
{
    static Int_t colors[1000];
    static Bool_t initialized = kFALSE;
    Double_t Red[5] = {250. / 255., 178. / 255., 141. / 255., 104. / 255., 67. / 255.};
    Double_t Green[5] = {250. / 255., 201. / 255., 174. / 255., 146. / 255., 120. / 255.};
    Double_t Blue[5] = {250. / 255., 230. / 255., 216. / 255., 202. / 255., 188. / 255.};
    Double_t Length[5] = {0.00, 0.25, 0.50, 0.75, 1.00};
    if (!initialized)
    {
        Int_t FI = TColor::CreateGradientColorTable(5, Length, Red, Green, Blue, 1000);
        for (int i = 0; i < 1000; i++)
            colors[i] = FI + i;
        initialized = kTRUE;
        return;
    }
    gStyle->SetPalette(1000, colors);
}

void Pal_red()
{
    static Int_t colors[1000];
    static Bool_t initialized = kFALSE;
    Double_t Red[5] = {255. / 255., 255. / 255., 255. / 255., 205. / 255., 153. / 255.};
    Double_t Green[5] = {248. / 255., 185. / 255., 133. / 255., 81. / 255., 29. / 255.};
    Double_t Blue[5] = {248. / 255., 185. / 255., 133. / 255., 81. / 255., 29. / 255.};
    Double_t Length[5] = {0.00, 0.25, 0.50, 0.75, 1.00};
    if (!initialized)
    {
        Int_t FI = TColor::CreateGradientColorTable(5, Length, Red, Green, Blue, 1000);
        for (int i = 0; i < 1000; i++)
            colors[i] = FI + i;
        initialized = kTRUE;
        return;
    }
    gStyle->SetPalette(1000, colors);
}

void ReverseXAxis(TH1 *h)
{
    // Remove the current axis
    h->GetXaxis()->SetLabelOffset(999);
    //  h->GetXaxis()->SetTickLength(0);

    // Redraw the new axis
    gPad->Update();
    TGaxis *newaxis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmin(), gPad->GetUymin(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), 510, "-");
    newaxis->SetLabelOffset(-0.03);
    newaxis->SetTitleFont(42);
    newaxis->SetLabelFont(42);
    newaxis->SetTitleSize(0.05);
    newaxis->SetLabelSize(0.05);
    newaxis->Draw();
}

void readthreshold(std::string filename, std::map<std::string, std::pair<double, double>> &thrsmap)
{
    std::ifstream infile(filename);
    std::string line;

    // Regular expression to match the line format
    std::regex pattern(R"((\w+):\s*([\d.]+)\s+([\d.]+))");

    if (!infile.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    while (std::getline(infile, line))
    {
        // remove "$", "e", "^", "-", and "\pm" from the line
        line.erase(std::remove(line.begin(), line.end(), 'e'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '$'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '^'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '-'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\\'), line.end());
        line.erase(std::remove(line.begin(), line.end(), 'p'), line.end());
        line.erase(std::remove(line.begin(), line.end(), 'm'), line.end());
        std::smatch match;
        if (std::regex_search(line, match, pattern))
        {
            std::string stave = match[1];
            double threshold = std::stod(match[2]);
            double stddev = std::stod(match[3]);

            thrsmap[stave] = std::make_pair(threshold, stddev);
        }
        else
        {
            std::cerr << "No match found in line: " << line << std::endl;
        }
    }
    infile.close();
    std::cout << "Read " << thrsmap.size() << " thresholds from " << filename << " done!" << std::endl;
}

void plotMvtxThreshold(std::string inputfiledir)
{
    // get the runnumber and date from the input file directory
    std::string runnumber = inputfiledir.substr(inputfiledir.find_last_of('/') + 1);
    runnumber = runnumber.substr(0, runnumber.find('_'));
    std::string date = inputfiledir.substr(inputfiledir.find_last_of('/') + 1);
    date = date.substr(date.find('_') + 1);

    // format date to YYYY-MM-DD
    if (date.length() >= 8)
    {
        date = date.substr(0, 4) + "-" + date.substr(4, 2) + "-" + date.substr(6, 2);
    }
    else
    {
        std::cerr << "Warning: date string '" << date << "' is too short to format as YYYY-MM-DD." << std::endl;
    }
    std::cout << "Run number: " << runnumber << ", Date: " << date << std::endl;

    std::map<std::string, std::pair<double, double>> thrsmap;
    std::map<std::string, std::pair<double, double>> rmsmap;

    // histogram for the thresholds
    TH2Poly *h = new TH2Poly();
    createPoly(h);
    h->SetMinimum(80.);
    h->SetMaximum(110.);

    // histogram for the rms (noise)
    TH2Poly *h_rms = new TH2Poly();
    createPoly(h_rms);
    h_rms->SetMinimum(3.5);
    h_rms->SetMaximum(10.5);

    for (int flx = 0; flx < 6; flx++)
    {
        // read the thresholds from the file
        std::string fname = inputfiledir + "/mvtx" + std::to_string(flx) + "_thresholds.txt";
        readthreshold(fname, thrsmap);
        // read the rmss from the file
        std::string fname_rms = inputfiledir + "/mvtx" + std::to_string(flx) + "_tnoise.txt";
        readthreshold(fname_rms, rmsmap);
    }

    // print out the map
    for (const auto &pair : thrsmap)
    {
        std::cout << pair.first << ": " << pair.second.first << " " << pair.second.second << std::endl;
        std::string stave = pair.first;
        double threshold = pair.second.first;
        double stddev = pair.second.second;
        int layer = std::stoi(stave.substr(1, 1));
        int staveid = std::stoi(stave.substr(3, 2));
        int index = layer_starting_index[layer] + staveid;

        double *px = new double[4];
        double *py = new double[4];
        getStavePoint(layer, staveid, px, py);
        for (int i = 0; i < 4; i++)
        {
            px[i] /= 10.;
            py[i] /= 10.;
        }
        double avg_x = 0.;
        double avg_y = 0.;
        for (int i = 0; i < 4; i++)
        {
            avg_x += px[i] / 4;
            avg_y += py[i] / 4;
        }
        h->SetBinContent(mapstave[layer][staveid], threshold);
        h->SetBinError(mapstave[layer][staveid], stddev);
        delete[] px;
        delete[] py;
    }

    for (const auto &pair : rmsmap)
    {
        std::cout << "rmss - " << pair.first << ": " << pair.second.first << " " << pair.second.second << std::endl;
        std::string stave = pair.first;
        double rmss = pair.second.first;
        double stddev = pair.second.second;
        int layer = std::stoi(stave.substr(1, 1));
        int staveid = std::stoi(stave.substr(3, 2));
        int index = layer_starting_index[layer] + staveid;
        double *px = new double[4];
        double *py = new double[4];
        getStavePoint(layer, staveid, px, py);
        for (int i = 0; i < 4; i++)
        {
            px[i] /= 10.;
            py[i] /= 10.;
        }
        double avg_x = 0.;
        double avg_y = 0.;
        for (int i = 0; i < 4; i++)
        {
            avg_x += px[i] / 4;
            avg_y += py[i] / 4;
        }
        // h_rms->SetBinContent(mapstave[layer][staveid], std::sqrt(rmss));
        h_rms->SetBinContent(mapstave[layer][staveid], rmss);
        h_rms->SetBinError(mapstave[layer][staveid], stddev);
        delete[] px;
        delete[] py;
    }

    TCanvas *c = new TCanvas("c", "c", 830, 700);
    gStyle->SetPalette(kLightTemperature);
    // TColor::InvertPalette();
    gStyle->SetPaintTextFormat(".2f");
    gPad->SetLeftMargin(0.1);
    gPad->SetBottomMargin(0.1);
    gPad->SetRightMargin(0.22);
    gPad->SetTopMargin(0.07);
    h->GetZaxis()->SetTitle("Threshold [e^{-}]");
    h->GetZaxis()->SetTitleOffset(1.7);
    h->Draw("L COLZ TEXT");
    // static TExec *ex1 = new TExec("ex1", "Pal_red();");
    // ex1->Draw();
    // h->Draw("L colz text same");
    ReverseXAxis(h);
    // Print the run number and date
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.05);
    latex->SetTextAlign(kHAlignLeft + kVAlignCenter);
    latex->DrawLatex(gPad->GetLeftMargin(), (1 - gPad->GetTopMargin()) + 0.025, Form("Run %s, %s", runnumber.c_str(), date.c_str()));
    c->SaveAs((inputfiledir + "/thresholds_staves.png").c_str());
    c->SaveAs((inputfiledir + "/thresholds_staves.pdf").c_str());

    c->Clear();
    gPad->SetLeftMargin(0.1);
    gPad->SetBottomMargin(0.1);
    gPad->SetRightMargin(0.22);
    gPad->SetTopMargin(0.07);
    h_rms->GetZaxis()->SetTitle("Average chip temporal noise [e^{-}]");
    h_rms->GetZaxis()->SetTitleOffset(1.7);
    h_rms->Draw("L COLZ TEXT");
    ReverseXAxis(h_rms);
    latex->DrawLatex(gPad->GetLeftMargin(), (1 - gPad->GetTopMargin()) + 0.025, Form("Run %s, %s", runnumber.c_str(), date.c_str()));
    c->SaveAs((inputfiledir + "/rms_staves.png").c_str());
    c->SaveAs((inputfiledir + "/rms_staves.pdf").c_str());
}
