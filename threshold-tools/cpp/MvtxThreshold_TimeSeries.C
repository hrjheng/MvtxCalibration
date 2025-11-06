std::string inputdir = "./results";

#include <map>
#include <string>
#include <vector>

bool verbose = false; // set to true to print debug information
std::string prelimtext = "Internal";

std::map<std::string, std::vector<std::string>> stave_map = {
    {"0", {"L0_03", "L0_04", "L2_01", "L2_02", "L0_05", "L2_03", "L2_04", "L2_15"}}, //
    {"1", {"L1_00", "L1_01", "L2_06", "L2_07", "L1_02", "L1_03", "L2_08", "L2_09"}}, //
    {"2", {"L0_00", "L0_01", "L1_04", "L2_10", "L0_02", "L1_05", "L1_06", "L1_07"}}, //
    {"3", {"L0_09", "L0_10", "L2_05", "L2_11", "L0_11", "L2_12", "L2_13", "L2_14"}}, //
    {"4", {"L0_06", "L0_07", "L1_12", "L2_00", "L0_08", "L1_13", "L1_14", "L1_15"}}, //
    {"5", {"L1_08", "L1_09", "L2_16", "L2_17", "L1_10", "L1_11", "L2_18", "L2_19"}}  //
};

std::vector<std::string> colors_mvtx = {"#77aadd", "#44bb99", "#ddcc77", "#FB9E3A", "#cc3311", "#aa4499"};
std::vector<int> markers = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 40, 48, 47, 46, 45, 44, 43, 41, 39};
std::vector<std::string> colors_layer = {"#143D60", "#89AC46", "#EB5B00"}; // colors for the layers
std::vector<int> markers_layer = {20, 21, 47};                             // markers for the layers

float thrs_lower_bound = 84.5;  // lower bound for the threshold
float thrs_upper_bound = 110.5; // upper bound for the threshold
float rms_lower_bound = 3.5;    // lower bound for the rms
float rms_upper_bound = 10.5;   // upper bound for the rms

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
    if (verbose)
        std::cout << "Read " << thrsmap.size() << " thresholds from " << filename << " done!" << std::endl;
}

void MvtxThreshold_TimeSeries()
{
    TDatime ini_datepoint(2025, 6, 1, 0, 0, 0);
    TGraph *dummy_graph_layer0 = new TGraph();
    TGraph *dummy_graph_layer1 = new TGraph();
    TGraph *dummy_graph_layer2 = new TGraph();

    // initialize a map of TGraph* to store the threshold time series for each stave
    std::map<std::string, TGraph *> stave_thrs_timeseries;
    std::map<std::string, TGraph *> stave_rms_timeseries;
    for (int i = 0; i < 6; ++i)
    {
        std::string mvtx_key = std::to_string(i);
        for (const auto &stave : stave_map[mvtx_key])
        {
            TGraph *graph = new TGraph();
            graph->SetName(stave.c_str());
            stave_thrs_timeseries[stave] = graph;
            TGraph *rms_graph = new TGraph();
            rms_graph->SetName((stave + "_rms").c_str());
            stave_rms_timeseries[stave] = rms_graph;
        }
    }

    std::map<std::string, std::vector<std::vector<float>>> map_layer_thresholds; // to store thresholds for each layer for each time point
    std::map<std::string, std::vector<std::vector<float>>> map_layer_rms;        // to store rms (noise) for each layer for each time point

    // loop over all directories in inputdir
    for (const auto &entry : std::filesystem::directory_iterator(inputdir))
    {
        if (entry.is_directory())
        {
            std::string dir_name = entry.path().filename().string();
            // get the run number and datetime stamp from the directory name
            std::string run_number = dir_name.substr(0, dir_name.find('_'));
            std::string datetime_stamp = dir_name.substr(dir_name.find('_') + 1);
            // get the year, month, day from the datetime stamp
            std::string year = datetime_stamp.substr(0, 4);
            std::string month = datetime_stamp.substr(4, 2);
            std::string day = datetime_stamp.substr(6, 2);
            // print the run number and datetime stamp
            if (verbose)
                std::cout << "Processing directory: " << dir_name << " run number: " << run_number << ", date: " << year << "-" << month << "-" << day << std::endl;

            std::map<std::string, std::pair<double, double>> thrsmap;
            std::map<std::string, std::pair<double, double>> rmsmap;

            TDatime dataPnt(std::stoi(year), std::stoi(month), std::stoi(day), 0, 0, 0);
            dummy_graph_layer0->AddPoint(dataPnt.Convert(), -1);
            dummy_graph_layer1->AddPoint(dataPnt.Convert(), -1);
            dummy_graph_layer2->AddPoint(dataPnt.Convert(), -1);

            if (map_layer_thresholds.find(datetime_stamp) == map_layer_thresholds.end())
            {
                map_layer_thresholds[datetime_stamp] = {{}, {}, {}}; // initialize the map for this datetime_stamp
            }

            if (map_layer_rms.find(datetime_stamp) == map_layer_rms.end())
            {
                map_layer_rms[datetime_stamp] = {{}, {}, {}}; // initialize the map for this datetime_stamp
            }

            // extract the thresholds for all staves from the files in the directory
            for (int mvtx_id = 0; mvtx_id < 6; ++mvtx_id)
            {
                // conditions to skip
                if ((run_number == "67472" && (mvtx_id == 3))                                                                    // skip mvtx_id 3 for run 67472
                    || (run_number == "68756" && (mvtx_id == 1 || mvtx_id == 5))                                                 // skip mvtx_id 1 and 5 for run 68756
                    || (run_number == "68757" && (mvtx_id == 0 || mvtx_id == 2 || mvtx_id == 3 || mvtx_id == 4))                 // skip mvtx_id 0-4 for run 68757
                    || (run_number == "68937" && (mvtx_id == 1))                                                                 // skip mvtx_id 1 for run 68937
                    || (run_number == "68938" && (mvtx_id == 0 || mvtx_id == 2 || mvtx_id == 3 || mvtx_id == 4 || mvtx_id == 5)) // skip mvtx_id 0-5 for run 68938
                    || (run_number == "69111" && (mvtx_id == 1 || mvtx_id == 2 || mvtx_id == 5))                                 // skip mvtx_id 1, 2, and 5 for run 69111
                    || (run_number == "69113" && (mvtx_id == 0 || mvtx_id == 3 || mvtx_id == 4))                                 // skip mvtx_id 0, 3, and 4 for run 69113
                    || (run_number == "69437" && (mvtx_id == 1))                                                                 // skip mvtx_id 1 for run 69437
                    || (run_number == "69438" && (mvtx_id == 0 || mvtx_id == 2 || mvtx_id == 3 || mvtx_id == 4 || mvtx_id == 5)) // skip mvtx_id 0-5 for run 69438
                    || (run_number == "69679" && (mvtx_id == 2 || mvtx_id == 5))                                                 // skip mvtx_id 2 and 5 for run 69679
                    || (run_number == "69680" && (mvtx_id == 0 || mvtx_id == 1 || mvtx_id == 2 || mvtx_id == 3 || mvtx_id == 4)) // skip mvtx_id 0, 1, 2, 3, and 4 for run 69680
                    || (run_number == "69683" && (mvtx_id == 0 || mvtx_id == 1 || mvtx_id == 3 || mvtx_id == 4 || mvtx_id == 5)) // skip mvtx_id 0, 1, 3, 4, and 5 for run 69683
                    || (run_number == "71240" && (mvtx_id == 1))                                                                 // skip mvtx_id 1 for run 71240
                    || (run_number == "71241" && (mvtx_id == 0 || mvtx_id == 2 || mvtx_id == 3 || mvtx_id == 4 || mvtx_id == 5)) // skip mvtx_id 0, 2, 3, 4, and 5 for run 71241
                    || (run_number == "71526")                                                                                   // skip run 71526 fully (since we have run 71528)
                    || (run_number == "71527")                                                                                   // skip run 71527 fully (since we have run 71528)
                    || (run_number == "71790")                                                                                   // skip run 71528 fully (since we have run 71529)
                    || (run_number == "71792")                                                                                   // skip run 71792 fully (since we have run 71793)
                    || (run_number == "71793")                                                                                   // skip run 71793 fully
                    || (run_number == "71877")                                                                                   // skip run 71877 fully
                    || (run_number == "71878")                                                                                   // skip run 71878 fully
                    || (run_number == "71879")                                                                                   // skip run 71879 fully (since we have run 71880)
                    || (run_number == "72011")                                                                                   // skip run 72011 fully
                    || (run_number == "72012")                                                                                   // skip run 72012 fully
                    || (run_number == "72013")                                                                                   // skip run 72013 fully
                    || (run_number == "72014")                                                                                   // skip run 72014 fully
                    || (run_number == "72015")                                                                                   // skip run 72015 fully
                    || (run_number == "72321")                                                                                   // skip run 72321 fully
                    || (run_number == "72324")                                                                                   // skip run 72324 fully
                    || (run_number == "72555")                                                                                   //
                    || (run_number == "72556")                                                                                   //
                    || (run_number == "72558")                                                                                   //
                    || (run_number == "72803")                                                                                   //
                    || (run_number == "72805")                                                                                   //
                    || (run_number == "73072")                                                                                   //
                    || (run_number == "73073")                                                                                   //
                    || (run_number == "73075")                                                                                   //
                    || (run_number == "73489")                                                                                   //
                    || (run_number == "73491")                                                                                   //
                    || (run_number == "74426")                                                                                   //
                    || (run_number == "74428")                                                                                   //
                    || (run_number == "74583")                                                                                   //
                    || (run_number == "74585")                                                                                   //
                    || (run_number == "74649")                                                                                   //
                    || (run_number == "74651")                                                                                   //
                    || (run_number == "75038")                                                                                   //
                    || (run_number == "75039")                                                                                   //
                    || (run_number == "75761")                                                                                   //
                    || (run_number == "75763")                                                                                   //
                    || (run_number == "75853")                                                                                   //
                    || (run_number == "75857")                                                                                   //
                    || (run_number == "76388")                                                                                   //
                    || (run_number == "76390")                                                                                   //
                )
                {
                    continue;
                }

                std::string mvtx_key = std::to_string(mvtx_id);
                std::string mvtx_file = inputdir + "/" + dir_name + "/mvtx" + mvtx_key + "_thresholds.txt";
                readthreshold(mvtx_file, thrsmap);
                std::string mvtx_rms_file = inputdir + "/" + dir_name + "/mvtx" + mvtx_key + "_tnoise.txt";
                readthreshold(mvtx_rms_file, rmsmap);

                for (const auto &pair : stave_thrs_timeseries)
                {
                    const std::string &stave = pair.first;
                    if (thrsmap.find(stave) != thrsmap.end())
                    {
                        double threshold = thrsmap[stave].first;
                        pair.second->AddPoint(dataPnt.Convert(), threshold);
                        if (verbose)
                            std::cout << "Stave: " << stave << ", Threshold: " << threshold << " at " << dataPnt.AsString() << std::endl;

                        for (int layer = 0; layer < 3; ++layer)
                        {
                            if (stave.substr(1, 1) == std::to_string(layer))
                            {
                                // skip values that are below 50 or above 150. They are not valid entries
                                if (threshold > 50 && threshold < 150)
                                {
                                    map_layer_thresholds[datetime_stamp][layer].push_back(threshold);
                                    break;
                                }
                            }
                        }
                    }
                }

                for (const auto &pair : stave_rms_timeseries)
                {
                    const std::string &stave = pair.first;
                    if (rmsmap.find(stave) != rmsmap.end())
                    {
                        // double rms = std::sqrt(rmsmap[stave].first); // convert variance to rms
                        double rms = rmsmap[stave].first;
                        pair.second->AddPoint(dataPnt.Convert(), rms);
                        if (verbose)
                            std::cout << "Stave: " << stave << ", temporal noise: " << rms << " at " << dataPnt.AsString() << std::endl;

                        for (int layer = 0; layer < 3; ++layer)
                        {
                            if (stave.substr(1, 1) == std::to_string(layer))
                            {
                                // skip values that are below 0 or above 20. They are not valid entries
                                if (rms > 0 && rms < 20)
                                {
                                    map_layer_rms[datetime_stamp][layer].push_back(rms);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    float pad_height = 0.32;    // height of each pad
    float size_scale = 2.5;     // scale factor for the size of the pads
    int yaxis_ndivisions = 503; // number of divisions for the y-axis
    TCanvas *c = new TCanvas("c", "c", 1200, 700);
    TPad *pad_layer0 = new TPad("pad_layer0", "Layer 0", 0.0, 0.01, 0.81, 0.01 + pad_height);
    TPad *pad_layer1 = new TPad("pad_layer1", "Layer 1", 0.0, 0.01 + pad_height, 0.81, 0.01 + 2 * pad_height);
    TPad *pad_layer2 = new TPad("pad_layer2", "Layer 2", 0.0, 0.01 + 2 * pad_height, 0.81, 0.01 + 3 * pad_height);
    pad_layer0->SetTopMargin(0.0);
    pad_layer0->SetBottomMargin(0.15);
    pad_layer0->SetRightMargin(0.07);
    pad_layer1->SetTopMargin(0.0);
    pad_layer1->SetBottomMargin(0.0);
    pad_layer1->SetRightMargin(0.07);
    pad_layer2->SetBottomMargin(0.0);
    pad_layer2->SetRightMargin(0.07);
    pad_layer0->Draw();
    pad_layer1->Draw();
    pad_layer2->Draw();
    pad_layer0->cd();
    dummy_graph_layer0->GetXaxis()->SetTimeDisplay(1);
    dummy_graph_layer0->GetXaxis()->SetTimeOffset(0);
    dummy_graph_layer0->GetXaxis()->SetTimeFormat("%Y-%m-%d");
    dummy_graph_layer0->GetXaxis()->SetNdivisions(-4);
    dummy_graph_layer0->GetXaxis()->SetLabelSize(0.035 * size_scale);
    dummy_graph_layer0->GetXaxis()->SetLabelOffset(0.03);
    dummy_graph_layer0->GetYaxis()->SetTitleSize(0.05 * size_scale);
    dummy_graph_layer0->GetYaxis()->SetLabelSize(0.035 * size_scale);
    dummy_graph_layer0->GetYaxis()->SetRangeUser(thrs_lower_bound, thrs_upper_bound);
    dummy_graph_layer0->GetYaxis()->SetNdivisions(yaxis_ndivisions);
    dummy_graph_layer0->Draw("AP");
    TLatex *padlabel_layer0 = new TLatex();
    padlabel_layer0->SetNDC();
    padlabel_layer0->SetTextSize(0.035 * size_scale);
    padlabel_layer0->SetTextAlign(kHAlignLeft + kVAlignCenter);
    padlabel_layer0->DrawLatex(pad_layer0->GetLeftMargin() + 0.04, (1 - pad_layer0->GetTopMargin()) - 0.08, "Layer 0");
    pad_layer1->cd();
    dummy_graph_layer1->GetXaxis()->SetTimeDisplay(1);
    dummy_graph_layer1->GetXaxis()->SetTimeOffset(0);
    dummy_graph_layer1->GetXaxis()->SetNdivisions(-4);
    dummy_graph_layer1->GetYaxis()->SetTitleSize(0.05 * size_scale);
    dummy_graph_layer1->GetYaxis()->SetLabelSize(0.035 * size_scale);
    dummy_graph_layer1->GetYaxis()->SetRangeUser(thrs_lower_bound, thrs_upper_bound);
    dummy_graph_layer1->GetYaxis()->SetNdivisions(yaxis_ndivisions);
    dummy_graph_layer1->Draw("AP");
    TLatex *padlabel_layer1 = new TLatex();
    padlabel_layer1->SetNDC();
    padlabel_layer1->SetTextSize(0.035 * size_scale);
    padlabel_layer1->SetTextAlign(kHAlignLeft + kVAlignCenter);
    padlabel_layer1->DrawLatex(pad_layer1->GetLeftMargin() + 0.04, (1 - pad_layer1->GetTopMargin()) - 0.08, "Layer 1");
    pad_layer2->cd();
    dummy_graph_layer2->GetXaxis()->SetTimeDisplay(1);
    dummy_graph_layer2->GetXaxis()->SetTimeOffset(0);
    dummy_graph_layer2->GetXaxis()->SetNdivisions(-4);
    dummy_graph_layer2->GetYaxis()->SetTitleSize(0.05 * size_scale);
    dummy_graph_layer2->GetYaxis()->SetLabelSize(0.035 * size_scale);
    dummy_graph_layer2->GetYaxis()->SetRangeUser(thrs_lower_bound, thrs_upper_bound);
    dummy_graph_layer2->GetYaxis()->SetNdivisions(yaxis_ndivisions);
    dummy_graph_layer2->Draw("AP");
    TLatex *padlabel_layer2 = new TLatex();
    padlabel_layer2->SetNDC();
    padlabel_layer2->SetTextSize(0.035 * size_scale);
    padlabel_layer2->SetTextAlign(kHAlignLeft + kVAlignCenter);
    padlabel_layer2->DrawLatex(pad_layer2->GetLeftMargin() + 0.04, (1 - pad_layer2->GetTopMargin()) - 0.08, "Layer 2");

    // Loop over the stave_thrs_timeseries map and draw the graphs in the corresponding pads
    for (const auto &pair : stave_thrs_timeseries)
    {
        const std::string &stave = pair.first;
        int layer = std::stoi(stave.substr(1, 1));
        int stave_id = std::stoi(stave.substr(3, 2));
        std::string mvtx_key;
        for (const auto &kv : stave_map)
        {
            if (std::find(kv.second.begin(), kv.second.end(), stave) != kv.second.end())
            {
                mvtx_key = kv.first;
                break;
            }
        }

        if (layer == 0)
        {
            pad_layer0->cd();
        }
        else if (layer == 1)
        {
            pad_layer1->cd();
        }
        else if (layer == 2)
        {
            pad_layer2->cd();
        }
        else
        {
            std::cerr << "Error: Invalid layer " << layer << " for stave " << stave << ". Skipping." << std::endl;
            continue;
        }

        pair.second->SetMarkerStyle(markers[stave_id]);
        pair.second->SetMarkerSize(0.75);
        pair.second->SetMarkerColor(TColor::GetColor(colors_mvtx[std::stoi(mvtx_key)].c_str()));
        pair.second->SetLineColor(TColor::GetColor(colors_mvtx[std::stoi(mvtx_key)].c_str()));
        pair.second->Draw("P same");
    }
    c->RedrawAxis();
    c->cd();
    // Manually add the Y-axis title
    TLatex *ytitle = new TLatex();
    ytitle->SetNDC();
    ytitle->SetTextSize(0.04);
    ytitle->SetTextAlign(kHAlignCenter + kVAlignCenter);
    ytitle->SetTextAngle(90);
    ytitle->DrawLatex(0.05, 0.5, "Average chip threshold [e^{-}]");
    // Add a big legend at the right side
    float entry_height = 0.03;
    float legend_textsize = 0.02;
    TLegend *leg_l0 = new TLegend(0.81, 0.95 - 12 * entry_height, 0.87, 0.95);
    leg_l0->SetNColumns(1);
    leg_l0->SetTextSize(legend_textsize);
    leg_l0->SetFillColor(0);
    leg_l0->SetBorderSize(0);
    for (const auto &pair : stave_thrs_timeseries)
    {
        if (pair.first.substr(1, 1) != "0")
            continue; // Only add layer 0 staves to the legend
        leg_l0->AddEntry(pair.second, pair.first.c_str(), "p");
    }
    leg_l0->Draw();
    TLegend *leg_l1 = new TLegend(0.87, 0.95 - 16 * entry_height, 0.93, 0.95);
    leg_l1->SetNColumns(1);
    leg_l1->SetTextSize(legend_textsize);
    leg_l1->SetFillColor(0);
    leg_l1->SetBorderSize(0);
    for (const auto &pair : stave_thrs_timeseries)
    {
        if (pair.first.substr(1, 1) != "1")
            continue; // Only add layer 1 staves to the legend
        leg_l1->AddEntry(pair.second, pair.first.c_str(), "p");
    }
    leg_l1->Draw();
    TLegend *leg_l2 = new TLegend(0.93, 0.95 - 20 * entry_height, 0.99, 0.95);
    leg_l2->SetNColumns(1);
    leg_l2->SetTextSize(legend_textsize);
    leg_l2->SetFillColor(0);
    leg_l2->SetBorderSize(0);
    for (const auto &pair : stave_thrs_timeseries)
    {
        if (pair.first.substr(1, 1) != "2")
            continue; // Only add layer 2 staves to the legend
        leg_l2->AddEntry(pair.second, pair.first.c_str(), "p");
    }
    leg_l2->Draw();
    c->SaveAs(Form("%s/MvtxThreshold_TimeSeries.pdf", inputdir.c_str()));
    c->SaveAs(Form("%s/MvtxThreshold_TimeSeries.png", inputdir.c_str()));

    // Similar plot for RMS (temporal noise)
    TCanvas *c_rms_allstave = new TCanvas("c_rms_allstave", "RMS All Stave", 1200, 700);
    pad_layer0->Draw();
    pad_layer1->Draw();
    pad_layer2->Draw();
    pad_layer0->cd();
    dummy_graph_layer0->GetYaxis()->SetRangeUser(rms_lower_bound, rms_upper_bound);
    dummy_graph_layer0->Draw("AP");
    padlabel_layer0->DrawLatex(pad_layer0->GetLeftMargin() + 0.04, (1 - pad_layer0->GetTopMargin()) - 0.08, "Layer 0");
    pad_layer1->cd();
    dummy_graph_layer1->GetYaxis()->SetRangeUser(rms_lower_bound, rms_upper_bound);
    dummy_graph_layer1->Draw("AP");
    padlabel_layer1->DrawLatex(pad_layer1->GetLeftMargin() + 0.04, (1 - pad_layer1->GetTopMargin()) - 0.08, "Layer 1");
    pad_layer2->cd();
    dummy_graph_layer2->GetYaxis()->SetRangeUser(rms_lower_bound, rms_upper_bound);
    dummy_graph_layer2->Draw("AP");
    padlabel_layer2->DrawLatex(pad_layer2->GetLeftMargin() + 0.04, (1 - pad_layer2->GetTopMargin()) - 0.08, "Layer 2");
    // Loop over the stave_rms_timeseries map and draw the graphs in the corresponding pads
    for (const auto &pair : stave_rms_timeseries)
    {
        const std::string &stave = pair.first;
        int layer = std::stoi(stave.substr(1, 1));
        int stave_id = std::stoi(stave.substr(3, 2));
        std::string mvtx_key;
        for (const auto &kv : stave_map)
        {
            if (std::find(kv.second.begin(), kv.second.end(), stave) != kv.second.end())
            {
                mvtx_key = kv.first;
                break;
            }
        }

        if (layer == 0)
        {
            pad_layer0->cd();
        }
        else if (layer == 1)
        {
            pad_layer1->cd();
        }
        else if (layer == 2)
        {
            pad_layer2->cd();
        }
        else
        {
            std::cerr << "Error: Invalid layer " << layer << " for stave " << stave << ". Skipping." << std::endl;
            continue;
        }

        pair.second->SetMarkerStyle(markers[stave_id]);
        pair.second->SetMarkerSize(0.75);
        pair.second->SetMarkerColor(TColor::GetColor(colors_mvtx[std::stoi(mvtx_key)].c_str()));
        pair.second->SetLineColor(TColor::GetColor(colors_mvtx[std::stoi(mvtx_key)].c_str()));
        pair.second->Draw("P same");
    }

    c_rms_allstave->RedrawAxis();
    c_rms_allstave->cd();
    ytitle->SetNDC();
    ytitle->SetTextSize(0.04);
    ytitle->SetTextAlign(kHAlignCenter + kVAlignCenter);
    ytitle->SetTextAngle(90);
    ytitle->DrawLatex(0.05, 0.5, "Average chip temporal noise (RMS) [e^{-}]");
    leg_l0->Draw();
    leg_l1->Draw();
    leg_l2->Draw();
    c_rms_allstave->SaveAs(Form("%s/MvtxRMS_TimeSeries.pdf", inputdir.c_str()));
    c_rms_allstave->SaveAs(Form("%s/MvtxRMS_TimeSeries.png", inputdir.c_str()));

    // For the average threshold graphs, draw them in a new canvas
    // Calculate the average threshold for each layer for each time point
    // TGraphErrors for the layer average threshold
    std::vector<TGraphErrors *> grapherr_layer_avg;
    std::vector<TGraphErrors *> grapherr_layer_rms_avg;
    for (int i = 0; i < 3; ++i)
    {
        TGraphErrors *grapherr = new TGraphErrors();
        grapherr->SetName(Form("grapherr_layer%d_avg", i));
        grapherr_layer_avg.push_back(grapherr);
        TGraphErrors *grapherr_rms = new TGraphErrors();
        grapherr_rms->SetName(Form("grapherr_layer%d_rms_avg", i));
        grapherr_layer_rms_avg.push_back(grapherr_rms);
    }
    std::vector<TGraph *> graph_layer_avg;
    std::vector<TGraph *> graph_layer_rms_avg;
    for (int i = 0; i < 3; ++i)
    {
        TGraph *graph = new TGraph();
        graph->SetName(Form("graph_layer%d_avg", i));
        graph_layer_avg.push_back(graph);
        TGraph *graph_rms = new TGraph();
        graph_rms->SetName(Form("graph_layer%d_rms_avg", i));
        graph_layer_rms_avg.push_back(graph_rms);
    }

    for (const auto &pair : map_layer_thresholds)
    {
        const std::string &datetime_stamp = pair.first;
        int year = std::stoi(datetime_stamp.substr(0, 4));
        int month = std::stoi(datetime_stamp.substr(4, 2));
        int day = std::stoi(datetime_stamp.substr(6, 2));
        TDatime dataPnt(year, month, day, 0, 0, 0);
        const std::vector<std::vector<float>> &layer_thresholds = pair.second;

        for (int layer = 0; layer < 3; ++layer)
        {
            if (layer_thresholds[layer].empty())
                continue;
            double avg = std::accumulate(layer_thresholds[layer].begin(), layer_thresholds[layer].end(), 0.0) / layer_thresholds[layer].size();
            double sq_sum = std::accumulate(layer_thresholds[layer].begin(), layer_thresholds[layer].end(), 0.0, [avg](double acc, double x) { return acc + (x - avg) * (x - avg); });
            double stdev = std::sqrt(sq_sum / layer_thresholds[layer].size());
            grapherr_layer_avg[layer]->AddPoint(dataPnt.Convert(), avg);
            grapherr_layer_avg[layer]->SetPointError(grapherr_layer_avg[layer]->GetN() - 1, 0.0, stdev);
            graph_layer_avg[layer]->AddPoint(dataPnt.Convert(), avg);
        }
    }

    for (const auto &pair : map_layer_rms)
    {
        const std::string &datetime_stamp = pair.first;
        int year = std::stoi(datetime_stamp.substr(0, 4));
        int month = std::stoi(datetime_stamp.substr(4, 2));
        int day = std::stoi(datetime_stamp.substr(6, 2));
        TDatime dataPnt(year, month, day, 0, 0, 0);
        const std::vector<std::vector<float>> &layer_rms = pair.second;
        for (int layer = 0; layer < 3; ++layer)
        {
            if (layer_rms[layer].empty())
                continue;
            double avg = std::accumulate(layer_rms[layer].begin(), layer_rms[layer].end(), 0.0) / layer_rms[layer].size();
            double sq_sum = std::accumulate(layer_rms[layer].begin(), layer_rms[layer].end(), 0.0, [avg](double acc, double x) { return acc + (x - avg) * (x - avg); });
            double stdev = std::sqrt(sq_sum / layer_rms[layer].size());
            grapherr_layer_rms_avg[layer]->AddPoint(dataPnt.Convert(), avg);
            grapherr_layer_rms_avg[layer]->SetPointError(grapherr_layer_rms_avg[layer]->GetN() - 1, 0.0, stdev);
            graph_layer_rms_avg[layer]->AddPoint(dataPnt.Convert(), avg);
        }
    }

    gStyle->SetEndErrorSize(2);
    TCanvas *c_avg = new TCanvas("c_avg", "c_avg", 900, 600);
    c_avg->cd();
    gPad->SetRightMargin(0.07);
    dummy_graph_layer0->GetXaxis()->SetTimeDisplay(1);
    dummy_graph_layer0->GetXaxis()->SetTimeOffset(0);
    dummy_graph_layer0->GetXaxis()->SetTimeFormat("%Y-%m-%d");
    dummy_graph_layer0->GetXaxis()->SetNdivisions(-4);
    dummy_graph_layer0->GetXaxis()->SetLabelSize(0.04);
    dummy_graph_layer0->GetXaxis()->SetLabelOffset(0.03);
    dummy_graph_layer0->GetYaxis()->SetTitleSize(0.05);
    dummy_graph_layer0->GetYaxis()->SetLabelSize(0.05);
    dummy_graph_layer0->GetYaxis()->SetTitle("Average chip threshold [e^{-}]");
    dummy_graph_layer0->GetYaxis()->SetRangeUser(thrs_lower_bound, thrs_upper_bound);
    dummy_graph_layer0->GetYaxis()->SetNdivisions(yaxis_ndivisions);
    dummy_graph_layer0->Draw("AP");
    for (int i = 0; i < 3; ++i)
    {
        grapherr_layer_avg[i]->SetMarkerStyle(markers_layer[i]);
        grapherr_layer_avg[i]->SetMarkerSize(1.5);
        grapherr_layer_avg[i]->SetMarkerColor(TColor::GetColor(colors_layer[i].c_str()));
        grapherr_layer_avg[i]->SetLineColor(TColor::GetColor(colors_layer[i].c_str()));
        grapherr_layer_avg[i]->SetLineWidth(2);
        grapherr_layer_avg[i]->Draw("P same");
        graph_layer_avg[i]->SetLineStyle(2);
        graph_layer_avg[i]->SetLineWidth(1);
        graph_layer_avg[i]->SetLineColor(TColor::GetColor(colors_layer[i].c_str()));
        graph_layer_avg[i]->SetMarkerSize(0);
        graph_layer_avg[i]->Draw("L same");
    }
    c_avg->RedrawAxis();
    // Add a legend for the average thresholds
    TLegend *leg_avg = new TLegend(gPad->GetLeftMargin() + 0.05,      //
                                   (1 - gPad->GetTopMargin()) - 0.25, //
                                   gPad->GetLeftMargin() + 0.25,      //
                                   (1 - gPad->GetTopMargin()) - 0.06  //
    );
    leg_avg->SetTextSize(0.045);
    leg_avg->SetFillColor(0);
    leg_avg->SetBorderSize(0);
    for (int i = 0; i < 3; ++i)
    {
        leg_avg->AddEntry(grapherr_layer_avg[i], Form("Layer %d", i), "PEL");
    }
    leg_avg->Draw();
    TLegend *sphnxleg = new TLegend(1 - gPad->GetRightMargin() - 0.33, //
                                    (1 - gPad->GetTopMargin()) - 0.18, //
                                    1 - gPad->GetRightMargin() - 0.13, //
                                    (1 - gPad->GetTopMargin()) - 0.06  //
    );
    sphnxleg->SetTextAlign(kHAlignLeft + kVAlignCenter);
    sphnxleg->SetTextSize(0.04);
    sphnxleg->SetFillStyle(0);
    sphnxleg->AddEntry("", Form("#it{#bf{sPHENIX}} %s", prelimtext.c_str()), "");
    // sphnxleg->AddEntry("", "Au+Au #sqrt{s_{NN}}=200 GeV", "");
    sphnxleg->AddEntry("", "Run2025", "");
    sphnxleg->Draw();
    c_avg->SaveAs(Form("%s/MvtxThreshold_TimeSeries_LayerAverage.pdf", inputdir.c_str()));
    c_avg->SaveAs(Form("%s/MvtxThreshold_TimeSeries_LayerAverage.png", inputdir.c_str()));

    TCanvas *c_rms = new TCanvas("c_rms", "c_rms", 900, 600);
    c_rms->cd();
    gPad->SetRightMargin(0.07);
    dummy_graph_layer0->GetXaxis()->SetTimeDisplay(1);
    dummy_graph_layer0->GetXaxis()->SetTimeOffset(0);
    dummy_graph_layer0->GetXaxis()->SetTimeFormat("%Y-%m-%d");
    dummy_graph_layer0->GetXaxis()->SetNdivisions(-4);
    dummy_graph_layer0->GetXaxis()->SetLabelSize(0.04);
    dummy_graph_layer0->GetXaxis()->SetLabelOffset(0.03);
    dummy_graph_layer0->GetYaxis()->SetTitleSize(0.05);
    dummy_graph_layer0->GetYaxis()->SetLabelSize(0.05);
    dummy_graph_layer0->GetYaxis()->SetTitle("Average chip temporal noise (RMS) [e^{-}]");
    dummy_graph_layer0->GetYaxis()->SetRangeUser(rms_lower_bound, rms_upper_bound);
    dummy_graph_layer0->GetYaxis()->SetNdivisions(yaxis_ndivisions);
    dummy_graph_layer0->Draw("AP");
    for (int i = 0; i < 3; ++i)
    {
        grapherr_layer_rms_avg[i]->SetMarkerStyle(markers_layer[i]);
        grapherr_layer_rms_avg[i]->SetMarkerSize(1.5);
        grapherr_layer_rms_avg[i]->SetMarkerColor(TColor::GetColor(colors_layer[i].c_str()));
        grapherr_layer_rms_avg[i]->SetLineColor(TColor::GetColor(colors_layer[i].c_str()));
        grapherr_layer_rms_avg[i]->SetLineWidth(2);
        grapherr_layer_rms_avg[i]->Draw("P same");
        graph_layer_rms_avg[i]->SetLineStyle(2);
        graph_layer_rms_avg[i]->SetLineWidth(1);
        graph_layer_rms_avg[i]->SetLineColor(TColor::GetColor(colors_layer[i].c_str()));
        graph_layer_rms_avg[i]->SetMarkerSize(0);
        graph_layer_rms_avg[i]->Draw("L same");
    }
    c_rms->RedrawAxis();
    // can use the same legend as above
    leg_avg->Draw();
    sphnxleg->Draw();
    c_rms->SaveAs(Form("%s/MvtxRMS_TimeSeries_LayerAverage.pdf", inputdir.c_str()));
    c_rms->SaveAs(Form("%s/MvtxRMS_TimeSeries_LayerAverage.png", inputdir.c_str()));
}
