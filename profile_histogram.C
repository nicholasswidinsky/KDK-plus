#include <fstream>
#include <sstream>
#include <iostream>
#include<string>
#include <filesystem>
#include "TCanvas.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCutG.h"
#include <iomanip>   // for std::setprecision
#include <fstream>   // for std::ofstream
#include "TMath.h"   // for TMath::Prob
#include "TSystem.h" // for gSystem->mkdir

// void fit2D_cut() {
void profile_histogram() {

    // TCanvas *c = new TCanvas("c","2D Histogram with Banana Cut",1000,800);
    // c->Divide(1,1); 

    std::vector<TH2D*> histograms;

    


    //This Tries to open the file. If it fails it returns an error. 
    // std::string inputPath("/home/nick/PhD/KDK+/Annulus_Compton_scatter_V1/2026_01_28/NaI_annulus_LS_2_Cs137_All_NaI_higher_HV/RAW/coinc_sorted/SDataR_NaI_annulus_LS_2_Cs137_All_NaI_higher_HV_coinc_4_5_8.txt");
    std::string inputPath("/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/2026_04_28/2026_04_28_Daily_LSC_calibration_Cs137_coinc_2/RAW/coinc_sorted/SDataR_2026_04_28_Daily_LSC_calibration_Cs137_coinc_2_coinc_4_5_8.txt");
    
    std::filesystem::path filePath(inputPath);

    std::string inputDir      = filePath.parent_path().string();
    std::string inputBasename = filePath.stem().string();
    std::string inputFilename = filePath.filename().string();

    std::string outputDir = inputDir + "/" + inputBasename + "_fit_results";
    gSystem->mkdir(outputDir.c_str(), kTRUE);  // kTRUE creates intermediate dirs

    // Open file using the same variable
    std::ifstream file(inputPath);


    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    std::string line; //Initialize the line variable for when we read in lines. 

    // Skip first 3 lines
    for (int i = 0; i < 3; i++) {
        std::getline(file, line);
    }

    // Read the 4th line to extract n
    std::getline(file, line);
    int n = 0;
    size_t colonPos = line.rfind(':');
    if (colonPos != std::string::npos) {
        std::istringstream(line.substr(colonPos + 1)) >> n; //Read in the number of channels from the 4th line and save it as n.
    } else {
        std::cerr << "Error: could not parse number of channels!" << std::endl;
        return;
    }
    std::cout << "Number of channels: " << n << std::endl;

    // Read first data line to extract channel IDs
    std::string firstRowData;
    std::vector<int> channelIDs(n, -1);
    if (std::getline(file, line) && !line.empty()) {
        std::stringstream ss(line);
        std::string token;
        std::vector<double> firstRow;
        while (std::getline(ss, token, ';')) {
            try { firstRow.push_back(std::stod(token)); }
            catch (...) { firstRow.push_back(0.0); }
        }
        for (int i = 0; i < n; i++) {
            int idIdx = 6 * i + 1;
            if (idIdx < (int)firstRow.size()) {
                channelIDs[i] = (int)firstRow[idIdx];
            }
        }
        // Store first row so it isn't skipped in the main loop
        firstRowData = line;
    }

    // Now create histograms using real channel IDs
    // std::vector<TH2D*> histograms;
    std::vector<std::pair<int,int>> pairs;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            std::string chI = std::to_string(channelIDs[i]);
            std::string chJ = std::to_string(channelIDs[j]);

            std::string name  = "h_ch" + chI + "_vs_ch" + chJ;
            std::string title = "Channel " + chI + " vs Channel " + chJ +
                                ";Channel " + chI + " Energy" +
                                ";Channel " + chJ + " Energy";
            histograms.push_back(new TH2D(name.c_str(), title.c_str(), 100, 0, 4000, 100, 0, 4000));
            pairs.push_back({i, j});
        }
    }

    // Main data reading loop - prepend the first row back in
    std::istringstream firstStream(firstRowData);
    auto processRow = [&](const std::string& rowLine) {
        std::stringstream ss(rowLine);
        std::string token;
        std::vector<double> rowValues;
        while (std::getline(ss, token, ';')) {
            try { rowValues.push_back(std::stod(token)); }
            catch (...) { rowValues.push_back(0.0); }
        }
        std::vector<double> channelVals(n, 0.0);
        for (int i = 0; i < n; i++) {
            int dataIdx = 6 * i + 3;
            if (dataIdx < (int)rowValues.size()) {
                channelVals[i] = rowValues[dataIdx];
            }
        }
        for (int k = 0; k < (int)pairs.size(); k++) {
            histograms[k]->Fill(channelVals[pairs[k].first], channelVals[pairs[k].second]);
        }
    };

    // Process first row then continue with rest of file
    processRow(firstRowData);
    while (std::getline(file, line)) {
        if (!line.empty()) processRow(line);
    }


    file.close();


for (int k = 0; k < (int)pairs.size(); k++) {
    int i = pairs[k].first;
    int j = pairs[k].second;

    std::string chI = std::to_string(channelIDs[i]);
    std::string chJ = std::to_string(channelIDs[j]);

    std::string suffix      = "ch" + chI + "_vs_ch" + chJ;
    std::string canvasName  = "c_"   + suffix;
    std::string canvasTitle = "Channel " + chI + " vs Channel " + chJ;
    std::string cutName     = "banana_" + suffix;
    std::string proName     = "pro_"    + suffix;
    std::string fitName     = "fit_"    + suffix;

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 800, 600);
    c->cd(1);
    histograms[k]->Draw("COLZ");

    TCutG *cut = new TCutG(cutName.c_str(), 4);
    cut->SetPoint(0, 2500,    0);
    cut->SetPoint(1,    0,    0);
    cut->SetPoint(2,    0, 1000);
    cut->SetPoint(3, 2500,    0);
    cut->SetLineColor(kRed+2);
    cut->SetLineWidth(3);
    cut->Draw("same");

    TProfile *pro = new TProfile(proName.c_str(),
        "Profile Histogram;NaI scintillator;Liquid scintillator",
        100, 0, 4000, 0, 4000, "");

    for (int bi = 1; bi <= histograms[k]->GetNbinsX(); bi++) {
        for (int bj = 1; bj <= histograms[k]->GetNbinsY(); bj++) {
            double x = histograms[k]->GetXaxis()->GetBinCenter(bi);
            double y = histograms[k]->GetYaxis()->GetBinCenter(bj);
            double z = histograms[k]->GetBinContent(bi, bj);
            if (!cut->IsInside(x, y) && z > 0) {
                pro->Fill(x, y, z);
            }
        }
    }
    pro->Draw("same");

    TF1 *fit = new TF1(fitName.c_str(), "pol1", 300, 1500);
    pro->Fit(fit, "R");
    fit->SetLineColor(kRed);
    fit->Draw("same");    
    
    std::string outputDir = inputDir + "/" + inputBasename + "_fit_results";
    gSystem->mkdir(outputDir.c_str(), kTRUE);  // kTRUE creates intermediate dirs



    // // ofstream out(seglist[0] + "_fit_results/Fit_results");
    // for (int i = 0; i < fit->GetNpar(); i++){
    //     out << fit->GetParName(i) << " " << fit->GetParameter(i)<< " " << fit->GetParError(i) << endl;
    // }
    // out.close();

    Double_t chi2 = fit->GetChisquare();
    Double_t p0 = fit->GetParameter(0);
    Double_t p1 = fit->GetParameter(1);
    Double_t e0 = fit->GetParError(0);
    Double_t e1 = fit->GetParError(1);



    c->Update();
    std::string pdfPath = outputDir + "/" + inputBasename + "_" + suffix + ".pdf";
    c->SaveAs(pdfPath.c_str());
    std::string pngPath = outputDir + "/" + inputBasename + "_" + suffix + ".png";

    c->SaveAs(pngPath.c_str());

}


    // Extract directory and filename from input path
    // std::string inputPath = "/home/nick/Downloads/data1.csv";
    // std::string inputPath = file;
    size_t lastSlash = inputPath.rfind('/');
    // std::string inputDir      = inputPath.substr(0, lastSlash);
    // std::string inputFilename = inputPath.substr(lastSlash + 1);
    // Strip extension from filename
    size_t dotPos = inputFilename.rfind('.');
    // std::string inputBasename = (dotPos != std::string::npos) ? inputFilename.substr(0, dotPos) : inputFilename;

    // Create output folder path and make the directory


    // Open output file
    std::string outputFile = outputDir + "/" + inputBasename + "_fit_results.txt";
    std::ofstream results(outputFile);
    if (!results.is_open()) {
        std::cerr << "Error opening results file!" << std::endl;
        return;
    }

    // Write header
    results << "==================================================\n";
    results << "Fit Results for: " << inputFilename << "\n";
    results << "==================================================\n\n";

    // Re-run fits and write results for each pair
    for (int k = 0; k < (int)pairs.size(); k++) {
        int i = pairs[k].first;
        int j = pairs[k].second;
        std::string chI = std::to_string(channelIDs[i]);
        std::string chJ = std::to_string(channelIDs[j]);
        std::string suffix  = "ch" + chI + "_vs_ch" + chJ;
        std::string proName = "pro_" + suffix;
        std::string fitName = "fit_" + suffix;

        // Retrieve the profile and fit by name from ROOT's directory
        TProfile *pro = (TProfile*)gDirectory->FindObject(proName.c_str());
        TF1      *fit = (TF1*)gDirectory->FindObject(fitName.c_str());

        results << "--------------------------------------------------\n";
        results << "Channel Pair: " << chI << " vs " << chJ << "\n";
        results << "--------------------------------------------------\n";

        if (!pro || !fit) {
            results << "  WARNING: Could not retrieve fit objects for this pair.\n\n";
            continue;
        }

        // Fit quality
        double chi2   = fit->GetChisquare();
        int    ndf    = fit->GetNDF();
        double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
        double pValue = TMath::Prob(chi2, ndf);

        results << std::fixed << std::setprecision(6);
        results << "  Fit Function : pol1 (y = p0 + p1*x)\n";
        results << "  Fit Range    : " << fit->GetXmin() << " to " << fit->GetXmax() << "\n\n";

        results << "  Chi^2        : " << chi2     << "\n";
        results << "  NDF          : " << ndf      << "\n";
        results << "  Chi^2 / NDF  : " << chi2ndf  << "\n";
        results << "  P-Value      : " << pValue   << "\n\n";

        // Parameter results
        results << "  Parameters:\n";
        for (int p = 0; p < fit->GetNpar(); p++) {
            results << "    p" << p << " (" << fit->GetParName(p) << ")"
                    << "  = " << fit->GetParameter(p)
                    << " +/- " << fit->GetParError(p) << "\n";
        }

        // Profile histogram info
        results << "\n  Profile Histogram Info:\n";
        results << "    Entries : " << pro->GetEntries()  << "\n";
        results << "    Mean X  : " << pro->GetMean(1)    << "\n";
        results << "    Mean Y  : " << pro->GetMean(2)    << "\n";
        results << "    RMS X   : " << pro->GetRMS(1)     << "\n";
        results << "    RMS Y   : " << pro->GetRMS(2)     << "\n\n";
    }

    results << "==================================================\n";
    results << "End of Results\n";
    results << "==================================================\n";
    results.close();

    std::cout << "Fit results written to: " << outputFile << std::endl;



   

}
