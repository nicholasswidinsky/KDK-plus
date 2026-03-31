#include <fstream>
#include <sstream>
#include <iostream>
#include "TCanvas.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCutG.h"

void fit2D_cut() {

    TCanvas *c = new TCanvas("c","2D Histogram with Banana Cut",1000,800);
    c->Divide(1,1); 

    TH2D *h2 = new TH2D("h2"," ;NaI scintillator;Liquid scintillator", 100,0,4000,100,0,4000);

    TProfile *pro = new TProfile("pro","Profile Histogram;NaI scintillator;Liquid scintillator", 100,0,4000,0,4000,"");

    std::ifstream file("/home/nick/Downloads/data1.csv");
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        double x, y;
        char comma;
        if (ss >> x >> comma >> y) {
            h2->Fill(x,y);
        }
    }
    file.close();

    c->cd(1);
    h2->Draw("COLZ");

    TCutG *cut = new TCutG("banana",3);  
    cut->SetPoint(0, 1000, 200);
    cut->SetPoint(1, 000, 2500);
    cut->SetPoint(2, 000, 200);
    cut->SetPoint(3, 1000, 200);
    //cut->SetPoint(4, 3000, 1600);
    //cut->SetPoint(5, 3500, 1800);
    //cut->SetPoint(6, 3500, 2000);
    //cut->SetPoint(7, 1000, 2000);
    //cut->SetPoint(8, 1000, 500); 
    cut->SetLineColor(kRed+2);
    cut->SetLineWidth(3);
    cut->Draw("same");


    for (int i = 1; i <= h2->GetNbinsX(); i++) {
        for (int j = 1; j <= h2->GetNbinsY(); j++) {
            double x = h2->GetXaxis()->GetBinCenter(i);
            double y = h2->GetYaxis()->GetBinCenter(j);
            double z = h2->GetBinContent(i,j);
            if (!cut->IsInside(x,y) && z > 0) {
                pro->Fill(x,y,z);
            }
            //pro->Fill(x,y,z);
        }
    }

    pro->Draw("same");


    TF1 *fit = new TF1("fit","pol1",500,1100);  
    pro->Fit(fit,"R");
    fit->SetLineColor(kRed);
    
    fit->Draw("same");


}
