#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>

void make_histogram(char* filename, int column_number, char* xtitle, char* pngname, double xmin=-999, double xmax=-999) {
    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    //c1->SetFillColor(42);
    c1->SetGrid();
    //c1->GetFrame()->SetFillColor(21);
    //c1->GetFrame()->SetBorderSize(12);
    std::ifstream file(filename);

    std::string line;
    std::vector< float > firsts;
    std::vector< float > seconds;
    while (std::getline(file, line)) {
        // Splitting each line by tab
        std::istringstream iss(line);
        float firstValue, secondValue;
        iss >> firstValue >> secondValue;
        firsts.push_back(firstValue);
        seconds.push_back(secondValue);
    }

    file.close();
    
    double thisxmin = *std::min_element(std::begin(seconds), std::end(seconds));
    double thisxmax = *std::max_element(std::begin(seconds), std::end(seconds));
    int N = seconds.size();
    if(xmin!=-999) thisxmin=xmin;
    if(xmax!=-999) thisxmax=xmax;
    TH1F *histo = new TH1F("statistics",Form("%s distribution",xtitle),500,xmin,xmax);
    for(int i = 0; i < N; i++) histo->Fill(seconds[i]);
    histo->SetBinContent(1,0.0);
    histo->GetXaxis()->SetTitle(xtitle);
    histo->GetYaxis()->SetTitle("N [a.u.]");

    histo->SetMarkerStyle(24);
    histo->SetMarkerSize(0.5);
//    histo->SetMarkerColorAlpha(kBlue, 0.35);
    histo->Sumw2();
    histo->Rebin(2);
    histo->Draw("EP");
    histo->Print("all");
    c1->SetLogy(1);
    c1->Print(Form("figs/%s",pngname));
}



int main(int argc, char* argv[]) {

    make_histogram(argv[1], atoi(argv[2]), argv[3], argv[4], atof(argv[5]), atof(argv[6]));
    return 0;
}

