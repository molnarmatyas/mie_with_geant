#include <fstream>
#include <sstream>
#include <string>
#include <TH2D.h>
#include <TFile.h>

void txtToHist() {
    // Input text file path
    std::string filename = "../build/output3500_990_0_6.txt";
    
    TH1D* dhTheta = new TH1D("dhTheta", "dhTheta", 1000, 0, TMath::Pi());
    TH1D* dhGenTheta = new TH1D("dhGenTheta", "dhGenTheta", 1000, 0, TMath::Pi());
    TH1D* dhXSect = new TH1D("dhXSect", "dhXSect", 1000, 0, TMath::Pi());
    TH1D* dhR = new TH1D("dhR", "dhR", 1000, 0, 9 * 200);
    TH1D* dhPosX = new TH1D("dhPosX", "dhPosX", 1000, 0, 9);
    TH1D* dhPosY = new TH1D("dhPosY", "dhPosY", 1000, 0, 9);
    TH1D* dhPosZ = new TH1D("dhPosZ", "dhPosZ", 1000, 0, 9);
    TH2D* dh2D_yz = new TH2D("dh2D_yz", "dh2D_yz", 1000, -5, 5, 1000, -5, 5);
    TH2D* dhR_alpha = new TH2D("dhR_alpha", "dhR_alpha", 1000, 0, TMath::Pi(), 1000, 0, 9 * 200);
    
    // Open the input file
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    // Read file line by line
    std::string line;
    while (std::getline(infile, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        std::replace(line.begin(), line.end(), ',', ' ');

        // Parse the line
        std::istringstream iss(line);
        double theta3, genTheta, R, postX, postY, postZ;
        
        // Read all 6 columns
        if (!(iss >> theta3 >> genTheta >> R >> postX >> postY >> postZ)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue;
        }
        
        // Fill the histogram with the last two columns
        dhTheta->Fill(theta3);
        dhGenTheta->Fill(genTheta);
        dhR->Fill(R * 200.0);
        dhPosX->Fill(postX);
        dhPosY->Fill(postY);
        dhPosZ->Fill(postZ);
        dhXSect->Fill(theta3 / std::sin(theta3));
        dhR_alpha->Fill(theta3, R * 200.0);
        dh2D_yz->Fill(postY, postZ);
    }
//    dh2D_yz->GetZaxis()->SetRangeUser(0,10);
    
    // Create a ROOT output file to save the histogram
    TFile* outfile = new TFile("output.root", "RECREATE");
    
    // Optional: draw and save as image
    TCanvas* c1 = new TCanvas("c1", "", 800, 800);
    gStyle->SetCanvasDefH(550);
    gStyle->SetCanvasDefW(650);
    c1->SetLogz(1);
    dh2D_yz->Draw("COLZ");
    c1->SaveAs("dh2D_yz_output3500_um_pointsource_lowstat_pld_fixed_10_1592.png");





    dhTheta->Write();
    dhGenTheta->Write();
    dhXSect->Write();
    dhR->Write();
    dhPosX->Write();
    dhPosY->Write();
    dhPosZ->Write();
    dh2D_yz->Write();
    dhR_alpha->Write();
    outfile->Close();
}
