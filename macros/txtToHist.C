#include <fstream>
#include <sstream>
#include <string>
#include <TH2D.h>
#include <TFile.h>

void txtToHist() {
    // Input text file path
    std::string filename = "../build/output3500_990_0_149.txt";
    std::string outputprefix = "shield_02mm";
    
    TH1D* dhTheta = new TH1D("dhTheta", "dhTheta", 1000, 0, TMath::Pi());
    TH1D* dhGenTheta = new TH1D("dhGenTheta", "dhGenTheta", 1000, 0, TMath::Pi());
    TH1D* dhXSect = new TH1D("dhXSect", "dhXSect", 1000, 0, TMath::Pi());
    TH1D* dhR = new TH1D("dhR", "dhR", 1000, 0, 10 * 200);
    TH1D* dhPosX = new TH1D("dhPosX", "dhPosX", 1000, 0, 10);
    TH1D* dhPosY = new TH1D("dhPosY", "dhPosY", 1000, 0, 10);
    TH1D* dhPosZ = new TH1D("dhPosZ", "dhPosZ", 1000, 0, 10);
    TH2D* dh2D_yz = new TH2D("dh2D_yz", "; Y [mm]; Z [mm]", 1000, -5, 5, 1000, -5, 5);
    TH2D* dhR_alpha = new TH2D("dhR_alpha", "#alpha vs R; #alpha; R [pixel]", 1000, 0, TMath::Pi()/2.0, 1000, 0, 15 * 200);
    TH2D* dhtheta_alpha = new TH2D("dhtheta_alpha", "#theta vs #alpha; #theta; #alpha", 1000, 0, TMath::Pi(), 1000, 0, TMath::Pi());
    TH2D* dh2D_rx_tantheta = new TH2D("dh2D_rx_tantheta", "dh2D_rx_tantheta", 1000, 0, 3, 1000, 0, 3);
    
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
        double theta3, genTheta, R, postX, postY, postZ, alpha, alpha_man;
        
        // Read all 6 columns
        if (!(iss >> theta3 >> genTheta >> R >> postX >> postY >> postZ >> alpha >> alpha_man)) {
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
        dhXSect->Fill(theta3);
        dhR_alpha->Fill(alpha, R * 200.0);
        dhtheta_alpha->Fill(theta3, alpha);
        dh2D_yz->Fill(postZ, postY);
        dh2D_rx_tantheta->Fill(R / postX, std::tan(theta3));
    }
//  dh2D_yz->GetZaxis()->SetRangeUser(0,10);
    TF1* oneoversin = new TF1("oneoversin","1/sin(x)");
    dhXSect->Multiply(oneoversin);

    TF1* Ltanalpha = new TF1("Ltanalpha","9.0*tan(x) * 200");
    Ltanalpha->SetLineWidth(1);
    Ltanalpha->SetLineColorAlpha(kRed,0.75);
    //setup for cross section
    /*
    int inpNevent = 1000000;
    int outNevent = dhTheta->GetEntries();
    double cellAreamm = 1.539e-4; 
    double Lmm = 10.0;
    double dSigmadOmega;
    for(int ievent = 0; ievent < dhTheta->GetEntries(); ievent++) {

       dSigmadOmega = outNevent * cellAreamm / inpNevent / (2*TMath::Pi() / Lmm * TMath::Cos(
    }
    */
    
    // Create a ROOT output file to save the histogram
    TFile* outfile = new TFile(Form("%s_alpha_output.root", outputprefix.c_str()), "RECREATE");
    
    // Optional: draw and save as image
    TCanvas* c1 = new TCanvas("c1", "", 800, 800);
    gStyle->SetOptStat(0);
    dhR_alpha->Draw("COLZ");
    Ltanalpha->Draw("same");
    c1->SaveAs(Form("%s_mirror_dhR_alpha_output3500_um_pointsource_1M_pld_1592_ver_thetafix.pdf", outputprefix.c_str()));

    //gStyle->SetCanvasDefH(550);
    //gStyle->SetCanvasDefW(650);
    c1->SetLogz(1);


    dh2D_yz->SetTitle("2D scattering, n = 1.592, distance = 9 mm, d = 7 #mu m");
    dh2D_yz->GetXaxis()->SetTitle("Z [mm]");
    dh2D_yz->GetYaxis()->SetTitle("Y [mm]");
    dh2D_yz->Draw("COLZ");
    c1->SaveAs(Form("%s_dh2D_yz_output3500_um_pointsource_1M_pld_1592_ver_thetafix.pdf", outputprefix.c_str()));




    dhTheta->Write();
    dhGenTheta->Write();
    dhXSect->Write();
    dhR->Write();
    dhPosX->Write();
    dhPosY->Write();
    dhPosZ->Write();
    dh2D_yz->Write();
    dhR_alpha->Write();
    dhtheta_alpha->Write();
    dh2D_rx_tantheta->Write();

    outfile->Close();
}
