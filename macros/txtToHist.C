#include <fstream>
#include <sstream>
#include <string>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>

void txtToHist() {
    // Input text file path
    std::string filename = "../build/output3500_990_0.txt";
    std::string outputprefix = "3D_modell_camera_test_1";
    
    TH1D* dhTheta = new TH1D("dhTheta", "dhTheta", 1000, 0, TMath::Pi());
    TH1D* dhGenTheta = new TH1D("dhGenTheta", "dhGenTheta", 1000, 0, TMath::Pi());
    TH1D* dhXSect = new TH1D("dhXSect", "dhXSect", 1000, 0, TMath::Pi());
    TH1D* dhR = new TH1D("dhR", "dhR", 1000, 0, 10 * 200);
    TH1D* dhPosX = new TH1D("dhPosX", "dhPosX", 1000, 0, 10);
    TH1D* dhPosY = new TH1D("dhPosY", "dhPosY", 1000, 0, 10);
    TH1D* dhPosZ = new TH1D("dhPosZ", "dhPosZ", 1000, 0, 10);
    TH2D* dh2D_yz = new TH2D("dh2D_yz", "; Y [mm]; Z [mm]", 100, -10, 10, 100, -10, 10);
    TH2D* dh2D_xy = new TH2D("dh2D_xy", "; X [mm]; Y [mm]", 150, -6, 6, 120, -4, 4);
    TH2D* dhR_alpha = new TH2D("dhR_alpha", "#alpha vs R; #alpha; R [pixel]", 1000, 0, TMath::Pi() / 2.0, 1000, 0, 7.5 * 200);
    TH2D* dhtheta_alpha = new TH2D("dhtheta_alpha", "#theta vs #alpha; #theta; #alpha", 1000, 0, TMath::Pi(), 1000, 0, TMath::Pi());
    TH2D* dh2D_rx_tantheta = new TH2D("dh2D_rx_tantheta", "dh2D_rx_tantheta", 1000, 0, 3, 1000, 0, 3);
    TH3D* dh3D_xyz = new TH3D("dh3D_xyz", "; X [mm]; Y [mm]; Z [mm]", 100,5,20, 100,90,100, 100,-110,-100);
    
    // Open the input file
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    double x_min = 8.478684;
    double y_min = 92.734001;
    double z_min = -107.531036;
    double x_max = 18.715816;
    double y_max = 99.765999;
    double z_max = -101.050804;
    double x_center = (x_min + x_max ) / 2.0;
    double y_center = (y_min + y_max ) / 2.0;
    double z_center = (z_min + z_max ) / 2.0;

    // Define rotation angle (convert to radians)
    //double M_PI = TMath::Pi();
    double theta = 25.0 * M_PI / 180.0; // FIXME this is not the real angle; how on earth can I measure angle in FreeCAD??? 
    double cosTheta = TMath::Cos(theta);
    double sinTheta = TMath::Sin(theta);
    
    // Read file line by line
    std::string line;
    while (std::getline(infile, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        std::replace(line.begin(), line.end(), ',', ' ');

        // Parse the line
        std::istringstream iss(line);
        double theta3, genTheta, R, postX, postY, postZ, alpha, alpha_man;
        double localR, localpostX, localpostY, localpostZ;
        
        // Read all 6 columns
        if (!(iss >> theta3 >> genTheta >> R >> postX >> postY >> postZ >> alpha >> alpha_man)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue;
        }

        // Translate to local origin
        double shiftedX = postX - x_center;
        double shiftedZ = postZ - z_center;

        // Apply 2D rotation in the XZ plane
        localpostX = cosTheta * shiftedX + sinTheta * shiftedZ;
        localpostZ = -sinTheta * shiftedX + cosTheta * shiftedZ;
        // Y is not affected by rotation
        localpostY = postY - y_center;

        //localpostX = (x_center - postX);
        //localpostY = (y_center - postY);
        //localpostZ = (z_center - postZ);
        
        localR = sqrt(localpostY * localpostY + localpostZ * localpostZ);
        
        // Fill the histogram with the last two columns
        dhTheta->Fill(theta3);
        dhGenTheta->Fill(genTheta);
        dhR->Fill(localR * 200.0);
        dhPosX->Fill(localpostX);
        dhPosY->Fill(localpostY);
        dhPosZ->Fill(localpostZ);
        dhXSect->Fill(theta3);
        dhR_alpha->Fill(genTheta, localR * 200.0);
        dhtheta_alpha->Fill(theta3, alpha);
        dh2D_yz->Fill(localpostZ, localpostY);
        dh2D_xy->Fill(localpostX, localpostY);
        dh3D_xyz->Fill(postX,postY,postZ);
        dh2D_rx_tantheta->Fill(localR / localpostX, std::tan(theta3));
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
    //Ltanalpha->Draw("same");
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
    dh2D_xy->Write();
    dhR_alpha->Write();
    dhtheta_alpha->Write();
    dh2D_rx_tantheta->Write();
    dh3D_xyz->Write();

    outfile->Close();
}
