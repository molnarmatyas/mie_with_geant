#include <fstream>
#include <sstream>
#include <string>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>

#define NDEG 1
double pixel = 0.02256; // in mm
pixel = 4*5.86 / 1000.0; // from notes, in mm


void SetGrayscalePalette(TCanvas* c, TH2D* hCCD) {
    const Int_t NCont = 255; // Number of colors in the palette
    gStyle->SetNumberContours(NCont);
    
    const Int_t nColors = 2; // Only black and white
    Double_t stops[nColors] = {0.0, 1.0}; // Position of colors in the range
    Double_t red[nColors]   = {0.0, 1.0}; // 1.0, 0.0: White to Black <--> 0.0, 1.0 black to white
    Double_t green[nColors] = {0.0, 1.0};
    Double_t blue[nColors]  = {0.0, 1.0};
    
    //Int_t palette = TColor::CreateGradientColorTable(nColors, stops, red, green, blue, NCont);
    //gStyle->SetPalette(palette);
    TColor::CreateGradientColorTable(nColors, stops, red, green, blue, NCont);

    c->SetFillColor(kBlack);  // Set canvas background color to black
    gPad->SetFrameFillColor(kBlack); // Set frame (plot background) to black
    //gStyle->SetMinimum(1e-6);  // Ensures empty bins are shown in the color scale
    
    hCCD->GetXaxis()->SetLabelColor(kWhite);
    hCCD->GetXaxis()->SetTitleColor(kWhite);
    hCCD->GetXaxis()->SetAxisColor(kWhite);
    hCCD->GetYaxis()->SetLabelColor(kWhite);
    hCCD->GetYaxis()->SetTitleColor(kWhite);
    hCCD->GetYaxis()->SetAxisColor(kWhite);
    hCCD->GetZaxis()->SetLabelColor(kWhite);
    hCCD->GetZaxis()->SetTitleColor(kWhite);
}
// --- MAIN ---
void txtToHist(std::string geantoutputname = "no_cell_measurement2_backgrond_extended_source.txt", std::string outputprefix = "cell_in_cell_beamprofiler_all_sensors_background_removed") {
  // Input text file path
  for(int ideg = 0; ideg < NDEG; ideg++)
  {
    std::string background_file =
        "../no_cell_measurement_backgrond_extended_source.txt";
    TH2D* dh2D_xz[2];
    //std::string filename =
    //    Form("../fs4_results/outputcrosssection_radians_poli_15_10000.txt7500_990_0_%i.txt", ideg);
    std::string filename = geantoutputname;
    //std::string outputprefix =
    //    Form("1M_3D_modell_homogeneous_cell_beamprofiler_source_y_shift_%i", ideg);

    TH1D* dhTheta = new TH1D("dhTheta", "dhTheta", 1000, 0, TMath::Pi());
    TH1D* dhGenTheta = new TH1D("dhGenTheta", "dhGenTheta", 1000, 0, TMath::Pi()/2);
    TH1D* dhXSect = new TH1D("dhXSect", "dhXSect", 1000, 0, TMath::Pi());
    TH1D* dhR = new TH1D("dhR", "dhR", 1000, 0, 10);
    TH1D* dhPosX = new TH1D("dhPosX", "dhPosX", 1000, 0, 10);
    TH1D* dhPosY = new TH1D("dhPosY", "dhPosY", 1000, 0, 10);
    TH1D* dhPosZ = new TH1D("dhPosZ", "dhPosZ", 1000, 0, 10);
    TH2D* dh2D_yz = new TH2D("dh2D_yz", "; Y [mm]; Z [mm]", 500, -10, 10, 500, -10, 10);
    dh2D_xz[0] = new TH2D("dh2D_xz_det_1", "; X [mm]; Z [mm]", 500, -2, 2, 500, -2, 2);
    dh2D_xz[1] = new TH2D("dh2D_xz_det_2", "; X [mm]; Z [mm]", 500, -2, 2, 500, -2, 2);
    TH2D* dh2D_xy = new TH2D("dh2D_xy", "; CCD X [pixel]; Y [pixel]", 480, -11.25/2.0/pixel, 11.25/2.0/pixel, 300, -7.03/2.0/pixel, 7.03/2.0/pixel);
    TH2D* dh2D_xy_ccd_background = new TH2D("dh2D_xy_ccd_background", "; CCD background X [pixel]; Y [pixel]", 480, -11.25/2.0/pixel, 11.25/2.0/pixel, 300, -7.03/2.0/pixel, 7.03/2.0/pixel);
    TH2D* dhR_alpha = new TH2D("dhR_alpha", "#alpha vs R; #alpha [deg] ; R [pixel]", 1000, -1, 180, 1000/pixel, 0, 7.5/pixel);
    TH2D* dhR_alpha_det_1 = new TH2D("dhR_alpha_det_1", "vbpw34s_1 #alpha vs R; #alpha [deg] ; R [pixel]", 100, -1, 30, 1000, 0, 3);
    TH2D* dhR_alpha_det_2 = new TH2D("dhR_alpha_det_2", "vbpw34s_2 #alpha vs R; #alpha [deg] ; R [pixel]", 100, -1, 30, 1000, 0, 3);
    TH2D* dhphi_alpha_det_1 = new TH2D("dhphi_alpha_det_1", "#phi vs #alpha; #phi; #alpha", 1000, -1.0*180, 180, 1000, 0.0, 180);
    TH2D* dhphi_alpha_det_2 = new TH2D("dhphi_alpha_det_2", "#phi vs #alpha; #phi; #alpha", 1000, -1.0*180, 180, 1000, 0.0, 180);
    TH2D* dh2D_rx_tantheta = new TH2D("dh2D_rx_tantheta", "dh2D_rx_tantheta", 1000, 0, 3, 1000, 0, 3);
    TH3D* dh3D_xyz = new TH3D("dh3D_xyz", "; CCD X [mm]; Y [mm]; Z [mm]", 100,5,20, 100,90,100, 100,-110,-100);

    // Open the input file
    std::ifstream infile(filename);
    if (!infile.is_open()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
    }
    std::ifstream bkg_file(background_file);
    if (!bkg_file.is_open()) {
      std::cerr << "Error opening file: " << background_file << std::endl;
      return;
    }
    /*
      Simple model sensor location
    double x_min = 8.478684;
    double y_min = 92.734001;
    double z_min = -107.531036;
    double x_max = 18.715816;
    double y_max = 99.765999;
    double z_max = -101.050804;
    */
    /*
      complete model sensor location
    */
    double x_min = 8.22955;
    double y_min = 92.734;
    double z_min = -106.676;
    double x_max = 18.2223;
    double y_max = 99.764;
    double z_max = -100.618;
    double x_center = (x_min + x_max ) / 2.0;
    double y_center = (y_min + y_max ) / 2.0;
    double z_center = (z_min + z_max ) / 2.0;

      
    double det_1_x_min = 48.3132;
    double det_1_y_min = 106.42;
    double det_1_z_min = -139.75;
    double det_1_x_max = 51.3132;
    double det_1_y_max = 106.56;
    double det_1_z_max = -136.75;
    double det_1_x_center = (det_1_x_min + det_1_x_max ) / 2.0;
    double det_1_y_center = (det_1_y_min + det_1_y_max ) / 2.0;
    double det_1_z_center = (det_1_z_min + det_1_z_max ) / 2.0;

    double det_2_x_min = 61.3132;
    double det_2_y_min = 106.42;
    double det_2_z_min = -139.75;
    double det_2_x_max = 64.3132;
    double det_2_y_max = 106.56;
    double det_2_z_max = -136.75;
    double det_2_x_center = (det_2_x_min + det_2_x_max ) / 2.0;
    double det_2_y_center = (det_2_y_min + det_2_y_max ) / 2.0;
    double det_2_z_center = (det_2_z_min + det_2_z_max ) / 2.0;

    // Define rotation angle (convert to radians)
    //double M_PI = TMath::Pi();
    double theta = 30.0 * M_PI / 180.0; 
    double cosTheta = TMath::Cos(theta);
    double sinTheta = TMath::Sin(theta);

    std::string bkg_line;
    while (std::getline(bkg_file, bkg_line)) {
      // Skip empty lines
      if (bkg_line.empty()) continue;

      std::replace(bkg_line.begin(), bkg_line.end(), ',', ' ');

      // Parse the line
      std::istringstream iss(bkg_line);
      double theta3, genTheta, R, postX, postY, postZ, alpha, alpha_man, phi3;
      double localR, localpostX, localpostY, localpostZ;
      double shiftedX, shiftedY, shiftedZ;
      int det_num;

      // Read all 6 columns
      if (!(iss >> theta3 >> genTheta >> R >> postX >> postY >> postZ >> alpha >> alpha_man >> phi3 >> det_num)) {
        std::cerr << "Error parsing line: " << bkg_line << std::endl;
        continue;
      }
      if(det_num == 2) //vbpw34s_2
      {
        // Translate to local origin
        double shiftedX = postX - det_2_x_center;
        double shiftedY = postY - det_2_y_center;
        double shiftedZ = postZ - det_2_z_center;

        continue;
      }
      else if(det_num == 1) //vbpw34s_1
      {
        // Translate to local origin
        double shiftedX = postX - det_1_x_center;
        double shiftedY = postY - det_1_y_center;
        double shiftedZ = postZ - det_1_z_center;

        continue;
      }
      else  // CCD
      {
        // Translate to local origin
        double shiftedX = postX - x_center;
        double shiftedY = postY - y_center;
        double shiftedZ = postZ - z_center;
        // Apply 2D rotation in the XZ plane
        localpostX = cosTheta * shiftedX + sinTheta * shiftedZ;
        localpostZ = -sinTheta * shiftedX + cosTheta * shiftedZ;
        // Y is not affected by rotation
        localpostY = postY - y_center;

        //localpostX = (x_center - postX);
        //localpostY = (y_center - postY);
        //localpostZ = (z_center - postZ);
        //std::cout << "x min: " << cosTheta * x_min + sinTheta * z_min << std::endl;
        //std::cout << "x max: " << cosTheta * x_max + sinTheta * z_max << std::endl;

        //std::cout << "z min: " << cosTheta * z_min - sinTheta * x_min << std::endl;
        //std::cout << "z max: " << cosTheta * z_max - sinTheta * x_max << std::endl;

        //std::cout << "y min: " << y_min << std::endl;
        //std::cout << "y max: " << y_max << std::endl;
        //sleep(5);

        dh2D_xy_ccd_background->Fill(localpostX/pixel, localpostY/pixel);

      }
    }

    // Read file line by line
    std::string line;
    while (std::getline(infile, line)) {
      // Skip empty lines
      if (line.empty()) continue;

      std::replace(line.begin(), line.end(), ',', ' ');

      // Parse the line
      std::istringstream iss(line);
      double theta3, genTheta, R, postX, postY, postZ, alpha, alpha_man, phi3;
      double localR, localpostX, localpostY, localpostZ;
      double shiftedX, shiftedY, shiftedZ;
      int det_num;

      // Read all 6 columns
      if (!(iss >> theta3 >> genTheta >> R >> postX >> postY >> postZ >> alpha >> alpha_man >> phi3 >> det_num)) {
        std::cerr << "Error parsing line: " << line << std::endl;
        continue;
      }
      if(det_num == 2) //vbpw34s_2
      {
        // Translate to local origin
        double shiftedX = postX - det_2_x_center;
        double shiftedY = postY - det_2_y_center;
        double shiftedZ = postZ - det_2_z_center;
        localR = sqrt(shiftedZ * shiftedZ + shiftedX * shiftedX);
        dhR_alpha_det_2->Fill(genTheta*180.0 / TMath::Pi(), localR);
        dh2D_xz[1]->Fill(shiftedZ, shiftedX);
        dhphi_alpha_det_2->Fill(phi3, genTheta*180.0 / TMath::Pi());

        continue;
      }
      else if(det_num == 1) //vbpw34s_1
      {
        // Translate to local origin
        double shiftedX = postX - det_1_x_center;
        double shiftedY = postY - det_1_y_center;
        double shiftedZ = postZ - det_1_z_center;
        localR = sqrt(shiftedZ * shiftedZ + shiftedX * shiftedX);
        dhR_alpha_det_1->Fill(genTheta*180.0 / TMath::Pi(), localR);
        dh2D_xz[0]->Fill(shiftedZ, shiftedX);
        dhphi_alpha_det_1->Fill(phi3, genTheta*180.0 / TMath::Pi());

        continue;
      }
      else  // CCD
      {
        // Translate to local origin
        double shiftedX = postX - x_center;
        double shiftedY = postY - y_center;
        double shiftedZ = postZ - z_center;
        // Apply 2D rotation in the XZ plane
        localpostX = cosTheta * shiftedX + sinTheta * shiftedZ;
        localpostZ = -sinTheta * shiftedX + cosTheta * shiftedZ;
        // Y is not affected by rotation
        localpostY = postY - y_center;

        //localpostX = (x_center - postX);
        //localpostY = (y_center - postY);
        //localpostZ = (z_center - postZ);

        localR = sqrt(localpostY * localpostY + localpostX * localpostX) / pixel;
        dhR_alpha->Fill(genTheta*180.0 / TMath::Pi(), localR);
        dhR->Fill(localR);
        dhPosX->Fill(localpostX);
        dhPosY->Fill(localpostY);
        dhPosZ->Fill(localpostZ);
        dh2D_yz->Fill(localpostZ, localpostY);
        dh2D_xy->Fill(localpostX/pixel, localpostY/pixel);
        dh3D_xyz->Fill(postX,postY,postZ);
        dh2D_rx_tantheta->Fill(localR / localpostX, std::tan(theta3));

      }



      // Fill the histogram with the last two columns
      dhTheta->Fill(theta3);
      dhGenTheta->Fill(genTheta);
      dhXSect->Fill(theta3);
    }
    //  dh2D_yz->GetZaxis()->SetRangeUser(0,10);
    TF1* oneoversin = new TF1("oneoversin","1/sin(x)");
    dhXSect->Multiply(oneoversin);

    TF1* Ltanalpha = new TF1("Ltanalpha","9.0*tan(x)");
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
    TFile* outfile = new TFile(Form("%s_%i_alpha_output.root",
                                    outputprefix.c_str(), ideg), "RECREATE");

    // Optional: draw and save as image
    TCanvas* c1 = new TCanvas("c1", "", 800, 800);
    gStyle->SetOptStat(0);
    //dhR_alpha->Draw("COLZ");
    //Ltanalpha->Draw("same");
    //c1->SaveAs(Form("figs/%s_mirror_dhR_alpha_output3500_um_pointsource_1M_pld_1592_ver_thetafix.pdf", outputprefix.c_str()));

    //gStyle->SetCanvasDefH(550);
    //gStyle->SetCanvasDefW(650);
    c1->SetLogz(1);


    dh2D_xy->SetTitle(Form("2D scattering, deg=%i, n = 1.592, 3D model, d = 7 #mu m", ideg));
    dh2D_xy->GetXaxis()->SetTitle("X [mm]");
    dh2D_xy->GetYaxis()->SetTitle("Y [mm]");
    //dh2D_xy->Draw("COLZ");
    //c1->SaveAs(Form("figs/%s_dh2D_xy_colorful.png", outputprefix.c_str()));
    dh2D_xy->Write();

    dh2D_xz[0]->SetTitle(Form("2D scattering, deg=%i, vbpw34s_1", ideg));
    dh2D_xz[0]->GetXaxis()->SetTitle("X [mm]");
    dh2D_xz[0]->GetYaxis()->SetTitle("Z [mm]");
    //dh2D_xz[0]->Draw("COLZ");

    dh2D_xz[1]->SetTitle(Form("2D scattering, deg=%i, vbpw34s_2", ideg));
    dh2D_xz[1]->GetXaxis()->SetTitle("X [mm]");
    dh2D_xz[1]->GetYaxis()->SetTitle("Z [mm]");
    //dh2D_xz[1]->Draw("COLZ");
    //c1->SaveAs(Form("../figs/%s_dh2D_xy_discrete.pdf", outputprefix.c_str()));

    // Black & white image of CCD screen (black=min, white=max)
    TCanvas *c2 = new TCanvas("c2", "CCD Image", 800, 600);
    SetGrayscalePalette(c2, dh2D_xy); // Apply grayscale palette
    //c2->SetLogz(1);

    //dh2D_xy_ccd_background->Scale(1.0/dh2D_xy_ccd_background->Integral());
    //dh2D_xy->Scale(1.0/dh2D_xy->Integral());


    dh2D_xy->GetXaxis()->SetTitle("X [pixel]");
    dh2D_xy->GetYaxis()->SetTitle("Y [pixel]");
    dh2D_xy->Draw("COLZ");
    //c2->SaveAs(Form("figs/%s_dh2D_xy_bw.png", outputprefix.c_str()));
    // Subtract background
    dh2D_xy->Add(dh2D_xy_ccd_background, -1.0);


    // NORMALISATION
    // Option 1: Set negative bins to zero
    for (int ix = 1; ix <= dh2D_xy->GetNbinsX(); ++ix) {
      for (int iy = 1; iy <= dh2D_xy->GetNbinsY(); ++iy) {
        double val = dh2D_xy->GetBinContent(ix, iy);
        if (val < 0.0) dh2D_xy->SetBinContent(ix, iy, 0.0);
      }
    }

    // Write non-negative intensity map to CSV (pixels coordinates)
    {
      std::string csvname = Form("%s_%i_intensity_map_non_normalized.csv", outputprefix.c_str(), ideg);
      std::ofstream csv(csvname);
      if (!csv.is_open()) {
        std::cerr << "Error opening CSV file: " << csvname << std::endl;
      } else {
        csv << "x_pixel,y_pixel,intensity\n";
        for (int ix = 1; ix <= dh2D_xy->GetNbinsX(); ++ix) {
          double x = dh2D_xy->GetXaxis()->GetBinCenter(ix);
          for (int iy = 1; iy <= dh2D_xy->GetNbinsY(); ++iy) {
            double y = dh2D_xy->GetYaxis()->GetBinCenter(iy);
            double val = dh2D_xy->GetBinContent(ix, iy);
            csv << std::fixed << std::setprecision(6) << x << "," << y << "," << val << "\n";
          }
        }
        csv.close();
        std::cout << "Wrote non-negative intensity map of "<<dh2D_xy->GetNbinsX()<<"x"<<dh2D_xy->GetNbinsY()<<" CCD to " << csvname << std::endl;
      }
    }

    // Write non-negative pixel map to CSV (pixel indices corresp. to table indices)
    {
      std::string csvname = Form("%s_%i_pixel_map_non_normalized.csv", outputprefix.c_str(), ideg);
      std::ofstream csv(csvname);
      if (!csv.is_open()) {
        std::cerr << "Error opening CSV file: " << csvname << std::endl;
      } else {
        for (int iy = dh2D_xy->GetNbinsY(); iy>=1; iy--) { // Reverse Y for image-like orientation
          for (int ix = 1; ix <= dh2D_xy->GetNbinsX(); ++ix) {
            double val = dh2D_xy->GetBinContent(ix, iy);
            csv << std::fixed << std::setprecision(6) << val;
            if (ix < dh2D_xy->GetNbinsX()) csv << ",";
          }
          csv << "\n";
        }
        csv.close();
        std::cout << "Wrote non-negative pixel map of "<<dh2D_xy->GetNbinsX()<<"x"<<dh2D_xy->GetNbinsY()<<" CCD to " << csvname << std::endl;
      }
    }

    // Option 2: Normalise so that 0.0 is the smallest value
    // (Uncomment if you want this option)
    /*
    double minVal = dh2D_xy->GetMinimum();
    if (minVal != 0.0) {
      for (int ix = 1; ix <= dh2D_xy->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= dh2D_xy->GetNbinsY(); ++iy) {
          double val = dh2D_xy->GetBinContent(ix, iy);
          dh2D_xy->SetBinContent(ix, iy, val - minVal);
        }
      }
    }
    */
    c2->SaveAs(Form("figs/%s_dh2D_xy_bw_bg_removed.png", outputprefix.c_str()));


    // R vs alpha
    TProfile* prof = dhR_alpha->ProfileX("_prof_max", 1, -1, "");
    prof->GetYaxis()->SetTitle("R [pixel]");
    prof->GetXaxis()->SetTitle("#alpha [deg]");



    prof->Write();
    dhTheta->Write();
    dhGenTheta->Write();
    dhXSect->Write();
    dhR->Write();
    dhPosX->Write();
    dhPosY->Write();
    dhPosZ->Write();
    dh2D_yz->Write();
    dhR_alpha->Write();
    dhR_alpha_det_1->Write();
    dh2D_xz[0]->Write();
    dh2D_xz[1]->Write();
    dhR_alpha_det_2->Write();
    dhphi_alpha_det_1->Write();
    dhphi_alpha_det_2->Write();
    dh2D_rx_tantheta->Write();
    dh3D_xyz->Write();

    outfile->Close();
  }
}
