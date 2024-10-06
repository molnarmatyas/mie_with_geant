void plotCrossSection(const char* filename, const char* xtitle, const char* pngname, double xmin=-999, double xmax=-999) {
    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    //c1->SetFillColor(42);
    c1->SetGrid();
    //c1->GetFrame()->SetFillColor(21);
    //c1->GetFrame()->SetBorderSize(12);

    TFile* file = new TFile(filename);
    TTree *tree = (TTree*)file->Get("Hits");

    double value;
    tree->SetBranchAddress("fTheta",&value);

    int nevents = tree->GetEntries();
    TH1D *histo = new TH1D("theta",Form("%s distribution",xtitle),200,xmin,xmax);
    for (int i=0; i< nevents; i++) {
        tree->GetEntry(i);
        histo->Fill(value);
    }
    histo->GetXaxis()->SetTitle(xtitle);
    histo->GetYaxis()->SetTitle("N [a.u.]");
    double smallAngleMin = 0.02635;
    double smallAngleMax = 0.1068;
    double smallAngleiPhotonCount = histo->Integral(histo->FindBin(smallAngleMin), histo->FindBin(smallAngleMax));

    double bigAngleMin = 0.2016;
    double bigAngleMax = 0.2815;
    double bigAnglePhotonCount = histo->Integral(histo->FindBin(bigAngleMin), histo->FindBin(bigAngleMax));

    std::cout << "Small (1.12°,6.51°) angle photon count: " << smallAngleiPhotonCount << "\t big (11.55°,16.13°) angle photon count: " << bigAnglePhotonCount << std::endl; 

    histo->SetMarkerStyle(24);
    histo->SetMarkerSize(0.5);

//    histo->Draw("EP");
//    histo->Print("all");
    c1->SetLogy(1);

    TF1* oneoversin = new TF1("oneoversin","1/sin(x)");
    double dtheta = (xmax-xmin)/200;
    double jbe = nevents/0.001;
    histo->SetBinContent(1, 0.0);
    histo->Multiply(oneoversin);
    histo->Scale(1/((2*M_PI)*dtheta*jbe));
    histo->SetTitle("d#sigma/d#Omega");
    histo->GetYaxis()->SetTitle("d#sigma/d#Omega [cm^{2}]");
    gStyle->SetOptStat(0);
    histo->Draw("HIST P");

    TLegend *leg = new TLegend(0.52,0.65,0.88,0.88);
    leg->AddEntry(histo,"GEANT4 simulation","EP");
    leg->Draw();
    c1->Print("figs/crosssection.png");
}
