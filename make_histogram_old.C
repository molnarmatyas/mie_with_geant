void make_histogram_old(const char* filename, int column_number, const char* xtitle, const char* pngname, double xmin=-999, double xmax=-999) {
    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    //c1->SetFillColor(42);
    c1->SetGrid();
    //c1->GetFrame()->SetFillColor(21);
    //c1->GetFrame()->SetBorderSize(12);

    TTree *tree = new TTree("tree","Tree");
    double value;
    tree->Branch("value",&value);

    ifstream file(filename);
    string line;
    while (getline(file,line)) {
        stringstream ss(line);
        for (int i=0; i<column_number; i++) {
            ss >> value;
        }
        tree->Fill();
    }

    double thisxmin = tree->GetMinimum("value");
    double thisxmax = tree->GetMaximum("value");
    if(xmin!=-999) thisxmin=xmin;
    if(xmax!=-999) thisxmax=xmax;
    TH1F *histo = new TH1F("statistics",Form("%s distribution",xtitle),100,xmin,xmax);
    histo->GetXaxis()->SetTitle(xtitle);
    histo->GetYaxis()->SetTitle("N [a.u.]");

    tree->Draw("value>>statistics");
    histo->SetMarkerStyle(20);
    histo->Draw("EP");
    histo->Print("all");
    c1->SetLogy(1);
    c1->Print(Form("figs/%s",pngname));

    if(strcmp(xtitle,"#theta")==0)
    {
      TF1* oneoversin = new TF1("oneoversin","1/sin(x)");
      double dtheta = (xmax-xmin)/100;
      double jbe = 10000/1.;
      histo->Multiply(oneoversin,1/((2*M_PI)*dtheta*jbe));
      histo->SetTitle("d#sigma/d#Omega");
      histo->GetYaxis()->SetTitle("d#sigma/d#Omega [{#mu}m^2]");
      gStyle->SetOptStat(0);
      histo->Draw("EP");
      TGraph* dsigmadOmega_refraction = new TGraph("dsdOtrue.txt"); 
      dsigmadOmega_refraction->SetLineColor(kRed);
      dsigmadOmega_refraction->SetLineWidth(2);
      dsigmadOmega_refraction->Draw("L");
      TLegend *leg = new TLegend(0.52,0.65,0.88,0.88);
      leg->AddEntry(histo,"GEANT4 simulation","EP");
      leg->AddEntry(dsigmadOmega_refraction,"Refraction calculation","L");
      leg->Draw();
      c1->Print("figs/crosssection.png");
    }
}
