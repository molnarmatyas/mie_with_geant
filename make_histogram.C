// root.exe -b -q 'make_histogram.C("build/output100M.root", 100, "#theta", "thetadist.png")'
void make_histogram(const char* filename, int column_number, const char* xtitle, const char* pngname, double xmin=-999, double xmax=-999) {
    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
    //c1->SetFillColor(42);
    c1->SetGrid();
    //c1->GetFrame()->SetFillColor(21);
    //c1->GetFrame()->SetBorderSize(12);

    TFile* file = new TFile(filename);
    TTree *tree = (TTree*)file->Get("Hits");
    double value;
    tree->SetBranchAddress("fTheta",&value);

    int Nbins = 200;
    TH1F *histo = new TH1F("statistics",Form("%s distribution",xtitle),Nbins,xmin,xmax);
    histo->GetXaxis()->SetTitle(xtitle);
    histo->GetYaxis()->SetTitle("N [a.u.]");

    for (int i=0; i<tree->GetEntries(); i++) 
    {
        tree->GetEntry(i);
        histo->Fill(value);
        cout<<i<<endl;
    }

    c1->SetLogy(1);
    histo->GetYaxis()->SetTitle("d#sigma/d#Omega [cm^{2}]");
    histo->SetMarkerStyle(6);
    gStyle->SetOptStat(0);
    histo->Draw("HIST P");
    c1->Print("figs/original.png");


    if(strcmp(xtitle,"#theta")==0)
    {
      TF1* oneoversin = new TF1("oneoversin","1/sin(x)");
      double dtheta = (xmax-xmin)/100;
      double jbe = 100000000/0.001;
      histo->Multiply(oneoversin);
      histo->Scale(1/((2*M_PI)*dtheta*jbe));
      histo->SetTitle("d#sigma/d#Omega");
      histo->GetYaxis()->SetTitle("d#sigma/d#Omega [cm^{2}]");
      histo->SetMarkerStyle(6);
      gStyle->SetOptStat(0);
      histo->Draw("HIST P");
      c1->Print("figs/crosssection.png");
    }
    
}
