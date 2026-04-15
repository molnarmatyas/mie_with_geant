void plot_det3_beamprofiler()
{
    // ROOT macro to plot BeamProfiler images from output files
    // Plots dh2D_xz_det_3 histograms with COLZ option
    // Saves all plots to a single PDF file
    
    gStyle->SetOptStat(0);
    
    TString outDir = gSystem->WorkingDirectory();
    TSystemDirectory dir(".", outDir);
    TList *files = dir.GetListOfFiles();
    
    // Collect output files
    vector<TString> smallAngleFiles;  // 0-4
    vector<TString> largeAngleFiles;  // 6-15
    
    TIter next(files);
    TSystemFile *file;
    while ((file = (TSystemFile*)next())) {
        TString fname = file->GetName();
        if (fname.Contains("outputtestdhkm") && fname.EndsWith(".root")) {
            if (fname.Contains("0-4")) {
                smallAngleFiles.push_back(outDir + "/" + fname);
            } else if (fname.Contains("6-15")) {
                largeAngleFiles.push_back(outDir + "/" + fname);
            }
        }
    }
    
    // Sort files
    sort(smallAngleFiles.begin(), smallAngleFiles.end());
    sort(largeAngleFiles.begin(), largeAngleFiles.end());
    
    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "BeamProfiler Det3 Plots", 800, 600);
    
    // Open PDF
    TString pdfName = outDir + "/BeamProfiler_Det3_Plots.pdf";
    c1->Print(pdfName + "[", "pdf");
    
    // Plot small angle files
    for (const auto &fname : smallAngleFiles) {
        TFile *f = TFile::Open(fname);
        if (!f || f->IsZombie()) continue;
        
        TH2D *hist = (TH2D*)f->Get("dh2D_xz_det_3");
        if (!hist) {
            f->Close();
            continue;
        }
        
        c1->Clear();
        hist->Draw("COLZ");
        TString title = TString::Format("%s (Small Angle: 0-4)", fname.Data());
        hist->SetTitle(title);
        c1->Print(pdfName, "pdf");
        f->Close();
    }
    
    // Plot large angle files
    for (const auto &fname : largeAngleFiles) {
        TFile *f = TFile::Open(fname);
        if (!f || f->IsZombie()) continue;
        
        TH2D *hist = (TH2D*)f->Get("dh2D_xz_det_3");
        if (!hist) {
            f->Close();
            continue;
        }
        
        c1->Clear();
        hist->Draw("COLZ");
        TString title = TString::Format("%s (Large Angle: 6-15)", fname.Data());
        hist->SetTitle(title);
        c1->Print(pdfName, "pdf");
        f->Close();
    }
    
    // Close PDF
    c1->Print(pdfName + "]", "pdf");
    
    cout << "PDF created: " << pdfName << endl;
}