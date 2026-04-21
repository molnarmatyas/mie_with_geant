void plot_det3_beamprofiler(bool degreeByDegree = true)
{
    // ROOT macro to plot BeamProfiler images from output files
    // Plots dh2D_xz_det_3 histograms with COLZ option
    // Saves all plots to a single PDF file
    //
    // Parameters:
    //   degreeByDegree = true : plot individual degree-by-degree files (0-1, 1-2, ..., 14-15)
    //   degreeByDegree = false: plot 0-4 and 6-15 aggregate files (legacy mode)
    
    gStyle->SetOptStat(0);
    
    TString outDir = gSystem->WorkingDirectory();
    TSystemDirectory dir(".", outDir);
    TList *files = dir.GetListOfFiles();
    
    // Collect output files
    vector<TString> smallAngleFiles;  // 0-4 degrees
    vector<TString> largeAngleFiles;  // 5+ degrees

    auto extractAngleRange = [](const TString &fullPath, int &startAngle, int &endAngle) {
        TString name = gSystem->BaseName(fullPath);

        // New/clean style: output_X-Y.root
        if (sscanf(name.Data(), "output_%d-%d.root", &startAngle, &endAngle) == 2) {
            return true;
        }

        // txtToHist-style long prefix: ..._output_X-Y_0_alpha_output.root
        const Ssiz_t markerPos = name.Index("_output_");
        if (markerPos != kNPOS) {
            TString tail = name(markerPos + 8, name.Length() - (markerPos + 8));
            if (sscanf(tail.Data(), "%d-%d", &startAngle, &endAngle) == 2) {
                return true;
            }
        }

        return false;
    };
    
    TIter next(files);
    TSystemFile *file;
    while ((file = (TSystemFile*)next())) {
        TString fname = file->GetName();
        
        if (degreeByDegree) {
            if (fname.EndsWith(".root")) {
                int startAngle = -1;
                int endAngle = -1;
                if (extractAngleRange(fname, startAngle, endAngle)) {
                    if (startAngle < 5) {
                        smallAngleFiles.push_back(outDir + "/" + fname);
                    } else {
                        largeAngleFiles.push_back(outDir + "/" + fname);
                    }
                }
            }
        } else {
            // Legacy mode: look for outputtestdhkm files with 0-4 or 6-15
            if (fname.Contains("outputtestdhkm") && fname.EndsWith(".root")) {
                if (fname.Contains("0-4")) {
                    smallAngleFiles.push_back(outDir + "/" + fname);
                } else if (fname.Contains("6-15")) {
                    largeAngleFiles.push_back(outDir + "/" + fname);
                }
            }
        }
    }
    
    // Sort files: numeric ordering for degree-by-degree mode, lexicographic for legacy mode.
    if (degreeByDegree) {
        auto degreeSort = [&](const TString &a, const TString &b) {
            int aStart = -1;
            int aEnd = -1;
            int bStart = -1;
            int bEnd = -1;

            const bool aOk = extractAngleRange(a, aStart, aEnd);
            const bool bOk = extractAngleRange(b, bStart, bEnd);

            if (aOk && bOk) {
                if (aStart != bStart) {
                    return aStart < bStart;
                }
                if (aEnd != bEnd) {
                    return aEnd < bEnd;
                }
            } else if (aOk != bOk) {
                // Parsed files first, unparsable ones last.
                return aOk;
            }
            return a < b;
        };

        sort(smallAngleFiles.begin(), smallAngleFiles.end(), degreeSort);
        sort(largeAngleFiles.begin(), largeAngleFiles.end(), degreeSort);
    } else {
        sort(smallAngleFiles.begin(), smallAngleFiles.end());
        sort(largeAngleFiles.begin(), largeAngleFiles.end());
    }
    
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
        cout << "Plotted (small) angle file: " << fname << endl;
        
        TString title;
        if (degreeByDegree) {
            int startAngle = -1;
            int endAngle = -1;
            if (extractAngleRange(fname, startAngle, endAngle)) {
                title = TString::Format("Degree Range %d-%d (Small Angle: 0-4 degrees)", startAngle, endAngle);
            } else {
                title = TString::Format("%s (Small Angle: 0-4 degrees)", fname.Data());
            }
        } else {
            title = TString::Format("%s (Small Angle: 0-4 degrees)", fname.Data());
        }
        
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
        cout << "Plotted (large) angle file: " << fname << endl;
        
        TString title;
        if (degreeByDegree) {
            int startAngle = -1;
            int endAngle = -1;
            if (extractAngleRange(fname, startAngle, endAngle)) {
                title = TString::Format("Degree Range %d-%d (Large Angle: 5-15 degrees)", startAngle, endAngle);
            } else {
                title = TString::Format("%s (Large Angle: 5-15 degrees)", fname.Data());
            }
        } else {
            title = TString::Format("%s (Large Angle: 6-15)", fname.Data());
        }
        
        hist->SetTitle(title);
        c1->Print(pdfName, "pdf");
        f->Close();
    }
    
    // Close PDF
    c1->Print(pdfName + "]", "pdf");
    
    cout << "PDF created: " << pdfName << endl;
}