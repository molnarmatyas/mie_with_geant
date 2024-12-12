#include <fstream>
#include <sstream>
#include <string>
#include <TH2D.h>
#include <TFile.h>

void txtToHist() {
    // Input text file path
    std::string filename = "../build/output3000_990_0_0.txt";
    
    // Create a 2D histogram (you may want to adjust binning)
    TH2D* hist = new TH2D("hist", "Histogram from Text File", 
                           10000, -0.2, 0.2, // x-axis bins, min, max 
                           10000, -0.2, 0.2  // y-axis bins, min, max
    );
    
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
        double val1, val2, val3, val4, val5, val6;
        
        // Read all 6 columns
        if (!(iss >> val1 >> val2 >> val3 >> val4 >> val5 >> val6)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue;
        }
        
        // Fill the histogram with the last two columns
        hist->Fill(val5, val6);
    }
    
    // Create a ROOT output file to save the histogram
    TFile* outfile = new TFile("output.root", "RECREATE");
    hist->Write();
    outfile->Close();
    
    // Optional: draw and save as image
    TCanvas* c1 = new TCanvas("c1", "Histogram", 1920, 1080);
    c1->SetLogz(1);
    hist->Draw("COLZ");
    c1->SaveAs("output3000_um_pointsource.pdf");
}
