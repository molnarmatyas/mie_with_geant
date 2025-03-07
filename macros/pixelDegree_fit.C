
TH2D* dhR_alpha; 
TGraph* dhR_alpha_max;
void pixelDegree_fit()
{ 
  TFile *hist_infile;
  hist_infile = new TFile("100M_3D_complete_model_shifted_cell_polistirol_15um__alpha_output.root","read");
  std::string outputprefix = "100M_3D_complete_model_shifted_cell_polistirol_15um_point_source_";
  TFile* outfile = new TFile(Form("%s_ralpha_output.root", outputprefix.c_str()), "RECREATE");

  
  dhR_alpha = (TH2D*)hist_infile->Get("dhR_alpha");
  int nbins = dhR_alpha->GetXaxis()->GetNbins();
  dhR_alpha_max = new TGraph(nbins);
  dhR_alpha_max->SetName("dhR_alpha_max");
  Int_t binx = 0;
  for(int ibin = 0; ibin < nbins; ibin++)
  { 
    double deg = dhR_alpha->GetXaxis()->GetBinCenter(ibin);
    TH1D* projection = dhR_alpha->ProjectionY("_py",ibin, ibin);
    std::cout << projection->GetMaximum() << " " << projection->GetBinWithContent(projection->GetMaximum(), binx) << " " << binx << " "<< projection->GetBinCenter(binx) << std::endl;

    projection->GetBinWithContent(projection->GetMaximum(), binx);
    double max_hit = projection->GetBinCenter(binx);
    if(max_hit < 1.0) continue;
    std::cout << "deg: " << deg << " max hit " << max_hit << std::endl;
    dhR_alpha_max->SetPoint(ibin, deg, max_hit);
  }
  dhR_alpha_max->Write();
  outfile->Close();

}
