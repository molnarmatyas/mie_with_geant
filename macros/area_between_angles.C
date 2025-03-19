
TH2D* dhphi_alpha_det; 
void area_between_angles()
{ 
  TFile *hist_infile;
  hist_infile = new TFile("1M_3D_updated_atan2_swapped_phi_complete_model_polistirol_15um__alpha_output.root","read");

  
  dhphi_alpha_det = (TH2D*)hist_infile->Get("dhphi_alpha_det_2");
  int lower_angle_1 = 3;
  int upper_angle_1 = 6;
  int lower_angle_2 = 6;
  int upper_angle_2 = 10;
  double intensity_1 = dhphi_alpha_det->Integral(0, dhphi_alpha_det->GetXaxis()->GetNbins(),
                                                 dhphi_alpha_det->GetYaxis()->FindBin(lower_angle_1),
                                                 dhphi_alpha_det->GetYaxis()->FindBin(upper_angle_1));
  double intensity_2 = dhphi_alpha_det->Integral(0, dhphi_alpha_det->GetXaxis()->GetNbins(),
                                                 dhphi_alpha_det->GetYaxis()->FindBin(lower_angle_2),
                                                 dhphi_alpha_det->GetYaxis()->FindBin(upper_angle_2));
  double intensity_full = dhphi_alpha_det->Integral();

  std::cout << "Intensity in region " << lower_angle_1 << "-" << upper_angle_1 << " : " << intensity_1 << std::endl;
  std::cout << "Intensity in region " << lower_angle_2 << "-" << upper_angle_2 << " : " << intensity_2 << std::endl;
  std::cout << "Intensity in total: " << intensity_full << std::endl;

  std::cout << lower_angle_1 << "-" << upper_angle_1 << " / " << lower_angle_2 << "-" << upper_angle_2 << " = "<< intensity_1 / intensity_2 << std::endl;
  std::cout << lower_angle_1 << "-" << upper_angle_1 << " / full = "<< intensity_1 / intensity_full << std::endl;
}
