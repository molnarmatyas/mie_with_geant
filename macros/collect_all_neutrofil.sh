#!/bin/bash

# Run txtToHist.C for first neutrofil output
root -b -q txtToHist.C\(\"outputneutrofil_rad.txt6500_990_0\"\)
mv figs/100M_3D_modell_cell_in_cell_beamprofiler_source_all_sensors_neutrofil_13um_normalized_background_removed_dh2D_xy_bw.png figs/1_full.png && mv figs/100M_3D_modell_cell_in_cell_beamprofiler_source_all_sensors_neutrofil_13um_normalized_background_removed_dh2D_xy_bw.png figs/1_full.png
mv figs/100M_3D_modell_cell_in_cell_beamprofiler_source_all_sensors_neutrofil_13um_normalized_background_removed_dh2D_xy_bw_bg_removed.png figs/1_nobg.png
# Loop to run txtToHist.C for all neutrofil outputs
for i in {0..18}
do
  root -b -q txtToHist.C\(\"outputneutrofil_rad.txt6500_990_0_$i\"\)
  mv figs/100M_3D_modell_cell_in_cell_beamprofiler_source_all_sensors_neutrofil_13um_normalized_background_removed_dh2D_xy_bw.png figs/$((i+2))_full.png
  mv figs/100M_3D_modell_cell_in_cell_beamprofiler_source_all_sensors_neutrofil_13um_normalized_background_removed_dh2D_xy_bw_bg_removed.png figs/$((i+2))_nobg.png
done
echo "All neutrofil histograms created and moved to figs/ directory."