/////////////////////////////
// Auto PeakFit for Olivia //
/////////////////////////////



// Specify here the relevant times... i.e. for (i=0; i<5; i++) means time_0, time_1... time_5
for (i=1; i<2; i++) {
path = "C:/Users/ojd34/OneDrive - University of Cambridge/Deskto/static_precision_MLE_fit/_" + i + "/";
output = path + "Results/";
File.makeDirectory(output);



// Specify here the numbering for the surface_N files to be PeakFitted... i.e. for (j=1; j<101; j++) will fit surface_1 to surface_100
for (j=1; j<2; j++) {
fileName = "numPhotons_" + j + ".tif";
open(path + fileName);
selectWindow(fileName);
run("Peak Fit", "template=[None] camera_type=EMCCD calibration=266 exposure_time=30 psf=[Circular Gaussian 2D] spot_filter_type=Difference spot_filter=Mean smoothing=0.32 search_width=0.90 border_width=0.88 fitting_width=2.14 fit_solver=MLE fail_limit=3 pass_rate=0.50 neighbour_height=0.30 residuals_threshold=1 duplicate_distance=0.50 shift_factor=1.20 signal_strength=3 min_photons=100 min_width_factor=0.40 width_factor=1.20 precision=70 show_results_table image=[Localisations (width=precision)] results_format=Text results_directory=[C:/Users/ojd34/OneDrive - University of Cambridge/Desktop/static_precision_MLE_fit/1] save_to_memory camera_bias=400 gain=5.9000 read_noise=7.1600 psf_parameter_1=1.700 precision_method=Mortensen table_distance_unit=[pixel (px)] table_intensity_unit=photon table_angle_unit=[degree (°)] table_precision=0 image_size_mode=Scaled image_scale=10 image_size=0 image_pixel_size=5 lut=Red-Hot file_distance_unit=[pixel (px)] file_intensity_unit=photon file_angle_unit=[degree (°)] spot_filter2=Mean smoothing2=2 camera_bias=400 read_noise=7.16 quantum_efficiency=0.95 em-ccd search_method=Powell relative_threshold=1.0E-6 absolute_threshold=1.0E-16 max_iterations=20 max_function_evaluations=2000 stack");
  <<warning: the options string contains one or more non-ascii characters>>

run("Close All"); }