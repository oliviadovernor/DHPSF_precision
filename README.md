# 3D DH Precision 
 
 1. Run peakFit automation script = PeakFit_Automation.ijm on selected beads. Run macro in ImageJ - change fitting parameters accordingly. 
 2. Run excel files through MATLAB code to extract z position  from xy. Script saves .3d and .csv output files in the same parent folder.
 3. Run the main.py script which calculates the precision of DH PSF signals as the error between a fitted 2D guassian for each axis and the data outputted from MATLAB. Results can be shown visulaised by a histogram which includes 2 standard deviations either side of the mean.

Automations have been completed for PeakFit and MATLAB scripts. 
Both run through subfolders within a parent folder so that a complete dataset can be analysed at any one time. 