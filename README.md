# Gratton2020_NI_HFmotion
The included MATLAB scripts were used in the study: Gratton et al. (2020) Removal of high frequency contamination from motion estimates in single-band fMRI saves data without biasing functional connectivty. Neuroimage.

These scripts will allow users to make measurements of the frequency content of motion metrics, relate them to demographic characteristics, and use different motion metrics to detect motion-functional connectivity links following Ciric et al., 2017.

The scripts should be run in the following order:
1. motion_HFcharacteristics.m
2. benchmarking_analysis.m

These scripts require some supporting files that can be found on Matlab Central. In addition, the scripts assume you have already have:
1. The 6 motion paramaters for each subject (in a specified order and provided in mm - see heading)
2. Subject demographic values in a structure (this can be altered as needed)
3. Subject correlation matrices and a distance matrix specifying distances between regions

The code is provided "as is", without any warranty. We can not provide techinical support for the code. For general questions or bugs found, please email cgratton@northwestern.edu. Thank you!
