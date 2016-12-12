# raman-fitting
A script for automatically analysing Raman spectra of carbonaceous material

This script implements the procedure for fitting Raman Spectra as described in: 

Sparkes et al., 2013. Automated Analysis of Carbon in Powdered Geological and Environmental Samples by Raman Spectroscopy. Applied Spectroscopy, 67, 7, 779-788. DOI: 10.1366/12-06826. 

Publisher's version of the paper (subscription required): 
http://asp.sagepub.com/content/67/7/779

Final manuscript version archived free-of-charge according to publisher's open access policy:
http://e-space.mmu.ac.uk/613531/

This script analyses Raman spectra of Carbonaceous Material, by fitting Lorentzian distributions
to the G, D1, D2, D3 and D4 peaks, as well as correcting for a linear background.
The input is taken from a series of `x y` text files as produced by a Renishaw Raman spectrometer using
Wire software, or from other Raman spectrometry software. 

Proprietary .wxd files can be converted into two-column space-separated text files 
(wavenumber intensity) using the "Wire Batch Convert" program provided by Renishaw. 

The text files should be contained within one single folder, or grouped into sub-folders.

The script outputs three graphs, containing the raw spectra with linear background identified,
raw spectra with overalll fit superimposed and a residual shown, and the spectra following the fitting,
showing the fitted peaks after the background has been removed. The fitting parameters (peak locations,
amplitudes, widths and areas, as well as characteristic area ratios) are outputted to a summary file 
for further analysis

The script requires the following software to run:
 - A Unix / Linux environment (tested with Ubuntu)
 - Bash terminal program
 - Dos2unix text file conversion software
 - Gnuplot graphing software. 
     - Version 4.5 or above is required
 - Ghostscript PostScript and PDF manipulation software
 - The script "prepraman.sh" should be run before the first files in a given folder are analysed.
   This script creates folders and initiates some datafiles for the subsequent fits
 - Both this fitting script and "prepraman.sh" require permission to execute as programs

The script executes from the command line, in the form
~~~~
raman-fitting.sh [options] [input files]
~~~~

The options are:

`-q` Quiet mode - graphs appear on screen but immediately disappear

`-d` Delete - removes previous files from "acombinedresults.txt"

`-t[value]` Threshold - the signal-to-noise ratio below which a peak is too noisy to process


Input files can be listed individually, or selected all at once using a wildcard (e.g. \*.txt).
After analysis the results are written to a file entitled "acombinedresults.txt". Any filename
already in this file will be ignored and not re-fitted, hence the delete option.

Example code to prepare for and then analyse all samples with "taiwan" in the file name:
~~~~
prepraman.sh 
raman-fitting.sh -d -q -t 5 taiwan*.txt
~~~~

Note: The included script "cropraman.sh" will take a file and crop to certain wavenumbers. The
fitting procedure is less accurate if files extend too far beyond 2000 cm-1 as the assumption
of a linear background is no longer valid.
