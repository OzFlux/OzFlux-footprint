README.md

Footprint calculation for OzFlux netCDF data files, generally L3.

Calculation of Kormann and Meixner, 2001 and Kljun et al., 2015 footprint climatologies.
Additionally a windrose climatology can be plotted for the same time periods.
============================================================================================================
Kljun et al., 2015:
Kljun, N., Calanca, P., Rotach, M. W., and Schmid, H. P., 2015: A simple two-dimensional parameterisation
  for Flux Footprint Prediction (FFP), Geosci. Model Dev., 8, 3695-3713.
  For details on the script see Natascha Kljun's website: http://footprint.kljun.net/
Kormann and Meixner, 2001:
Neftel et al. (2008) ART_footprint tool: https://zenodo.org/record/816236#.W2eqUXBx3VM 
  (http://doi.org/10.5281/zenodo.816236), which 
  Kormann, R. and Meixner, F.X., 2001: An analytical footprint model for non-neutral stratification.
    Boundary-Layer Meteorology 99: 207. https://doi.org/10.1023/A:1018991015119 and 
    Neftel, A., Spirig, C., Ammann, C., 2008: Application and test of a simple tool for operational footprint 
    evaluations. Environmental Pollution 152, 644-652.
Windrose:
  uses Lionel Roubeyrie's windrose.py V1.4. It works, so why should I change it.
  Newer version available at https://zenodo.org/record/1406384#.W6r6D3Bx3VM
============================================================================================================

footprint - controlfile:
========================================================================================================
 To convert this template to a version that works for data output by the standard OzFlux
 data logger program:
  1) Replace the following place holders with values for your site:
     <site_name> - the name of the site
     <year>      - the year of the data

[Files] - contains paths, file names and location of data in file
[Options] - contains the climatology to be calculated.
========================================================================================================
[Files]
    file_path = ../Sites/<site_name>/Data/Processed/<year>/
    in_filename  = <site_name>_<year>_L3.nc
    plot_path = ../Sites/<site_name>/Plots/
[Options]
    Climatology = 'Special' #'Annual' #'Monthly' #'Daily' #'Single'
    #StartDate = '2009-02-07 00:30'
    #EndDate   = '2009-02-08 00:00'

windrose - controlfile (optional, you can run windrose using the footprint controlfile!):
[Files]
    file_path = ../Sites/<site_name>/Data/Processed/<year>/
    in_filename  = <site_name>_<year>_L3.nc
    plot_path = ../Sites/<site_name>/Plots/

[Options]
    Climatology = 'Special' #'Annual' #'Monthly' #'Daily' #'Single'
    #StartDate = '2009-02-07 00:30'
    #EndDate   = '2009-02-08 00:00'

For a super explanation how to install python, git, and all necessary tools please refer to 
https://github.com/OzFlux/PyFluxPro by Peter Isaac