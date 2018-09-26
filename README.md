README.md

Footprint calculation for OzFlux netCDF data files, generally L3.

Calculation of Kormann and Meixner, 2001 and Kljun et al., 2015 footprint climatologies.
Additionally a windrose climatology can be plotted for the same time periods.

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

You need:
1) footprint - controlfile
2) windrose  - controlfile (optional, you can run windrose using the footprint controlfile!)

For a super explanation how to install python, git, and all necessary tools please refer to 
https://github.com/OzFlux/PyFluxPro by Peter Isaac

How to use it: 
Just run python footprint_GUI.py in your FootPrint directory. 
Setup your own controlfile. For FootPrints you need the height of your measurements and an average height 
of the canopy. You also need to specify how large the area of the footprint may be and give 
the number of cells for your grid. The larger the size and the larger the cell number is the 
longer it takes to calculate the footprint. 
The size of the footprint area depends on the atmospheric conditions, at night times the 
footprint is generally larger than over the day when turbulence is well developed. For measurements 
of about 5 to 10 m above the displacement height a footprint area of 500 by 500 m2 is generally
sufficient. However, test it on a small protion of your data set before you do long timeseries.
Climatologies are done for either "daily", "monthly" and "yearly", additionally you can specify
a single time or a time period setting Start- and End-time.
Note: Do daily only on a short timeseries. It creates a file for every single day (365 per year).
Plotting windroses is simple, you can use the same controlfile as you did for footprints (all info is included).

Please let me know of any problems you encounter. It works best for L3 netCDF files using the 
standard OzFlux names. 
