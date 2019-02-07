# standard modules
import datetime
import logging
import os
import zipfile
# 3rd party modules
import matplotlib.pyplot as plt
import numpy
import pandas
#from scipy.misc.pilutil import imread
# PFP modules
import constants as c
import footprint_io
import footprint_utils

# Kljun, N., Calanca, P., Rotach, M. W., and Schmid, H. P., 2015:
# A simple two-dimensional parameterisation for Flux Footprint Prediction (FFP),
# Geosci. Model Dev., 8, 3695-3713.
import footprint_FFP_climatology as calcfootNK
# The following script cacl_footprint_FKM_climatology is based on the Neftel et al. (2008) ART_footprint tool:
# https://zenodo.org/record/816236#.W2eqUXBx3VM (http://doi.org/10.5281/zenodo.816236), which 
# Kormann, R. and Meixner, F.X., 2001: An analytical footprint model for non-neutral stratification.
# Boundary-Layer Meteorology 99: 207. https://doi.org/10.1023/A:1018991015119 and Neftel, A., Spirig, C.,
# Ammann, C., 2008: Application and test of a simple tool for operational footprint evaluations. Environmental
# Pollution 152, 644-652.
import footprint_FKM_climatology as calcfootKM

logger = logging.getLogger("footprint_log")

# constant for converting degrees to radiant
c_d2r = numpy.pi/180.0
# constant to convert the footprint area from x,y in lon,lat coordinates at the tower site
onedegree = 6378100.0 * c_d2r # distance in m for 1 deg of latitude

# Coordinate steps in footprint process
def footprint_main(cf, mode):
    """
    This script reads data from a PyFluxPro .nc file and processes the data for:
    (1) Kormann&Meixner uses input (single time)        zm,z0,ustar,umean,L,sigmav,wind_dir
    (2) Kljun et al. uses input (vectors for all times) zm,z0,ustar,umean,L,sigmav,wind_dir,Habl
        so Natascha's FFP also needs the height of the boundary layer ===> currently ERAI data got Habl,
        and ACCESS got Habl00 ... Habl22
    === > input for Kormann & Meixner and Natascha Kljun's footprint climatology
    (a) PyFluxPro L3 netcdf file
    (b) ERAI/ACCESS netcdf file
    === > output for the climatology
    (a) daily footprint climatology
    (b) monthly footprint climatology
    (c) annual footprint climatology
    (d) special time set in controlfile for footprint climatology
    GOAL: Footprint climatology can be done on a set time in controlfile
          calculating Kljun et al., 2015 and Kormann and Meixner, 2001 footprint
    DONE: set time in controlfile, special, daily, monthly and annual
          Kljun et al. (2015) footprint
          Kormann and Meixner (2001) footprint
          save footprint fields in netcdf file
    Still to do: calculate Habl if not exist, better is set Habl (latter is done)
    C.M.Ewenz, 10 Jun 2018
               21 Jun 2018 (corrections to monthly indexing)
               29 Jun 2018 (kml file, single time stamp)
    P.R.Isaac,    Jul 2018 (re-wrote fp_data_in to get_footprint_data_in; configuration in get_footprint_cfg; time slicing; etc)
    C.M.Ewenz, 30 Jul 2018 (cleaned up printing of info, warning and error messages - include messages in logger)
    C.M.Ewenz, 22 Jan 2019 (included "Hourly" for plotting every timestep)
    C.M.Ewenz, 08 Feb 2019 (estimate cumulative footprint field)
    """
    logger.info(' Read input data files ...')
    # get the L3 data
    ds = get_footprint_data_in(cf, mode)
    ldt = ds.series["DateTime"]["Data"]
    # get the configuration data for the footprint
    d = get_footprint_cfg(cf, ds)
    logger.info(' Starting footprint calculation ...')
    list_StDate, list_EnDate = footprint_utils.create_index_list(cf, d, ldt)
    logger.info(' Starting footprint climatology calculation ...')
    # !!! Prepare Output netcdf file !!!
    # Set initial x,y Variables for output
    xout = numpy.linspace(d["xmin"], d["xmax"], d["nx"] + 1)
    yout = numpy.linspace(d["ymin"], d["ymax"], d["nx"] + 1)
    lat0 = float(d["latitude"])
    lon0 = float(d["longitude"])
    lat = lat0 + (yout / onedegree)
    lon = lon0 + (xout / (numpy.cos(lat0*c_d2r) * onedegree))
    lon_2d, lat_2d = numpy.meshgrid(lon, lat)
    # - Initialise output netcdf file and write x,y grid into file as xDistance and yDistance from the tower
    nc_name = d["file_out"]
    #print 'nc_name = ',nc_name
    nc_file = footprint_io.nc_open_write(nc_name)
    # create the x and y dimensions.
    nc_file.createDimension('longitude', len(lon))
    nc_file.createDimension('latitude', len(lat))
    # create time dimension (record, or unlimited dimension)
    nc_file.createDimension('time', None)
    # create number of footprints in climatology dimension (record, or unlimited dimension)
    nc_file.createDimension('dtime', None)
    nc_file.createDimension('num', None)
    # Define coordinate variables, which will hold the coordinate information, x and y distance from the tower location.
    X = nc_file.createVariable('longitude', "d", ('longitude',))
    Y = nc_file.createVariable('latitude', "d", ('latitude',))
    # Define time variable and number of footprints variable at each time
    tx = nc_file.createVariable('dtime', "d", ('dtime',))
    num = nc_file.createVariable('num', "d", ('num',))
    # Assign units attributes to coordinate var data, attaches text attribute to coordinate variables, containing units.
    X.units = 'degree'
    Y.units = 'degree'
    # write data to coordinate vars.
    X[:] = lon
    Y[:] = lat
    # create the sumphi variable
    phi = nc_file.createVariable('sumphi', "d", ('time', 'longitude', 'latitude'))
    # set the units attribute.
    phi.units = ' '
    # === General inputs for FFP
    zmt = d["zm_d"]
    domaint = [d["xmin"], d["xmax"], d["ymin"], d["ymax"]]
    nxt = d["nx"]
    rst = None #[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] #[90.] #None #[20.,40.,60.,80.]
    # if plotting to screen      is requested then iplot = 1
    # if plotting to googleEarth is requested then iplot = 2
    iplot = int(cf['General']['iplot'])
    # IF export images in kml format?
    if iplot == 2:  # write kml - format header
        kmlname = d["site_name"] + '_' + mode + '_fp' + '.kml'
        kml_name_path = d["plot_path"] +kmlname
        fi = open(kml_name_path, 'w')
        kml_initialise(d,fi,mode)
        # ------------------------
    # After deciding which climatology is done, let's do it!
    irun = -1
    for i in range(0, len(list_StDate)):
        irun = irun+1
        # get the start and end indices
        si = list_StDate[i]
        ei = list_EnDate[i]
        # get the series as masked arrays
        umeant, _, _ = footprint_utils.GetSeriesasMA(ds, "Ws", si=si, ei=ei)
        olt, _, _ = footprint_utils.GetSeriesasMA(ds, "L", si=si, ei=ei)
        sigmavt, _, _ = footprint_utils.GetSeriesasMA(ds, "V_Sd", si=si, ei=ei)
        ustart, _, _ = footprint_utils.GetSeriesasMA(ds, "ustar", si=si, ei=ei)
        wind_dirt, _, _ = footprint_utils.GetSeriesasMA(ds, "Wd", si=si, ei=ei)
        z0t, _, _ = footprint_utils.GetSeriesasMA(ds, "z0", si=si, ei=ei)
        ht, _, _ = footprint_utils.GetSeriesasMA(ds, "Habl", si=si, ei=ei)
        # get a composite mask over all variables
        mask_all = numpy.ma.getmaskarray(ustart)
        for item in [umeant, olt, sigmavt, wind_dirt, z0t, ht]:
            mask_item = numpy.ma.getmaskarray(item)
            mask_all = numpy.ma.mask_or(mask_all, mask_item)
        # and then apply the composite mask to all variables and remove masked elements
        umeant = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, umeant)))
        olt = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, olt)))
        sigmavt = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, sigmavt)))
        ustart = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, ustart)))
        wind_dirt = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, wind_dirt)))
        z0t = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, z0t)))
        ht = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, ht)))
        if len(umeant) == 0:
            msg = "No footprint input data for "+str(ldt[si])+" to "+str(ldt[ei])
            logger.warning(msg)
            num[irun]=0
        else:
            if mode == "kljun":
                FFP = calcfootNK.FFP_climatology (zm=zmt,z0=z0t,umean=umeant,h=ht,ol=olt,sigmav=sigmavt,ustar=ustart,\
                                                  wind_dir=wind_dirt,domain=domaint,dx=None,dy=None,nx=nxt,ny=None,\
                                                rs=rst,rslayer=1,smooth_data=1,crop=False,pulse=None,verbosity=2)
                x              = FFP['x_2d']
                y              = FFP['y_2d']
                f              = FFP['fclim_2d']
                num[irun]      = FFP['n']
                #tx[irun] = str(ldt[ei])
                phi[irun,:,:] = f
                fmax=numpy.amax(f)
                f=f/fmax
            elif mode == "kormei":
                FKM = calcfootKM.FKM_climatology(zm=zmt, z0=z0t, umean=umeant, ol=olt, sigmav=sigmavt, ustar=ustart,\
                                                 wind_dir=wind_dirt, domain=domaint, dx=None, dy=None, nx=nxt, ny=None, \
                         rs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], rslayer=0,\
                         smooth_data=1, crop=False, pulse=None, verbosity=2)
                x              = FKM['x_2d']
                y              = FKM['y_2d']
                f              = FKM['fclim_2d']
                num[irun]      = FKM['n']
                #tx[irun] = str(ldt[ei])
                phi[irun,:,:] = f
                fmax=numpy.amax(f)
                f=f/fmax
            else:
                msg = " Unrecognised footprint type " + str(mode)
                logger.error(msg)
                return
            i_cum = footprint_utils.get_keyvaluefromcf(cf,["Options"],'Cumulative')
            if i_cum:
               msg = "Caclulated cumulative footprint field"
               logger.info(msg)
               f_min = 0.05
               f_step = 0.05
               f = calc_cumulative(f,f_min,f_step)
            
        # ====================================================================================================
        # get the default plot width and height
        #clevs = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
        clevs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
        imagename = footprint_utils.get_keyvaluefromcf(cf,["General"],'OzFlux_area_image')
        if not num[irun] == 0:
            if iplot == 1:
                # plot on screen and in jpg
                plotphifield(x, y, ldt[si], ldt[ei], f, d["site_name"], mode, clevs, imagename,i_cum)
            elif iplot == 2:
                # plot on screen, in jpg and write kml (google earth) file
                #plotphifield(x, y, ldt[si], ldt[ei], f, d["site_name"], mode, clevs, imagename)
                kml_write(lon, lat, ldt[si], ldt[ei], f, d["site_name"], mode, clevs, fi, d["plot_path"],i_cum)
            plot_num = plt.gcf().number
            if  plot_num > 20:
                plt.close("all")
        # ====================================================================================================
        # Some stats:
        #  a) Possible total number of footprints per climatology = ei - si
        tot_fp = ei - si
        #  b) Remove each time step with "no value", number of times footprint is run = len(umeant)
        tot_fp_nv = len(umeant)
        #  c) Fianl number of valid footprints, removed all time steps where conditions not valid
        tot_valid = num[irun]
        msg = 'Total = ' + str(tot_fp) + ' Used = ' + str(tot_fp_nv) + ' Valid = ' + str(tot_valid) + ' footprints!'
        logger.info(msg)
        
        #progress = float(i+1)/float(len(list_StDate))
        #footprint_utils.update_progress(progress) 
    if iplot == 2:
        # Finish kml file and process a compressed kmz file including all images
        kml_finalise(d,fi,mode,kmlname)
    # ================================================================
    msg = " Finished " + str(mode) + " footprint writing"
    logger.info(msg)
    msg = " Closing netcdf file " + str(nc_name)
    logger.info(msg)
    nc_file.close()
    # ================================================================

def get_footprint_cfg(cf, ds):
    # Build dictionary of additional configs
    d={}
    # === Which climatology, either definded time, daily, monthly or annual
    d["Climatology"] = footprint_utils.get_keyvaluefromcf(cf,["Options"],"Climatology",default="Special")
    climfreq = d["Climatology"]
    #
    if "out_filename" in cf['Files']:
        file_out = os.path.join(cf['Files']['file_path'],cf['Files']['out_filename'])
    else:
        if climfreq == 'Annual':
            file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_y_fp.nc"))
        elif climfreq == 'Monthly':
            file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_m_fp.nc"))
        elif climfreq == 'Daily':
            file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_d_fp.nc"))
        elif climfreq == 'Hourly':
            file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_h_fp.nc"))
        elif climfreq == 'Single':
            file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_si_fp.nc"))
        elif climfreq == 'Special':
            file_out = os.path.join(cf['Files']['file_path'],cf['Files']['in_filename'].replace(".nc","_sp_fp.nc"))

    plot_path = "plots/"
    if "plot_path" in cf["Files"]: plot_path = os.path.join(cf["Files"]["plot_path"],"FP/")
    if not os.path.isdir(plot_path): os.makedirs(plot_path)

    results_path = cf['Files']['file_path']
    if not os.path.isdir(results_path): os.makedirs(results_path)

    d["tower_height"]   = float(cf["Tower"]["tower_height"])
    d["canopy_height"]  = float(cf["Tower"]["canopy_height"])
    d["footprint_size"] = int(cf["Tower"]["footprint_size"])
    d["zm_d"]           = d["tower_height"]-(2.0/3.0*d["canopy_height"])
    d["xTower"]         = 0 #int(cf['Tower']['xTower'])
    d["yTower"]         = 0 #int(cf['Tower']['yTower'])
    d["xmin"]           = -0.5*d["footprint_size"]
    d["xmax"]           =  0.5*d["footprint_size"]
    d["ymin"]           = -0.5*d["footprint_size"]
    d["ymax"]           =  0.5*d["footprint_size"]
    d["nx"]             = int(cf["Tower"]["num_cells"])

    d["flux_period"] = int(ds.globalattributes["time_step"])
    #d["timezone"] = int(ds.globalattributes["timezone"])
    d["site_name"] = ds.globalattributes["site_name"]
    if "Latitude" in cf["Tower"]:
        d["latitude"] = cf["Tower"]["Latitude"]
    else:
        d["latitude"] = ds.globalattributes["latitude"]
    if "Longitude" in cf["Tower"]:
        d["longitude"] = cf["Tower"]["Longitude"]
    else:
        d["longitude"] = ds.globalattributes["longitude"]

    d["call_mode"] = footprint_utils.get_keyvaluefromcf(cf,["Options"],"call_mode",default="interactive",mode="quiet")
    d["show_plots"] = footprint_utils.get_keyvaluefromcf(cf,["Options"],"show_plots",default=True,mode="quiet")
    d["file_out"] = file_out
    d["plot_path"] = plot_path

    return d

def get_footprint_data_in(cf, mode):
    import footprint_utils
    # read input data and prepare for input into Kormann and Meixner, 2001 or Kljun et al., 2015
    # python routines
    # ---------------------- Get input / output file name ------------------------------------
    # Set input file and output path and create directories for plots and results
    file_in = os.path.join(cf['Files']['file_path'], cf['Files']['in_filename'])
    # read the netcdf file
    msg = ' Reading netCDF file ' + str(file_in)
    logger.info(msg)
    ds = footprint_io.nc_read_series(file_in)
    nrecs = int(ds.globalattributes["nc_nrecs"])
    # array of 0s for QC flag
    f0 = numpy.zeros(nrecs, dtype=numpy.int32)
    # array of 1s for QC flag
    f1 = numpy.ones(nrecs, dtype=numpy.int32)
    # read the external file for Habl if mode = kljun
    # === check to see if we have Habl timeseries in imports ??? What if not? Botheration!
    if mode == "kljun":
        footprint_io.ImportSeries(cf, ds)
    else: # kormei does not need Habl
        Habl = footprint_utils.create_empty_variable("Habl", nrecs)
        Habl["Label"] = "Habl"
        Habl["Data"] = numpy.ma.array(numpy.full(nrecs, 1000))
        Habl["Flag"] = f0
        Habl["Attr"] = {"long_name":" Boundary-layer height", "units":"m",
                        "standard_name":"not defined"}
        footprint_utils.CreateVariable(ds, Habl)
    # check to see if Monin-Obukhov length is in the data structure
    if "L" not in ds.series.keys():
        # if not, calculate it
        footprint_utils.CalculateMoninObukhovLength(ds)
    # if the cross wind standard deviation is not in the data set (quite common) then use something else
    if "V_Sd" not in ds.series.keys():
        # could do better with:
        # 1) reprocess L3 and output variance of U, V and W
        # 2) estimated from standard deviation of wind direction (if available)
        # 3) estimate using MO relations (needs Habl)
        V_Sd = footprint_utils.create_empty_variable("V_Sd", nrecs)
        if "Uy_Sd" in ds.series.keys():
            #logger.warning("Stdev of cross wind component not in data structure, estimated as 0.5*Ws")
            Uy_Sd = footprint_utils.GetVariable(ds, "Uy_Sd")
            V_Sd["Data"] = Uy_Sd["Data"]
            V_Sd["Attr"]["height"] = Uy_Sd["Attr"]["height"]
        else:
            logger.warning("Stdev of cross wind component not in data structure, estimated as 0.5*Ws")
            Ws = footprint_utils.GetVariable(ds, "Ws")
            V_Sd["Data"] = 0.5*Ws["Data"]
            V_Sd["Attr"]["height"] = Ws["Attr"]["height"]
        V_Sd["Flag"] = numpy.where(numpy.ma.getmaskarray(V_Sd["Data"])==True, f1, f0)
        V_Sd["Attr"]["long_name"] = "Variance of cross-wind velocity component, estimated from Ws"
        V_Sd["Attr"]["units"] = "(m/s)2"
        footprint_utils.CreateVariable(ds, V_Sd)
    # === roughness length
    if "z0" not in ds.series.keys():
        z0 = footprint_utils.create_empty_variable("z0", nrecs)
        # check the global attriibutes first
        if "roughness_length" in ds.globalattributes.keys():
            roughness_length = float(ds.globalattributes["roughness_length"])
            z0["Data"] = numpy.ma.array(numpy.full(nrecs, roughness_length))
            z0["Attr"]["long_name"] = "Roughness length from global attributes"
        elif "roughness_length" in cf["Tower"]:
            roughness_length = float(cf["Tower"]["roughness_length"])
            z0["Data"] = numpy.ma.array(numpy.full(nrecs, roughness_length))
            z0["Attr"]["long_name"] = "Roughness length from footprint control file"
        else:
            zT = float(cf["Tower"]["tower_height"])
            zC = float(cf["Tower"]["canopy_height"])
            zm = zT-(2.0/3.0)*zC
            L = footprint_utils.GetVariable(ds, "L")
            ustar = footprint_utils.GetVariable(ds, "ustar")
            Ws = footprint_utils.GetVariable(ds, "Ws")
            z0["Data"] = footprint_utils.z0calc(zm, L["Data"], Ws["Data"], ustar["Data"])
            z0["Attr"]["long_name"] = "Roughness length calculated from u*, L, Ws and (z-d)"
        z0["Flag"] = numpy.where(numpy.ma.getmaskarray(z0["Data"])==True, f1, f0)
        z0["Attr"]["units"] = "m"
        footprint_utils.CreateVariable(ds, z0)
    return ds

def kml_initialise(d,fi,mode):
    # 
    #!#kmlname = d["site_name"] + '_' + mode + '_fp' + '.kml'
    #!#kml_name_path = d["plot_path"] +kmlname
    #!#fi = open(kml_name_path, 'w')
    fi.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    fi.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
    fi.write("<Folder>\n")
    fi.write("  <name>" + d["site_name"] + "</name>")
    # GE zooms in to the site location
    fi.write('  <LookAt>\n')
    fi.write('    <longitude>'+str(d["longitude"])+'</longitude>\n')
    fi.write('    <latitude>'+str(d["latitude"])+'</latitude>\n')
    fi.write('    <altitude>'+str(d["footprint_size"])+'</altitude>\n')
    fi.write('    <range>'+str(d["footprint_size"])+'</range>\n')
    fi.write('    <tilt>0</tilt>\n')
    fi.write('    <heading>0</heading>\n')
    fi.write('    <altitudeMode>relativeToGround</altitudeMode>\n')
    fi.write('  </LookAt>\n')
    # Define the legend in a screen overlay
    fi.write('  <ScreenOverlay>\n')
    fi.write('    <name>Legend: Footprint</name>\n')
    fi.write('    <Icon> <href>cbar.png</href></Icon>\n')
    fi.write('    <overlayXY x="0" y="0" xunits="fraction" yunits="fraction"/>\n')
    fi.write('    <screenXY x="25" y="95" xunits="pixels" yunits="pixels"/>\n')
    fi.write('    <rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>\n')
    fi.write('    <size x="0" y="0" xunits="pixels" yunits="pixels"/>\n')
    fi.write('  </ScreenOverlay>\n')
    # Adding our own icon for the placemark
    #fi.write('  <Style id="tower">\n')
    #fi.write('    <IconStyle>\n')
    #fi.write('      <scale>1.5</scale>\n')
    #fi.write('      <Icon>\n')
    #fi.write('        <href>ec_tower.png</href>\n') # !!! this file needs to be copied into the plot directory
    #fi.write('      </Icon>\n')
    #fi.write('    </IconStyle>\n')
    #fi.write('  </Style>\n')
    # Adding a placemark for the site
    fi.write('  <Placemark>\n')
    fi.write('      <name>'+ d["site_name"] +'</name>\n')
    #fi.write('      <styleUrl>#tower</styleUrl>')
    fi.write('      <Point>\n')
    fi.write('          <coordinates>'+str(d["longitude"])+','+str(d["latitude"])+',0</coordinates>\n')
    fi.write('      </Point>\n')
    fi.write('  </Placemark>\n')

def kml_write(lon, lat, zt1, zt2, data, station, mode, clevs, fi, plot_path,i_cum):
    plot_in='Footprint_'+ mode + zt1.strftime("%Y%m%d%H%M") +'.png'
    plotname=plot_path + plot_in
    width = 5
    height = width * data.shape[0]/data.shape[1]
    plt.ioff()
    plt.figure(figsize=(width,height))
    cs = plt.contourf(data,clevs,cmap=plt.get_cmap('hsv'),alpha=0.5)
    plt.axis('off')
    plt.savefig(plotname,transparent=True)
    #plt.clf()
    fn = plt.gcf().number
    plt.close(fn)
    # draw a new figure and replot the colorbar there
    fig,ax = plt.subplots(figsize=(width,height))
    cbar = plt.colorbar(cs,ax=ax)
    # =========================================================================
    #rlevs = [1 - clev for clev in clevs if clev is not None]
    #cbar.set_ticks(rlevs)
    cbar.set_ticks(clevs)
    if i_cum:
        cbar.set_label('Cumulative footprint contribution in percent')
    else:
        cbar.set_label('Percentage of footprint contribution')
    #cbar.set_label('Footprint in fraction')
    #cbar.set_label('Flux footprint contribution in fraction')
    ax.remove()
    plt.savefig(plot_path+'cbar.png',bbox_inches='tight') #, transparent=True)
    fn = plt.gcf().number
    #plt.clf()
    plt.close(fn)
    plt.ion()
    # get the lat/lon bounds of the area
    lon1 = lon[0]
    lon2 = lon[-1]
    lat1 = lat[0]
    lat2 = lat[-1]
    # Hopefully the file was opened properly and the header written
    fi.write('<GroundOverlay>\n')
    fi.write('  <name>'+station+zt2.strftime("%Y%m%d%H%M")+'</name>\n')
    fi.write('  <bgColor>8fffffff</bgColor>\n')
    fi.write('  <Icon>\n')
    fi.write('    <href>'+plot_in+'</href>\n')
    fi.write('  </Icon>\n')
    fi.write('  <TimeSpan>\n')
    fi.write('    <begin>'+zt1.strftime("%Y-%m-%dT%H:%M")+'</begin>\n')
    fi.write('    <end>'+zt2.strftime("%Y-%m-%dT%H:%M")+'</end>\n')
    fi.write('  </TimeSpan>\n')
    fi.write('  <altitude>0.0</altitude>\n')
    fi.write('  <altitudeMode>clampToGround</altitudeMode>\n')
    fi.write('  <LatLonBox>\n')
    fi.write('    <north>'+str(lat2)+'</north>\n')
    fi.write('    <south>'+str(lat1)+'</south>\n')
    fi.write('    <east>'+str(lon2)+'</east>\n')
    fi.write('    <west>'+str(lon1)+'</west>\n')
    fi.write('    <rotation>0.0</rotation>\n')
    fi.write('  </LatLonBox>\n')
    fi.write('</GroundOverlay>\n')

def kml_finalise(d,fi,mode,kmlname):
    # write the footer of the kml file and close the file
    fi.write("</Folder>\n")
    fi.write('</kml>\n')
    fi.close()
    # copy tower icon into the plot path directory to be added to the kmz file
    # create a kmz file out of the kml file
    cwd = os.getcwd()
    os.chdir(d["plot_path"])
    kmzname = kmlname.replace(".kml", ".kmz")
    msg = " Creating KMZ file " + kmzname
    logger.info(msg)
    plotlist = [p for p in os.listdir('.') if p.endswith(".png")]
    compression = zipfile.ZIP_DEFLATED
    zf = zipfile.ZipFile(kmzname, mode='w')
    zf.write(kmlname, compress_type=compression)
    os.remove(kmlname)
    for f in plotlist:
        zf.write(f, compress_type=compression)
        os.remove(f)
    zf.close()
    os.chdir(cwd)

def plotphifield(x, y, zt1, zt2, data, station, mode, clevs, imagename,i_cum):
    # plot footprint in 2-dim field; use x,y - coordinates
    text = 'Footprint ' + station + ' ' + zt1.strftime("%Y%m%d%H%M") + '  to  ' + zt2.strftime("%Y%m%d%H%M")
    plotname='plots/Footprint_'+ mode + zt1.strftime("%Y%m%d%H%M") + '.jpg'
    x_ll = x[0,0]   #xllcorner #-250
    x_ur = x[-1,-1] #xurcorner # 250
    y_ll = y[0,0]   #yllcorner #-250
    y_ur = y[-1,-1] #yurcorner # 250
    # create figure and axes instances
    plt.ion()
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    cs = plt.contourf(x,y,data,clevs,cmap=plt.get_cmap('hsv'))
    cbar = plt.colorbar(cs,location='right',pad=0.04,fraction=0.046)
    if i_cum:
        cbar.set_label('Cumulative footprint contribution in percent')
    else:
        cbar.set_label('Percentage of footprint contribution')
    # contour levels
    plt.title(text)
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    if imagename != 'None':
        #img = imread(imagename)
        img = plt.imread(imagename)
        plt.imshow(img, zorder=0, extent=[x_ll, x_ur, y_ll, y_ur])
    plt.savefig(plotname)
    plt.draw()
    plt.pause(1e-9)
    plt.ioff()
    
def calc_cumulative(f, f_min,f_step):
    # calculate the cumulative footprint values
    fcum05 = numpy.ma.masked_where(f <= f_min, f)
    fcum05 = numpy.ma.filled(fcum05,float(0))
    fcum  = numpy.sum(fcum05)
    num = int(round((1.0 - (f_min)) / f_step))
    ser1 = numpy.linspace(1.0-f_step,f_min,num)
    ser2 = numpy.linspace(1.0,f_min+f_step,num)
    ser3 = 0.5*(ser1+ser2)
    #print ser1,ser2
    cclevs = []
    stest = 0.0
    fmax=numpy.amax(f)
    for i in range(0,len(ser1)):
        test = numpy.ma.masked_where((f <= ser1[i]) | (f > ser2[i]), f)
        if test.count() > 0:
            test = numpy.sum(test)/fcum
        else:
            test = 0.0
        stest = stest + test
        cclevs.append(stest)
        #print ser3[i],stest #fmax,numpy.sum(fcum05),ser1[i],ser2[i],test,stest
    fcum_eq = numpy.polyfit(ser3,cclevs,3)
    #print fcum_eq
    fcum=fcum_eq[0]*f*f*f+fcum_eq[1]*f*f+fcum_eq[2]*f+fcum_eq[3]
    f=fcum
    
    return f
