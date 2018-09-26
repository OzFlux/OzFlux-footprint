import constants as c
import numpy
from footprint_utils import SeriestoMA, MAtoSeries

def absolutehumidityfromRH(Ta,RH):
    # convert to masked arrays
    RH, WasND = SeriestoMA(RH)
    Ta, dummy = SeriestoMA(Ta)
    # do the job
    VPsat = es(Ta)
    vp = RH * VPsat / float(100)
    ah = float(1000000) * vp / ((Ta + 273.15) * c.Rv)
    # convert back to ndarray if input is not a masked array
    if WasND: ah, WasMA = MAtoSeries(ah)
    return ah

def co2_ppmfrommgpm3(c_mgpm3,T,p):
    """
     Convert CO2 concentration units of mg/m3 to umol/mol (ppm)
        Usage:
         CO2_ppm = co2_ppmfrommgpm3(CO2_mgpm3, T, p)
         where
         CO2_mgpm3 (input) - CO2 concentration, mg/m3
         T (input) - air temperature, C
         p (input) - air pressure, kPa
        Returns the CO2 concentration in ppm.
    """
    # convert to masked array if required
    c_mgpm3, WasND = SeriestoMA(c_mgpm3)
    T, dummy = SeriestoMA(T)
    T = T + 273.15             # temperature in K
    p, dummy = SeriestoMA(p)
    p = p * float(1000)        # pressure in Pa
    # do the job
    c_ppm = (c_mgpm3/c.Mc)*c.R*T/p
    # convert back to ndarray if input is not a masked array
    if WasND: c_ppm, WasMA = MAtoSeries(c_ppm)
    return c_ppm

def co2_mgpm3fromppm(c_ppm,T,p):
    """
     Convert CO2 concentration units of umol/mol (ppm) to mg/m3
        Usage:
         CO2_mgpm3 = co2_mgpm3fromppm(CO2_ppm, T, p)
         where
         CO2_ppm (input) - CO2 concentration, umol/mol
         T (input) - air temperature, C
         p (input) - air pressure, kPa
        Returns the CO2 concentration in mg/m3.
    """
    # convert to masked array if required
    c_ppm, WasND = SeriestoMA(c_ppm)
    T, dummy = SeriestoMA(T)
    T = T + 273.15             # temperature in K
    p, dummy = SeriestoMA(p)
    p = p * float(1000)        # pressure in Pa
    # do the job
    c_mgpm3 = c_ppm*c.Mc*p/(c.R*T)
    # convert back to ndarray if input is not a masked array
    if WasND: c_mgpm3, WasMA = MAtoSeries(c_mgpm3)
    return c_mgpm3

def densitydryair(Ta,ps,vp):
    # Calculate density of dry air from temperature, pressure and vapour pressure
    #  Ta - air temperature, C
    #  ps - pressure, kPa
    #  vp - vapour pressure, kPa
    # Returns
    #  rhod - dry air density, kg/m3
    rhod = 1000*(ps-vp)/((Ta+273.15)*c.Rd)
    return rhod

def densitymoistair(Ta,ps,vp):
    # Calculate density of moist air from temperature, pressure and vapour pressure
    #  Ta - air temperature, C
    #  ps - pressure, kPa
    #  vp - vapour pressure, kPa
    # Returns
    #  rhom - moist air density, kg/m3
    rhow = vp*1000/((Ta+273.15)*c.Rv)
    rhod = (ps-vp)*1000/((Ta+273.15)*c.Rd)
    rhom = rhow + rhod
    return rhom

def densitywatervapour(Ta,vp):
    # Calculate partial density of water vapour from temperature and vapour pressure
    #  Ta - air temperature, C
    #  vp - vapour pressure, kPa
    # Returns
    #  rhow - partial density of water vapour, kg/m3
    rhow = vp*1000/((Ta+273.15)*c.Rv)
    return rhow

def es(T):
    # Saturation vapour pressure.
    #  T is the air temperature, C
    #  es is the saturation vapour pressure in kPa
    es = 0.6106 * numpy.exp(17.27 * T / (T + 237.3))
    return es

def Fc_umolpm2psfrommgpm2ps(Fc_mgpm2ps):
    """
    Convert Fc in units of mg/m2/s to units of umol/m2/s
    Usage:
     Fc_umolpm2ps = Fc_umolpm2psfrommgpm2ps(Fc_mgpm2ps)
     where:
      Fc_mgpm2ps (input) - CO2 flux in units of mg/m2/s
    Returns the CO2 flux in units of umol/m2/s
    """
    # convert to masked array
    Fc_mgpm2ps, WasND = SeriestoMA(Fc_mgpm2ps)
    # do the job
    Fc_umolpm2ps = Fc_mgpm2ps / c.Mc
    # convert back to ndarray if input is not a masked array
    if WasND: Fc_umolpm2ps, WasMA = MAtoSeries(Fc_umolpm2ps)
    return Fc_umolpm2ps

def Fc_mgpm2psfromumolpm2ps(Fc_umolpm2ps):
    """
    Convert Fc in units of umol/m2/s to units of mg/m2/s
    Usage:
     Fc_mgpm2ps = Fc_mgpm2psfromumolpm2ps(Fc_umolpm2ps)
     where:
      Fc_umolpm2ps (input) - CO2 flux in units of umol/m2/s
    Returns the CO2 flux in units of mg/m2/s
    """
    # convert to masked array
    Fc_umolpm2ps, WasND = SeriestoMA(Fc_umolpm2ps)
    # do the job
    Fc_mgpm2ps = Fc_umolpm2ps * c.Mc
    # convert back to ndarray if input is not a masked array
    if WasND: Fc_mgpm2ps, WasMA = MAtoSeries(Fc_mgpm2ps)
    return Fc_mgpm2ps

def h2o_mmolpmolfromgpm3(h_gpm3,T,p):
    """
     Convert H2O concentration units of g/m3 to mmol/mol.
        Usage:
         H2O_mmolpmol = h2o_mmolpmolfromgpm3(H2O_gpm3, T, p)
         where
         H2O_gpm3 (input) - H2O concentration, g/m3
         T (input) - air temperature, C
         p (input) - air pressure, kPa
        Returns the H2O concentration in mmol/mol.
    """
    # convert to masked arrays
    h_gpm3, WasND = SeriestoMA(h_gpm3)
    T, dummy = SeriestoMA(T)
    p, dummy = SeriestoMA(p)
    # do the job
    h_mmpm = (h_gpm3/c.Mv)*c.R*(T+273.15)/(p*1000)
    # convert to ndarray if input is not a masked array
    if WasND: h_mmpm, WasMA = MAtoSeries(h_mmpm)
    return h_mmpm

def h2o_gpm3frommmolpmol(h_mmpm,T,p):
    """
     Convert H2O concentration units of mmol/mol to g/m3.
        Usage:
         H2O_gpm3 = h2o_gpm3frommmolpmol(H2O_mmolpmol, T, p)
         where
         H2O_mmolpmol (input) - H2O concentration, mmol/mol
         T (input) - air temperature, C
         p (input) - air pressure, kPa
        Returns the H2O concentration in g/m3.
    """
    # convert to masked arrays
    h_mmpm, WasND = SeriestoMA(h_mmpm)
    T, dummy = SeriestoMA(T)
    p, dummy = SeriestoMA(p)
    # do the job
    h_gpm3 = (c.Mv*h_mmpm*p*1000)/(c.R*(T+273.15))
    # convert to ndarray if input is not a masked array
    if WasND: h_gpm3, WasMA = MAtoSeries(h_gpm3)
    return h_gpm3

def Lv(Ta):
    # Calculate Lv as a function of temperature, from Stull 1988
    #  Ta - air temperature, C
    # Returns
    #  Lv - latent heat of vapourisation, J/kg
    Lv = 2500800 - (2366.8 * Ta)
    return Lv

def mixingratio(ps,vp):
    # Calculate mixing ratio from vapour pressure and pressure
    #  ps - presure, kPa
    #  vp - vapour pressure, kPa
    # Returns
    #  mr - mixing ratio, kg/kg
    mr = 0.622*vp/(ps- vp)
    return mr

def molen(T,Ah,p,ustar,heatflux,fluxtype='sensible'):
    # Calculate the Monin-Obukhov length
    ustar = numpy.ma.sqrt(numpy.square(ustar))    # force the sign of ustar to be positive
    vp = vapourpressure(Ah,T)       # calculate the vapour pressure
    mr = mixingratio(p,vp)          # calculate the mixing ratio
    Tv = theta(T,p)                 # calculate potential temperature
    Tvp = virtualtheta(Tv,mr)
    if fluxtype=='sensible':
        L = -Tvp*densitydryair(T, p, vp)*c.Cp*(ustar**3)/(c.g*c.k*heatflux)
    elif fluxtype=='kinematic':
        L = -Tvp*(ustar**3)/(c.g*c.k*heatflux)
    else:
        raise Exception(" meteorologicalfunctions.molen: unkown value for fluxtype (="+str(fluxtype)+") encountered")
    return L

def qfromrh(RH, T, p):
    # Specific humidity (kg/kg) from relative humidity, temperature and pressure
    #  RH is the relative humidity, %
    #  T is the air temperature, C
    #  p is the atmospheric pressure, kPa
    qRH = (c.Mv / c.Md) * (0.01 * RH * es(T) / p)
    return qRH

def qsat(esat,ps):
    qsat = 0.622 * (esat / ps)
    return qsat

def RHfromabsolutehumidity(Ah,Ta):
    # Relative humidity from absolute humidity
    #  Ta is the air temperature, C
    #  Ah is the absolute humidity, g/m3
    #  RH is the relative humidity, %
    # convert to masked arrays
    Ah, WasND = SeriestoMA(Ah)
    Ta, dummy = SeriestoMA(Ta)
    # do the job
    VPsat = es(Ta)
    #vp = Ah * ((Ta+273.15)*c.Rv)/float(1000000)
    vp = vapourpressure(Ah,Ta)
    RH = float(100)*vp/VPsat
    # convert back to ndarray if input is not a masked array
    if WasND: RH, WasMA = MAtoSeries(RH)
    return RH

def RHfromdewpoint(Td,Ta):
    # Relative humidity from dew point temperature
    #  Ta is the air temperature, C
    #  Td is the dew point temperature, C
    #  RH is the relative humidity, %
    # convert to masked arrays
    Td, WasND = SeriestoMA(Td)
    Ta, dummy = SeriestoMA(Ta)
    # do the job
    RH = 100*10**(7.591386*(Td/(Td+240.7263)-Ta/(Ta+240.7263)))
    # convert back to ndarray if input is not a masked array
    if WasND: RH, WasMA = MAtoSeries(RH)
    return RH

def RHfromspecifichumidity(q,Ta,ps):
    # Relative humidity from specific humidity
    #  q is the specific humidity, kg/kg
    #  Ta is the air temperature, C
    #  ps is the pressure, kPa
    #  RH is the relative humidity, %
    # convert to masked arrays
    q, WasND = SeriestoMA(q)
    Ta, dummy = SeriestoMA(Ta)
    # do the job
    VPsat = es(Ta)
    RH = float(100)*q*(c.Md/c.Mv)*ps/VPsat
    # convert back to ndarray if input is not a masked array
    if WasND: RH, WasMA = MAtoSeries(RH)
    return RH

def densitytimesspecificheat(rhow,Cpw,rhoa,Cpa):
    '''
    Product of air density and specific heat capacity for moist air.
    '''
    return rhow*Cpw+rhoa*Cpa

def specificheatcapacitydryair(Tv):
    '''
    Specific heat capacity of air at constant pressure.
    USEAGE:
     cpd = mf.cpd(Tv)
    INPUT:
     Tv - virtual temperature (from sonic anemometer), C
    OUTPUT:
     cpd - specific heat capacity of dry air at constant pressure, J/kg/K
    SOURCE:
     EddyPro source code
    '''
    cpd = float(1005)+((Tv+23.12)**2)/float(3364)
    return cpd

def specificheatcapacitywatervapour(Ta, Ah):
    '''
    Specific heat capacity of water vapour at constant pressure.
    USEAGE:
     cpv = mf.cpv(Ta,Ah)
    INPUT:
     Ta - air temperature, C
     Ah - absolute humidity, %
    OUTPUT:
     cpv - specific heat capacity of water vapour at constant pressure, J/kg/K
    SOURCE:
     EddyPro source code
    '''
    RH = RHfromabsolutehumidity(Ah,Ta)
    cpv = float(1859)+0.13*RH+(0.193+0.0056*RH)*Ta+(0.001+0.00005*RH)*Ta**2
    return cpv

def specificheatmoistair(q):
    # Calculate Cp of moist air, from Stull 1988
    #  Cp - specific heat of dry air at constant pressure, J/kg-K
    #  q - specific humidity
    # Returns
    #  Cpm - specific heat of moist air at constant pressure, J/kg-K
    Cpm = c.Cpd * (1 + 0.84 * q)
    return Cpm

def specifichumidity(mr):
    # Calculate specific humidity from mixing ratio
    #  mr - mixing ration, kg/kg
    # Returns
    #  q = specific humidity, kg/kg
    q = mr/(1+mr)
    return q

def specifichumidityfromRH(RH, T, p):
    # Specific humidity (kg/kg) from relative humidity, temperature and pressure
    #  RH is the relative humidity, %
    #  T is the air temperature, C
    #  p is the atmospheric pressure, kPa
    # Returns
    #  q = specific humidity, kg/kg
    q = (c.Mv / c.Md) * (0.01 * RH * es(T) / p)
    return q

def tafromtv(Tv,q):
    # Calculate air temperature from virtual temperature using formula
    # from Campbell Scientific CSAT manual.
    # NOTE: this differs from the usual definition by using 0.51 not 0.61
    #  Tv - virtual temperature, C
    #  q - specific humidity, kg/kg
    # Returns
    #  Ta - air temperature, C
    Ta = ((Tv+273.15)/(1+0.51*q))-273.15
    return Ta

def theta(T,p):
    # Calculate potential temperature from air temperature and pressure
    #  T - air temperature, C
    #  p - pressure, kPa
    # Returns
    #  theta - potential temperature, K
    return (T+273.15)*(100/p)**0.286

def vapourpressure(Ah,Ta):
    # Calculate vapour pressure from absolute humidity and temperature
    #  Ah - absolute humidity, g/m3
    #  Ta - air temperature, C
    # Returns
    #  vp - vapour pressure, kPa
    vp = 0.000001*Ah*(Ta+273.15)*c.R/c.Mv
    return vp

def virtualtheta(theta,mr):
    # Calculate virtual potential temperature
    #  theta - potential temperature, K
    #  mr - mixing ratio, kg/kg
    # Returns
    #  Tvp - virtual potential temperature, K
    Tvp = theta * (1 + (0.61 * mr))
    return Tvp
