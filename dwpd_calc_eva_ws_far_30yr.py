'''
dwpd_calc_eva_ws_far_30yr.py
Loads annual maxima and exceedances.
Calculates and save the return levels for 2, 5, 10, 20, 30, 50, 70, 100, and 200 years.
Calculates this for GEV and GPD.
Repeat it for every point in the domain and save it.
Repeat it for every model.
Repeat it for the far period of the future.

Output will be stored as dill pickle which will be read by a separate file that will create NetCDF with metadata.
The data from here was not output as NetCDF initially so that this code could be running asap, 
while I work on the code to create NetCDF with all the releveant metadata.

Written by Stephen Outten June 2022
'''

import numpy as np
import dill
from scipy.stats import genextreme as gev
from scipy.stats import genpareto as gpd
from scipy.stats import chi2
import time
import warnings

begin_time = time.time()
# Turn off the warnings   (e.g. "all-nan slice encountered")
warnings.filterwarnings('ignore')

FolderName = # Path removed in online shared codes for security
DS_Name = 'CLMcom_CNRM', 'CLMcom_ICHEC', 'CLMcom_MOHC', 'CLMcom_MPI', 'CNRM_CNRM', 'CNRM_MOHC', 'CNRM_MPI', 'CNRM_NorESM', 'DMI_CNRM', 'DMI_ICHEC', 'DMI_IPSL', 'DMI_MOHC', 'DMI_MPI', 'DMI_NorESM', 'ETH_CNRM', 'ETH_ICHEC', 'ETH_MOHC', 'ETH_MPI', 'ETH_NorESM', 'GERICS_CNRM', 'GERICS_IPSL', 'GERICS_MPI', 'GERICS_NorESM', 'ICTP_CNRM', 'ICTP_ICHEC', 'ICTP_MOHC', 'ICTP_MPI', 'ICTP_NorESM', 'IPSL_CNRM', 'IPSL_ICHEC', 'IPSL_IPSL', 'IPSL_MOHC', 'IPSL_MPI', 'IPSL_NorESM', 'KNMI_CNRM', 'KNMI_ICHEC', 'KNMI_IPSL', 'KNMI_MOHC', 'KNMI_MPI', 'KNMI_NorESM', 'MOHC_CNRM', 'MOHC_ICHEC', 'MOHC_MOHC', 'MOHC_MPI', 'MOHC_NorESM', 'MPI_MPI', 'SMHI_CNRM', 'SMHI_ICHEC', 'SMHI_IPSL', 'SMHI_MOHC', 'SMHI_MPI', 'SMHI_NorESM'
# DS_Name = 'ICTP_CNRM', 'ICTP_ICHEC', 'ICTP_MOHC', 'ICTP_MPI', 'ICTP_NorESM'

Period_Name = 'near', 'mid', 'far'
period = 2   # Far 

for dsi in range(len(DS_Name)):
    FileName = 'ws_fut_' + Period_Name[0] + '_' + DS_Name[dsi] + '.dpkl'
    fn = FolderName + FileName
    lat = dill.load(open(fn,'rb'))['lat']
    a,b = lat.shape
    
    # Create empty arrays for return values
    GEV_R2 = np.ones([a,b])*np.nan
    GEV_R5 = np.ones([a,b])*np.nan
    GEV_R10 = np.ones([a,b])*np.nan
    GEV_R20 = np.ones([a,b])*np.nan
    GEV_R30 = np.ones([a,b])*np.nan
    GEV_R50 = np.ones([a,b])*np.nan
    GEV_R70 = np.ones([a,b])*np.nan
    GEV_R100 = np.ones([a,b])*np.nan
    GEV_R200 = np.ones([a,b])*np.nan
    GEV_location = np.ones([a,b])*np.nan
    GEV_scale = np.ones([a,b])*np.nan
    GEV_shape = np.ones([a,b])*np.nan

    GPD_R2 = np.ones([a,b])*np.nan
    GPD_R5 = np.ones([a,b])*np.nan
    GPD_R10 = np.ones([a,b])*np.nan
    GPD_R20 = np.ones([a,b])*np.nan
    GPD_R30 = np.ones([a,b])*np.nan
    GPD_R50 = np.ones([a,b])*np.nan
    GPD_R70 = np.ones([a,b])*np.nan
    GPD_R100 = np.ones([a,b])*np.nan
    GPD_R200 = np.ones([a,b])*np.nan
    GPD_location = np.ones([a,b])*np.nan
    GPD_scale = np.ones([a,b])*np.nan
    GPD_shape = np.ones([a,b])*np.nan
   

    # Load data from one of the models 
    FileName = 'ws_fut_' + Period_Name[period] + '_' + DS_Name[dsi] + '.dpkl'
    print('Loading...' + FileName)
    fn = FolderName + FileName
    lat = dill.load(open(fn,'rb'))['lat']
    lon = dill.load(open(fn,'rb'))['lon']
    exs = dill.load(open(fn,'rb'))['exceeds']
    exs_lambda = dill.load(open(fn,'rb'))['exc_lambda']
    threshold = dill.load(open(fn,'rb'))['threshold']
    amx = dill.load(open(fn,'rb'))['annual_max']

    # Loop over all locations
    start_time = time.time()
    for xx in range(a):
        # print('Row number ' + np.str(xx))
        for yy in range(b):
            # Extract single location
            vecamx = amx[:,xx,yy]
            vecexs = exs[xx,yy]
            veclmd = exs_lambda[xx,yy]
            vecth = threshold[xx,yy]

            # Fit GEV and calculate return values
            shp,loc,scl = gev.fit(vecamx)
            GEV_location[xx,yy] = loc
            GEV_scale[xx,yy] = scl
            GEV_shape[xx,yy] = shp
            GEV_R2[xx,yy] = gev.isf(1/(2),shp,loc=loc,scale=scl)
            GEV_R5[xx,yy] = gev.isf(1/(5),shp,loc=loc,scale=scl)
            GEV_R10[xx,yy] = gev.isf(1/(10),shp,loc=loc,scale=scl)
            GEV_R20[xx,yy] = gev.isf(1/(20),shp,loc=loc,scale=scl)
            GEV_R30[xx,yy] = gev.isf(1/(30),shp,loc=loc,scale=scl)
            GEV_R50[xx,yy] = gev.isf(1/(50),shp,loc=loc,scale=scl)
            GEV_R70[xx,yy] = gev.isf(1/(70),shp,loc=loc,scale=scl)
            GEV_R100[xx,yy] = gev.isf(1/(100),shp,loc=loc,scale=scl)
            GEV_R200[xx,yy] = gev.isf(1/(200),shp,loc=loc,scale=scl)

            # Fit GPD and calculate return values
            shp,loc,scl = gpd.fit(vecexs)
            GPD_location[xx,yy] = loc
            GPD_scale[xx,yy] = scl
            GPD_shape[xx,yy] = shp
            GPD_R2[xx,yy] = gpd.isf(1/(veclmd*2),shp,loc=loc,scale=scl) + vecth
            GPD_R5[xx,yy] = gpd.isf(1/(veclmd*5),shp,loc=loc,scale=scl) + vecth
            GPD_R10[xx,yy] = gpd.isf(1/(veclmd*10),shp,loc=loc,scale=scl) + vecth
            GPD_R20[xx,yy] = gpd.isf(1/(veclmd*20),shp,loc=loc,scale=scl) + vecth
            GPD_R30[xx,yy] = gpd.isf(1/(veclmd*30),shp,loc=loc,scale=scl) + vecth
            GPD_R50[xx,yy] = gpd.isf(1/(veclmd*50),shp,loc=loc,scale=scl) + vecth
            GPD_R70[xx,yy] = gpd.isf(1/(veclmd*70),shp,loc=loc,scale=scl) + vecth
            GPD_R100[xx,yy] = gpd.isf(1/(veclmd*100),shp,loc=loc,scale=scl) + vecth
            GPD_R200[xx,yy] = gpd.isf(1/(veclmd*200),shp,loc=loc,scale=scl) + vecth

    print('This took {0:.6} hours to calculate'.format([time.time()-start_time][0]/3600))


    # Store lat, lon, and return values as dill pickle
    fnout = # Path removed in online shared codes for security + DS_Name[dsi] + '_ws_' + Period_Name[period] + '_gpd_maps.dpkl'
    with open(fnout, 'wb') as f:
        dill.dump({'lat':lat, 'lon':lon, 'GEV_R2':GEV_R2, 'GEV_R5':GEV_R5, 'GEV_R10':GEV_R10, 'GEV_R20':GEV_R20, 'GEV_R30':GEV_R30, 'GEV_R50':GEV_R50, 'GEV_R70':GEV_R70, 'GEV_R100':GEV_R100, 'GEV_R200':GEV_R200, 'GEV_location':GEV_location, 'GEV_scale':GEV_scale, 'GEV_shape':GEV_shape, 'GPD_R2':GPD_R2, 'GPD_R5':GPD_R5, 'GPD_R10':GPD_R10, 'GPD_R20':GPD_R20, 'GPD_R30':GPD_R30, 'GPD_R50':GPD_R50, 'GPD_R70':GPD_R70, 'GPD_R100':GPD_R100, 'GPD_R200':GPD_R200, 'GPD_location':GPD_location, 'GPD_scale':GPD_scale, 'GPD_shape':GPD_shape},f)

    del amx, exs, exs_lambda, threshold, vecexs, veclmd, vecth, GEV_R2, GEV_R5, GEV_R10, GEV_R20, GEV_R30, GEV_R50, GEV_R70, GEV_R100, GEV_R200, GPD_R2, GPD_R5, GPD_R10, GPD_R20, GPD_R30, GPD_R50, GPD_R70, GPD_R100, GPD_R200, lat, lon, GEV_location, GEV_scale, GEV_shape, GPD_location, GPD_scale, GPD_shape



print('Final run time was {0:.6} hours'.format([time.time()-begin_time][0]/3600))





