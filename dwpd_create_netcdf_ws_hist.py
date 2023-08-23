'''
dwpd_create_netcdf_ws_hist.py
Load dill pickle containing return levels and parameters for historical period.
Output NetCDF with suitable metadata.

Written by Stephen Outten June 2022.
'''

import numpy as np
import dill
from netCDF4 import Dataset


# FolderName = # Path removed in online shared codes for security
FolderName = # Path removed in online shared codes for securit
FolderNameNetCDF = # Path removed in online shared codes for security
# FolderNameOut = # Path removed in online shared codes for security
FolderNameOut = # Path removed in online shared codes for security
Rx = np.array([2,5,10,20,30,50,70,100,200])
params = np.array(['location', 'scale', 'shape'])
# DS_Name = 'CLMcom_CNRM', 'CLMcom_ICHEC', 'CLMcom_MOHC', 'CLMcom_MPI', 'CNRM_CNRM', 'CNRM_MOHC', 'CNRM_MPI', 'CNRM_NorESM', 'DMI_CNRM', 'DMI_ICHEC', 'DMI_IPSL', 'DMI_MOHC', 'DMI_MPI', 'DMI_NorESM', 'ETH_CNRM', 'ETH_ICHEC', 'ETH_MOHC', 'ETH_MPI', 'ETH_NorESM', 'GERICS_CNRM', 'GERICS_IPSL', 'GERICS_MPI', 'GERICS_NorESM', 'IPSL_CNRM', 'IPSL_ICHEC', 'IPSL_IPSL', 'IPSL_MOHC', 'IPSL_MPI', 'IPSL_NorESM', 'KNMI_CNRM', 'KNMI_ICHEC', 'KNMI_IPSL', 'KNMI_MOHC', 'KNMI_MPI', 'KNMI_NorESM', 'MOHC_CNRM', 'MOHC_ICHEC', 'MOHC_MOHC', 'MOHC_MPI', 'MOHC_NorESM', 'MPI_MPI', 'SMHI_CNRM', 'SMHI_ICHEC', 'SMHI_IPSL', 'SMHI_MOHC', 'SMHI_MPI', 'SMHI_NorESM'
# End_Name = '20051231', '20051231', '20051230', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051230', '20051231', '20051231', '20051231', '20051231', '20051230', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051230', '20051231', '20051231', '20051231', '20051231', '20051130', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051230', '20051231', '20051231'

DS_Name = 'CLMcom_CNRM', 'CNRM_CNRM', 'DMI_CNRM'
End_Name = '20051231', '20051231', '20051231'


for dsi in range(len(DS_Name)):

    # Load return elvels and parameters from dill pickle
    FileName = DS_Name[dsi] + '_ws_hist_gpd_maps.dpkl'
    print('Loading...' + FileName)
    fn = FolderName + FileName
    lat = dill.load(open(fn,'rb'))['lat']
    lon = dill.load(open(fn,'rb'))['lon']
    GEV_R2 = dill.load(open(fn,'rb'))['GEV_R2']
    GEV_R5 = dill.load(open(fn,'rb'))['GEV_R5']
    GEV_R10 = dill.load(open(fn,'rb'))['GEV_R10']
    GEV_R20 = dill.load(open(fn,'rb'))['GEV_R20']
    GEV_R30 = dill.load(open(fn,'rb'))['GEV_R30']
    GEV_R50 = dill.load(open(fn,'rb'))['GEV_R50']
    GEV_R70 = dill.load(open(fn,'rb'))['GEV_R70']
    GEV_R100 = dill.load(open(fn,'rb'))['GEV_R100']
    GEV_R200 = dill.load(open(fn,'rb'))['GEV_R200']
    GEV_location = dill.load(open(fn,'rb'))['GEV_location']
    GEV_scale = dill.load(open(fn,'rb'))['GEV_scale']
    GEV_shape = dill.load(open(fn,'rb'))['GEV_shape']
    GEV_Return_Levels = np.stack((GEV_R2, GEV_R5, GEV_R10, GEV_R20, GEV_R30, GEV_R50, GEV_R70, GEV_R100, GEV_R200), axis=0)
    GEV_Parameters = np.stack((GEV_location, GEV_scale, GEV_shape), axis=0)
 
    GPD_R2 = dill.load(open(fn,'rb'))['GPD_R2']
    GPD_R5 = dill.load(open(fn,'rb'))['GPD_R5']
    GPD_R10 = dill.load(open(fn,'rb'))['GPD_R10']
    GPD_R20 = dill.load(open(fn,'rb'))['GPD_R20']
    GPD_R30 = dill.load(open(fn,'rb'))['GPD_R30']
    GPD_R50 = dill.load(open(fn,'rb'))['GPD_R50']
    GPD_R70 = dill.load(open(fn,'rb'))['GPD_R70']
    GPD_R100 = dill.load(open(fn,'rb'))['GPD_R100']
    GPD_R200 = dill.load(open(fn,'rb'))['GPD_R200']
    GPD_location = dill.load(open(fn,'rb'))['GPD_location']
    GPD_scale = dill.load(open(fn,'rb'))['GPD_scale']
    GPD_shape = dill.load(open(fn,'rb'))['GPD_shape']
    GPD_Return_Levels = np.stack((GPD_R2, GPD_R5, GPD_R10, GPD_R20, GPD_R30, GPD_R50, GPD_R70, GPD_R100, GPD_R200), axis=0)
    GPD_Parameters = np.stack((GPD_location, GPD_scale, GPD_shape), axis=0) 

    # Load relevant metadata from NetCDF
    fn = FolderNameNetCDF + 'sfcWindmax_' + DS_Name[dsi] + '_19760101_' + End_Name[dsi] + '.nc'
    print(fn)
    dataset = Dataset(fn)
    driving_model = dataset.driving_model_id
    model = dataset.model_id
    ensemble = dataset.driving_model_ensemble_member
    dataset.close()
    frequency = 'day'
    experiment = 'historical'
    variable = 'Daily Maximum Near-Surface Wind Speed'     


    ##########################################################################
    ###         Write data out to NetCDF                                   ###
    ##########################################################################

    # Determine size of dimensions
    xx, yy  = lat.shape
    nrx, = Rx.shape
    nparam, = params.shape
   
    # open a new netCDF file for writing.
    FileNameOut =  'sfcWindmax_params_rx_historical_' + DS_Name[dsi] + '.nc'
    fnout = FolderNameOut + FileNameOut
    ncfile = Dataset(fnout, 'w')

    # create the lat and lon dimensions.
    x = ncfile.createDimension('x', xx)
    y = ncfile.createDimension('y', yy)
    return_levels = ncfile.createDimension('return_levels', nrx)    
    parameters = ncfile.createDimension('parameters', nparam)

    # Define the coordinate variables.
    x = ncfile.createVariable('x', np.float32, ('x'))
    y = ncfile.createVariable('y', np.float32, ('y'))
    return_levels = ncfile.createVariable('return_levels', np.float32, ('return_levels'))
    parameters = ncfile.createVariable('parameters', str, ('parameters'))


    # Assign units attributes to coordinate variable data.
    x.units = ''
    y.units = ''
    return_levels.units = 'years'
    parameters.units = ''
    x.long_name = 'Number of grid points in x direciton'
    y.long_name = 'Number of grid points in y direction'
    return_levels.long_name = 'Return levels'
    parameters.long_name = 'Names of distribution parameters'

    # write data to coordinate vars.
    x[:] = np.arange(xx) + 1
    y[:] = np.arange(yy) + 1
    return_levels[:] = Rx
    parameters[:] = params
 
    # create main variables
    lat_out = ncfile.createVariable('lat', np.float32, ('x','y'))
    lon_out = ncfile.createVariable('lon', np.float32, ('x','y'))
    gev_parameters_out = ncfile.createVariable('gev_parameters', np.float32, ('parameters','x','y'))
    gev_rx_out = ncfile.createVariable('gev_return_levels', np.float32, ('return_levels','x','y'))
    gpd_parameters_out = ncfile.createVariable('gpd_parameters', np.float32, ('parameters','x','y'))
    gpd_rx_out = ncfile.createVariable('gpd_return_levels', np.float32, ('return_levels','x','y'))

    # set the units attribute.
    lat_out.units = 'degrees_north'
    lat_out.standard_name = 'latitude'
    lat_out.long_name = 'latitude coordinate'
    lon_out.units = 'degrees_east'
    lon_out.standard_name = 'longitude'
    lon_out.long_name = 'longitude coordinate'
    gev_parameters_out.units = ''
    gev_parameters_out.long_name = 'Parameters of fitted GEV distributions'
    gev_rx_out.units = 'm s-1'
    gev_rx_out.long_name = 'Return level wind speeds based on fitted GEV distribution'
    gpd_parameters_out.units = ''
    gpd_parameters_out.long_name = 'Parameters of fitted GPD distributions'
    gpd_rx_out.units = 'm s-1'
    gpd_rx_out.long_name = 'Return level wind speeds based on fitted GPD distribution'


    # write data to variables along record (unlimited) dimension.
    lat_out[:] = lat
    lon_out[:] = lon
    gev_parameters_out[:] = GEV_Parameters
    gev_rx_out[:] = GEV_Return_Levels
    gpd_parameters_out[:] = GPD_Parameters
    gpd_rx_out[:] = GPD_Return_Levels

    # Write global parameters
    ncfile.comment = "Return levels and parameters for GEV and GPD distributions fitted to daily maximum near-surface wind speed in EURO-CORDEX 0.11 simulation"
    ncfile.driving_model_id = driving_model
    ncfile.model_id = model
    ncfile.driving_model_ensemble_member = ensemble
    ncfile.frequency = frequency
    ncfile.experiment = experiment
    ncfile.variable_fitted = variable

    # close the file.
    ncfile.close()

