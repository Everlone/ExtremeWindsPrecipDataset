'''
dwpd_ws_hist_amax_exceed_timing_30yr.py
Loads max daily wind speed for CORDEX models and calculated GPD.
Performs the calculation using 30 years of data from 1976-2005.
Also stored the timing for the occurence of each peak.
This could be used for analysing the co-occurence of extreme wind and precip events.
Output is stored for each model separately in pickle.

NOTE:   Declustering originally relied on winds of specific temporal spacing not being exactly the same.
        e.g. declustering [0,0,0,4,3,4,3,0,0,0] gives no peaks.
        This is because the two 4s are exactly the same so there difference produces zero.
        They are within the spacing of each other and there is no 'peak' - i.e. they are equal.
        This is solves by adding random tiny number to one and removing it later.

Data is stored as dictionary in dpkl so variables can be laoded by name.
Load data into other .py scripts using the following:
*****************************************************
lat = dill.load(open('FILENAME','rb'))['lat']
lon = dill.load(open('FILENAME','rb'))['lon']
amx = dill.load(open('FILENAME','rb'))['annual_max']
exs = dill.load(open('FILENAME','rb'))['exceeds']
*****************************************************

NOTE: MOHC_MOHC needs to be dealt with separately as it only goes until 20051130, hence the calendar block
needs to be built differrently. This will most likely be done interactively by hand, possibly then coded in and commented out.

Written by Stephen Outten June 2022
'''

from netCDF4 import Dataset
import numpy as np
import dill

FolderName = # Path removed in online shared codes for security
# DS_Name = 'CLMcom_CNRM', 'CLMcom_ICHEC', 'CLMcom_MOHC', 'CLMcom_MPI', 'CNRM_CNRM', 'CNRM_MOHC', 'CNRM_MPI', 'CNRM_NorESM', 'DMI_CNRM', 'DMI_ICHEC', 'DMI_IPSL', 'DMI_MOHC', 'DMI_MPI', 'DMI_NorESM', 'ETH_CNRM', 'ETH_ICHEC', 'ETH_MOHC', 'ETH_MPI', 'ETH_NorESM', 'GERICS_CNRM', 'GERICS_IPSL', 'GERICS_MPI', 'GERICS_NorESM', 'IPSL_CNRM', 'IPSL_ICHEC', 'IPSL_IPSL', 'IPSL_MOHC', 'IPSL_MPI', 'IPSL_NorESM', 'KNMI_CNRM', 'KNMI_ICHEC', 'KNMI_IPSL', 'KNMI_MOHC', 'KNMI_MPI', 'KNMI_NorESM', 'MOHC_CNRM', 'MOHC_ICHEC', 'MOHC_MPI', 'MOHC_NorESM', 'MPI_MPI', 'SMHI_CNRM', 'SMHI_ICHEC', 'SMHI_IPSL', 'SMHI_MOHC', 'SMHI_MPI', 'SMHI_NorESM'
# End_Name = '20051231', '20051231', '20051230', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051230', '20051231', '20051231', '20051231', '20051231', '20051230', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051230', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051231', '20051230', '20051231', '20051231'

DS_Name = 'ICTP_CNRM', 'ICTP_ICHEC', 'ICTP_MOHC', 'ICTP_MPI', 'ICTP_NorESM'
End_Name = '20051231', '20051231', '20051230', '20051231', '20051231'



for dsi in range(len(DS_Name)):

    # Load data from selected model
    fn = FolderName + 'sfcWindmax_' + DS_Name[dsi] + '_19760101_' + End_Name[dsi] + '.nc'
    print(fn)
 
    dataset = Dataset(fn)
    lon = dataset.variables['lon'][:].data
    lat = dataset.variables['lat'][:].data
    # Determine annual blocks for calendar used by model and load data for 30 years only
    if dataset.variables['time'].calendar == '360_day':
        calblock = np.cumsum([360]*29)
        time = dataset.variables['time'][:].data
        wsm = dataset.variables['sfcWindmax'][:,:,:].data
    elif dataset.variables['time'].calendar == '365_day':
        calblock = np.cumsum([365]*29)
        time = dataset.variables['time'][:].data
        wsm = dataset.variables['sfcWindmax'][:,:,:].data        
    else:
        calblock = np.cumsum(([366,365,365,365]*8)[:-3])    # normal calendar
        time = dataset.variables['time'][:].data
        wsm = dataset.variables['sfcWindmax'][:,:,:].data



    # Split into years based on calendar blocks and take maxima of each year
    annual_max = []
    [annual_max.append(np.nanmax(x,0)) for x in np.split(wsm,calblock)]
    annual_max = np.rollaxis(np.dstack(annual_max),-1)


    # Create array for exceedances over threshold (smallest annual maxima)
    # Needs to be object type to allow assignment of different length arrays
    # Before declustering, some values may be exactly the same, even after
    # greater reducing the number of them after removing the threshold
    # Solve this by identfying which numbers are equal and not zero
    # Add tiny random numbers to these so declustering always sees peaks
    # and not plateus which it can't separate
    # Testing showed rand*0.0001 was the smallest value that guaranteed
    # the value changed given it is float32
    # The tiny value is at most 1/10,000th of a mm
    # Decluster by identifying difference from n to n+-2 and n+-1
    # If all positive (i.e. min(all) > 0), then this is local peak, if not set to 0
    # Convert to mask and multiply by timeseries
    [a,b,c] = wsm.shape
    exc_lambda = np.zeros([b,c])
    exceeds = np.array(np.zeros([b,c]), dtype=object)  
    exc_time = np.array(np.zeros([b,c]), dtype=object)
    threshold = np.zeros([b,c])
    winlen = 2  # +- size of windows around value, e.g. +-2 for wind, +-4 for pr5d
    for row in range(b):
        for col in range(c):
            threshold[row,col] = np.nanmin(annual_max[:,row,col])
            vec = wsm[:,row,col] - threshold[row,col]
            vec[vec<0] = 0

            # Check for duplicate peaks within window and add small random number
            eqpts = np.logical_or.reduce([np.equal(vec[:-winlen],vec[i+1:len(vec)-winlen+(i+1)]) for i in range(winlen)])
            eqpts = np.logical_and(np.not_equal(vec[:-winlen],0),eqpts)
            eqpts = eqpts * np.random.rand(len(eqpts)) * 0.0001
            eqpts = np.concatenate((eqpts,[0]*winlen))
            vec = vec + eqpts

            # Isolate peaks
            vec = np.concatenate(([0]*winlen, vec, [0]*winlen))
            win1 = [vec[winlen:-winlen]-vec[winlen-(i+1):-winlen-(i+1)] for i in range(winlen)]
            win2 = [vec[winlen:-winlen]-vec[winlen+(i+1):len(vec)-winlen+(i+1)] for i in range(winlen)]
            peaks = np.nanmin(np.vstack([win1,win2]),0)
            peaks[peaks<0] = 0
            peaks[peaks>0] = 1
            peak_time = time * peaks
            exc_time[row,col] = peak_time[np.where(peak_time>0)]

            # Count peaks
            cnt_peaks = []
            [cnt_peaks.append(np.nansum(x)) for x in np.split(peaks,calblock)]

            # Restore value of peaks and remove random numbeer, store results
            peaks = peaks * vec[winlen:-winlen]
            eqpts[peaks==0] = 0   # Remove random number from selected peak
            peaks = peaks-eqpts
            exceeds[row,col] = peaks[peaks>0]
            exc_lambda[row,col] = np.nanmean(cnt_peaks)    # identical to Poisson fit


    # Store lat, lon, annual maxima and exceedances as dill pickle 
    fnout = # Path removed in online shared codes for securit + 'ws_hist_' + DS_Name[dsi] + '.dpkl'
    with open(fnout, 'wb') as f:
        dill.dump({'lat':lat, 'lon':lon, 'annual_max':annual_max, 'exceeds':exceeds, 'exc_lambda':exc_lambda, 'threshold':threshold, 'exc_time':exc_time},f)


    # Clean up
    dataset.close()
    del dataset, wsm, lat, lon, annual_max, exceeds, calblock, cnt_peaks, exc_lambda, win1, win2, winlen, vec, peaks, exc_time, peak_time, time



