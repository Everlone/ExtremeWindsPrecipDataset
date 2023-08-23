'''
dwpd_ws_fut_amax_exceed_timing_30yr.py
Loads max daily wind speed for CORDEX models and calculated GPD.
Performs the calculation using 3 periods, each of 30 years of data.
The periods are 2011-204, 2041-2070, 2071-2100.
Also stored the timing for the occurence of each peak.
This could be used for analysing the co-occurence of extreme wind and precip events.
Output is stored for each model and for each period separately in pickle.

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

NOTE: Five models need to be handled separately as they each have different calendars for the final period.
The calendar block for the final period for each of these nees to be created separately.
This file will be run for the five models for the first two periods separately. 
For the final period, the calendar blocks will be created manually, probably interactively. 
These may be stored here in code and commented out, but possibly not. 
The five models and their end points are
    MOHC_MOHC   21001219
    MOHC_CNRM   21001230  only missing single last day
    IPSL_MOHC   20991201
    ETH_MOHC    20981230
    CNRM_MOHC   20991231

Written by Stephen Outten June 2022
'''

from netCDF4 import Dataset
import numpy as np
import dill

FolderName =  # Path removed in online shared codes for security
Period_Folder = '2011_2040/', '2041_2070/', '2071_2100/'
Period_Name = 'near', 'mid', 'far'
Start_Name = '20110101', '20410101', '20710101'

#DS_Name = 'CLMcom_CNRM', 'CLMcom_ICHEC', 'CLMcom_MOHC', 'CLMcom_MPI', 'CNRM_CNRM', 'CNRM_MPI', 'CNRM_NorESM', 'DMI_CNRM', 'DMI_ICHEC', 'DMI_IPSL', 'DMI_MOHC', 'DMI_MPI', 'DMI_NorESM', 'ETH_CNRM', 'ETH_ICHEC', 'ETH_MPI', 'ETH_NorESM', 'GERICS_CNRM', 'GERICS_IPSL', 'GERICS_MPI', 'GERICS_NorESM', 'IPSL_CNRM', 'IPSL_ICHEC', 'IPSL_IPSL', 'IPSL_MPI', 'IPSL_NorESM', 'KNMI_CNRM', 'KNMI_ICHEC', 'KNMI_IPSL', 'KNMI_MOHC', 'KNMI_MPI', 'KNMI_NorESM', 'MOHC_ICHEC', 'MOHC_MPI', 'MOHC_NorESM', 'MPI_MPI', 'SMHI_CNRM', 'SMHI_ICHEC', 'SMHI_IPSL', 'SMHI_MOHC', 'SMHI_MPI', 'SMHI_NorESM', 'MOHC_MOHC', 'MOHC_CNRM', 'IPSL_MOHC', 'ETH_MOHC', 'CNRM_MOHC'
#End_Name = []
#End_Name.append(('20401231', '20401231', '20401230', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401230', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401230', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401231', '20401230', '20401231', '20401231', '20401230', '20401231', '20401231', '20401230', '20401231'))
#End_Name.append(('20701231', '20701231', '20701230', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701230', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701230', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701231', '20701230', '20701231', '20701231', '20701230', '20701231', '20701231', '20701230', '20701231'))
#End_Name.append(('21001231', '21001231', '20991230', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '20991230', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '20991230', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '21001231', '20991230', '21001231', '21001231', '20991219', '21001230', '20991201', '20981230', '20991231'))


DS_Name = 'ICTP_CNRM', 'ICTP_ICHEC', 'ICTP_MOHC', 'ICTP_MPI', 'ICTP_NorESM'
End_Name = []
End_Name.append(('20401231', '20401231', '20401230', '20401231', '20401231'))
End_Name.append(('20701231', '20701231', '20701230', '20701231', '20701231'))
End_Name.append(('21001231', '21001231', '20991230', '21001231', '21001231'))



for period in range(len(Period_Folder)):
    for dsi in range(len(DS_Name)):

        # Load data from selected model
        fn = FolderName + Period_Folder[period] + 'sfcWindmax_' + DS_Name[dsi] + '_' + Start_Name[period] + '_' + End_Name[period][dsi] + '.nc'
        print(fn)
 
        dataset = Dataset(fn)
        lon = dataset.variables['lon'][:].data
        lat = dataset.variables['lat'][:].data
        time = dataset.variables['time'][:].data
        wsm = dataset.variables['sfcWindmax'][:].data

        # Determine annual blocks for calendar used by model and load data for 30 years only
        if dataset.variables['time'].calendar == '360_day':
            # The final period of MOHC models (calendar 360 days) is 29 years not 30 years long
            if period == 2:     
                calblock = np.cumsum([360]*28)
            else:
                calblock = np.cumsum([360]*29)
            if End_Name[period][dsi][0:4] == '2098':
                calblock = np.cumsum([360]*27)
        elif dataset.variables['time'].calendar == '365_day':
            calblock = np.cumsum([365]*29)
        else:
            # Set calendar depending upon when first leap year occurs in each period
            if period == 1:
                calblock = np.cumsum(([365,365,365,366]*8)[:-3])    # normal calendar, 4th year leap
            else:
                calblock = np.cumsum(([365,366,365,365]*8)[:-3])    # normal calendar, 2nd year leap
            if End_Name[period][dsi][0:4] == '2099':
                calblock = calblock[:-1]    # Clip the last year if it doesn't exist

        
    
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
        # fnout = '/Data/Projects/EMULATE/Processed/Future/' + 'ws_fut_' + Period_Name[period] + '_' + DS_Name[dsi] + '.dpkl'
        fnout =  # Path removed in online shared codes for security + 'ws_fut_' + Period_Name[period] + '_' + DS_Name[dsi] + '.dpkl'
        with open(fnout, 'wb') as f:
            dill.dump({'lat':lat, 'lon':lon, 'annual_max':annual_max, 'exceeds':exceeds, 'exc_lambda':exc_lambda, 'threshold':threshold, 'exc_time':exc_time},f)


        # Clean up
        dataset.close()
        del dataset, wsm, lat, lon, annual_max, exceeds, calblock, cnt_peaks, exc_lambda, win1, win2, winlen, vec, peaks, exc_time, peak_time, time



