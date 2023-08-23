'''
dwpd_plot_gpd_ws_hist_multiplot.py
Loads and plots maps of return values for the models.
Plots all 52 downscalings on one plot in grid with gaps.

Written by Stephen Outten August 2022
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import dill
from mpl_toolkits.basemap import Basemap

spx = 11 
spy = 6
iprint = False
itime = 30   # return period is years
imFolderName = # Path removed in online shared codes for security  
FolderName = # Path removed in online shared codes for security
DS_Name = 'CLMcom_CNRM', 'CLMcom_ICHEC', 'CLMcom_MOHC', 'CLMcom_MPI', 'CNRM_CNRM', 'CNRM_MOHC', 'CNRM_MPI', 'CNRM_NorESM', 'DMI_CNRM', 'DMI_ICHEC', 'DMI_IPSL', 'DMI_MOHC', 'DMI_MPI', 'DMI_NorESM', 'ETH_CNRM', 'ETH_ICHEC', 'ETH_MOHC', 'ETH_MPI', 'ETH_NorESM', 'GERICS_CNRM', 'GERICS_IPSL', 'GERICS_MPI', 'GERICS_NorESM', 'ICTP_CNRM', 'ICTP_ICHEC', 'ICTP_MOHC', 'ICTP_MPI', 'ICTP_NorESM', 'IPSL_CNRM', 'IPSL_ICHEC', 'IPSL_IPSL', 'IPSL_MOHC', 'IPSL_MPI', 'IPSL_NorESM', 'KNMI_CNRM', 'KNMI_ICHEC', 'KNMI_IPSL', 'KNMI_MOHC', 'KNMI_MPI', 'KNMI_NorESM', 'MOHC_CNRM', 'MOHC_ICHEC', 'MOHC_MOHC', 'MOHC_MPI', 'MOHC_NorESM', 'MPI_MPI', 'SMHI_CNRM', 'SMHI_ICHEC', 'SMHI_IPSL', 'SMHI_MOHC', 'SMHI_MPI', 'SMHI_NorESM'
# Number for location of downscaling in grid
# Grid_no = np.array([0, 11, 33, 44, 1, 34, 45, 56, 2, 13, 24, 35, 46, 57, 3, 14, 36, 47, 58, 4, 26, 48, 59, 5, 16, 38, 49, 60, 6, 17, 28, 39, 50, 61, 7, 18, 29, 40, 51, 62, 8, 19, 41, 52, 63, 53, 10, 21, 32, 43, 54, 65])
DS_no = np.array([0, 4, 8, 14, 19, 23, 28, 34, 40, -1, 46, 1, -1, 9, 15, -1, 24, 29, 35, 41, -1, 47, -1, -1, 10, -1, 20, -1, 30, 36, -1, -1, 48, 2, 5, 11, 16, -1, 25, 31, 37, 42, -1, 49, 3, 6, 12, 17, 21, 26, 32, 38, 43, 45, 50, -1, 7, 13, 18, 22, 27, 33, 39, 44, -1, 51])
imFileName1 = 'gpd_ws_hist_multiplot_R30_full.png'


# Set plotting style to semi-classic
matplotlib.style.use('classic')   # restores original matplotlib style after update to matplotlib 2.0
matplotlib.rcParamsOrig
matplotlib.rcParams.update({'font.size': 8})   # adjust plot font size
matplotlib.rcParams.update({'savefig.orientation': 'landscape'})
matplotlib.rcParams.update({'ps.papersize': 'a4'})
matplotlib.rcParams.update({'lines.linewidth': 0.7})
matplotlib.rcParams.update({'axes.titlepad': 1.5})
matplotlib.rcParams.update({'xtick.major.pad': 1.5})
matplotlib.rcParams.update({'ytick.major.pad': 1.5})
matplotlib.rcParams.update({'ytick.major.size': 3.0})
matplotlib.rcParams.update({'xtick.major.size': 3.0})
matplotlib.rcParams.update({'ytick.minor.size': 1.5})
matplotlib.rcParams.update({'xtick.minor.size': 1.5})
matplotlib.rcParams.update({'axes.labelpad': 1.0})


# PLot the map of return level
f1, axs = plt.subplots(spy, spx, sharex='col', sharey='row')
axs = axs.ravel()      # stops axs from being 4x5 and instead makes it 20 in size
f1.subplots_adjust(hspace = .17, wspace=.17)

# Loading grid data
FileName = DS_Name[0] + '_ws_hist_gpd_maps.dpkl'
fn = FolderName + FileName
lat = dill.load(open(fn,'rb'))['lat']
lon = dill.load(open(fn,'rb'))['lon']

# Identify corners of central area for plotting
# This avoids showing the curved domain
llat = np.max(lat[0,:])   # max lat of bottom curve
llon = lon[:,0][lat[:,0]<llat][-1]   # lon on left side which contains minlat
ulon = lon[:,-1][lat[:,-1]<llat][-1]   # lon on right side which contains minlat
ullat = lat[-1,:][lon[-1,:]<llon][-1]   # lat on top which contains lllon
urlat = lat[-1,:][lon[-1,:]>ulon][0]   # lat on top which contains lrlon
ulat = np.min([ullat,urlat])

'''
for dsi in range(len(DS_Name)):
    # Load data from one of the models
    print(DS_Name[dsi])
    print('Loading...')
    FileName = DS_Name[dsi] + '_ws_hist_gpd_maps.dpkl'
    fn = FolderName + FileName
    Rx = dill.load(open(fn,'rb'))['GPD_R'+str(itime)]
    lat = dill.load(open(fn,'rb'))['lat']
    lon = dill.load(open(fn,'rb'))['lon']
    lon[lon>180] = lon[lon>180] - 360

    print('Plotting...')
    grdno = Grid_no[dsi]
    m = Basemap(projection='mill', llcrnrlat=llat,urcrnrlat=ulat,llcrnrlon=llon
            ,urcrnrlon=ulon,resolution='l',ax=axs[grdno],fix_aspect=False)
    m.drawcoastlines()
    x, y = m(lon, lat)
#    cs = m.pcolor(x,y,Rx,vmin=0,vmax=45)
#    cs = m.pcolor(x,y,Rx,vmin=0,vmax=45,cmap='twilight_shifted')
    cs = m.pcolor(x,y,Rx,vmin=0,vmax=45,cmap='afmhot_r')
    axs[grdno].set_title(DS_Name[dsi].replace('_',' '),fontsize='5')
    if grdno/spx == grdno//spx:
        m.drawparallels(np.arange(0,90,10),labels=[True,False,False,False], linewidth=0, fontsize='5')
    if (grdno//spx) == spy-1:
        m.drawmeridians(np.arange(-50,90,10),labels=[False,False,False,True], linewidth=0, fontsize='5')
'''

for pi in range(spx*spy):
    print('Plotting...')
    m = Basemap(projection='mill', llcrnrlat=llat,urcrnrlat=ulat,llcrnrlon=llon
            ,urcrnrlon=ulon,resolution='l',ax=axs[pi],fix_aspect=False)
    x, y = m(lon, lat)  

    # Load data from one of the models
    if DS_no[pi]>-1:
    # Load data from one of the models
        print(DS_Name[DS_no[pi]])
        print('Loading...')
        FileName = DS_Name[DS_no[pi]] + '_ws_hist_gpd_maps.dpkl'
        fn = FolderName + FileName
        Rx = dill.load(open(fn,'rb'))['GPD_R'+str(itime)]
        lt = dill.load(open(fn,'rb'))['lat']
        ln = dill.load(open(fn,'rb'))['lon']
        ln[ln>180] = ln[ln>180] - 360
        x, y = m(ln, lt)
        m.drawcoastlines(linewidth=0.5)

    #    cs = m.pcolor(x,y,Rx,vmin=0,vmax=45)
    #    cs = m.pcolor(x,y,Rx,vmin=0,vmax=45,cmap='twilight_shifted')
        cs = m.pcolor(x,y,Rx,vmin=0,vmax=45,cmap='afmhot_r')
        axs[pi].set_title(DS_Name[DS_no[pi]].replace('_',' '),fontsize='5')
    
    if pi/spx == pi//spx:
        m.drawparallels(np.arange(0,90,10),labels=[True,False,False,False], linewidth=0, fontsize='4')
    if (pi//spx) == spy-1:
        m.drawmeridians(np.arange(-45,90,15),labels=[False,False,False,True], linewidth=0, fontsize='4')


f1.subplots_adjust(right=0.9)
cbar_ax = f1.add_axes([0.93, 0.15, 0.02, 0.7])
f1.colorbar(cs, cax=cbar_ax)

# Save the plot
if iprint:
    print("Saving Figures")
    imfn = imFolderName + imFileName1 
    plt.savefig(imfn, format='png', bbox_inches='tight',dpi=600)
#    plt.savefig(imfn, orientation='portrait', format='png', bbox_inches='tight',dpi=600)


plt.show()
