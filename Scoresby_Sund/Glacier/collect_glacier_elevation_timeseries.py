
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4

def get_glacier_elevation_timeseries(data_folder,glacier,source='ITS_LIVE Annual Mosaic Elevation Grids'):

    if source=='ICESat2 Elevation Points':
        timeseries_started = False
        folder = os.path.join(data_folder,'Data','Glaciers',glacier,'Elevation','Timeseries','Sample Point','Sources')
        for cycle in range(1,12):
            file_name = glacier+' '+source+' - Cycle '+str(cycle)+' Sample Point Elevation Timeseries.csv'
            if file_name in os.listdir(folder):
                cycle_timeseries = np.genfromtxt(os.path.join(folder,file_name), delimiter=',')
                if np.size(cycle_timeseries) == 2:
                    cycle_timeseries = np.reshape(cycle_timeseries, (1, 2))
                if not timeseries_started:
                    timeseries = cycle_timeseries
                    timeseries_started = True
                else:
                    timeseries = np.vstack([timeseries,cycle_timeseries])

    else:
        file_path = os.path.join(data_folder,'Data','Glaciers',glacier,'Elevation','Timeseries','Sample Point','Sources',
                                 glacier+' '+source+' Sample Point Elevation Timeseries.csv')
        timeseries = np.genfromtxt(file_path,delimiter=',')

    if np.size(timeseries)==2:
        timeseries = np.reshape(timeseries,(1,2))
    print(np.shape(timeseries))

    return(timeseries)

def save_timeseries_as_nc(output_file,sources,all_timeseries):

    ds = nc4.Dataset(output_file,'w')

    for s in range(len(sources)):
        if 'Orthophoto' in sources[s]:
            source_out = 'Orthophoto'
        else:
            source_out = sources[s].split()[0]
        grp = ds.createGroup(source_out)
        grp.createDimension('time',np.shape(all_timeseries[s])[0])
        tvar = grp.createVariable('time','f4',('time',))
        tvar[:] = all_timeseries[s][:,0]
        vvar = grp.createVariable('elevation', 'f4', ('time',))
        vvar[:] = all_timeseries[s][:, 1]

    ds.close()



data_folder = '/Volumes/upernavik/Research/Projects/Glacier Variability'

glacier = 'Daugaard Jensen'

sources = ['ArcticDEM Elevation Grids','GLISTIN Elevation Grids',
           'IceBridge ATM Elevation Points','Korsgaard Orthophoto Elevation Grid',
           'ICESat2 Elevation Points','GIMP Elevation Grid']

all_timeseries = []
for source in sources:
    timeseries = get_glacier_elevation_timeseries(data_folder,glacier,source)
    all_timeseries.append(timeseries)


output_file = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/Data/Remote Sensing/Elevation/' \
              'Daugaard Jensen Elevation Timeseries.nc'
save_timeseries_as_nc(output_file,sources,all_timeseries)





