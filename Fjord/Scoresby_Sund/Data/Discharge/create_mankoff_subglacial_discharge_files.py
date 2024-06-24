
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4


# step 1: find outlet ID points for a glacier
def find_glacier_ids_in_mankoff_data(glacier):
    if glacier == 'Daugaard Jensen':
        cats = [43939,43965,44065,44126]
    elif glacier == 'Jakobshavn Isbrae':
        cats = [61276, 61278, 61318, 61392, 61447, 61517, 61495, 61547, 61446, 61597, 61586, 61568]
    else:
        raise ValueError('Didnt finish this function for other glaciers yet')
    return(cats)


# step 2: read in the timeseries for the glacier from the Mankoff data
def read_subglacial_discharge_from_mankoff_data(mankoff_dir,glacier_ids,model='MAR'):

    years = np.arange(1992,2020)
    N = 0
    for year in years:
        if year%4==0:
            N+=366
        else:
            N+=365
    timeseries = np.zeros((N,2))

    days_counted = 0
    for year in years:
        print('    - ')
        file_path = os.path.join(mankoff_dir,'runoff','ice','runoff',model+'_'+str(year)+'.nc')
        ds = nc4.Dataset(file_path)
        station_ids = ds.variables['station'][:]
        # lon = ds.variables['lon'][:]
        # lat = ds.variables['lat'][:]
        time = ds.variables['time'][:]
        runoff = ds.variables['runoff'][:,:]
        year_time_series = np.zeros((len(time),2))
        year_time_series[:,0] = year+np.arange(len(time))/len(time)
        for id in glacier_ids:
            id_index = np.argmin(np.abs(station_ids-id))
            if np.sum(np.isnan(runoff[id_index, :])) == 0:
                year_time_series[:,1] += runoff[id_index,:]
        timeseries[days_counted:days_counted+len(time),:] = year_time_series
        days_counted += len(time)

    plt.plot(timeseries[:,0],timeseries[:,1])
    plt.show()

    return (timeseries)


# step 3: store timeseries as an nc file
def store_timeseries_as_nc(output_dir,glacier,timeseries):

    output_file = os.path.join(output_dir,glacier+' Subglacial Discharge.nc')
    ds = nc4.Dataset(output_file,'w')
    ds.createDimension('time',np.shape(timeseries)[0])
    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = timeseries[:,0]
    dvar = ds.createVariable('runoff', 'f4', ('time',))
    dvar[:] = timeseries[:, 1]
    ds.close()


mankoff_dir = '/Users/michwood/Documents/Research/Data Repository/Greenland/Runoff/Mankoff_liquid'

glacier = 'Daugaard Jensen'
glacier = 'Jakobshavn Isbrae'

output_dir = '/Users/michwood/Documents/Research/Projects/Scoresby Sund/Data/Modeling/MAR'

glacier_ids = find_glacier_ids_in_mankoff_data(glacier)

timeseries = read_subglacial_discharge_from_mankoff_data(mankoff_dir,glacier_ids)


store_timeseries_as_nc(output_dir,glacier,timeseries)


