


import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from datetime import datetime, timedelta


# step 1: find outlet ID points for a glacier
def read_solid_discharge_from_mankoff_data(mankoff_dir,glacier):

    ds = nc4.Dataset(os.path.join(mankoff_dir,'gate.nc'))
    names = ds.variables['name_Mouginot'][:].tolist()
    discharge = ds.variables['discharge'][:,:]
    error = ds.variables['err'][:,:]
    time = ds.variables['time'][:]
    ds.close()

    dec_yr = np.zeros((len(time),))
    for t in range(len(time)):
        date = datetime(1986,4,15)+timedelta(days=int(time[t]))
        yr = date.year
        yr_days = ((date-datetime(yr,1,1)).total_seconds())/(24*60*60)
        if yr%4==0:
            dec_yr[t] = yr + yr_days / 366
        else:
            dec_yr[t] = yr + yr_days / 365

    if glacier == 'Daugaard Jensen':
        index = names.index('DAUGAARD-JENSEN')
    elif glacier == 'Jakobshavn Isbrae':
        index = names.index('JAKOBSHAVN_ISBRAE')
    else:
        raise ValueError('Didnt finish this function for other glaciers yet')

    timeseries = np.column_stack([dec_yr.ravel(),discharge[index,:].ravel(),error[index,:].ravel()])

    return(timeseries)


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

    output_file = os.path.join(output_dir,glacier+' Solid Ice Discharge.nc')
    ds = nc4.Dataset(output_file,'w')
    ds.createDimension('time',np.shape(timeseries)[0])
    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = timeseries[:,0]
    dvar = ds.createVariable('discharge', 'f4', ('time',))
    dvar[:] = timeseries[:, 1]
    evar = ds.createVariable('error', 'f4', ('time',))
    evar[:] = timeseries[:, 2]
    ds.close()


mankoff_dir = '/Users/michwood/Documents/Research/Data Repository/Greenland/Discharge/Mankoff'

glacier = 'Daugaard Jensen'
# glacier = 'Jakobshavn Isbrae'

output_dir = '/Users/michwood/Documents/Research/Projects/Scoresby Sund/Data/Remote Sensing/Discharge'

timeseries = read_solid_discharge_from_mankoff_data(mankoff_dir,glacier)

plt.plot(timeseries[:,0],timeseries[:,1],'k.')
for t in range(len(timeseries)):
    plt.plot([timeseries[t,0],timeseries[t,0]],
             [timeseries[t,1]+0.5*timeseries[t,2],timeseries[t,1]-0.5*timeseries[t,2]],'k-')

plt.ylabel('Gt/yr')
plt.title('Jakobshavn Discharge')
plt.show()


store_timeseries_as_nc(output_dir,glacier,timeseries)


