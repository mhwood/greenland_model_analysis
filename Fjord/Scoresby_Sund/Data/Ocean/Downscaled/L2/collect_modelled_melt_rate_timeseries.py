





import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import interp1d

def collect_timeseries_from_iceplume_results(model_dir):

    years = np.arange(2000,2022).tolist()

    timeseries_started = False

    n_timesteps = 0
    for year in years:
        if year%4==0 and year!=1992:
            n_timesteps += 366*24
        else:
            n_timesteps += 365*24

    total_time = np.zeros((n_timesteps,))
    melt_timeseries = np.zeros((n_timesteps,))

    djg_number = 84

    timesteps_counted = 0
    for year in years:
        file_name = 'L3_CE_iceplume_profiles_'+str(year)+'.nc'
        if file_name in os.listdir(model_dir):
            print('    - Collecting data from '+file_name)
            ds = nc4.Dataset(os.path.join(model_dir,file_name))
            time = ds.variables['time'][:]
            total_time[timesteps_counted:len(time)+timesteps_counted] = time
            melt = ds.variables['ICEFRNTM'][:,:,:]
            melt = melt[:,:,djg_number]
            melt = np.max(melt,axis=0)
            melt_timeseries[timesteps_counted:len(time) + timesteps_counted] = melt
            timesteps_counted += len(time)

    indices = np.logical_and(melt_timeseries != 0, total_time!=0)
    total_time = total_time[indices]
    melt_timeseries = melt_timeseries[indices]

    # plt.plot(total_time,melt_timeseries)
    # plt.show()

    return(total_time,melt_timeseries)

def store_melt_timeseries_as_nc(output_dir,output_file,total_time,melt_timeseries):

    ds = nc4.Dataset(os.path.join(output_dir,output_file),'w')
    ds.createDimension('time',len(total_time))

    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = total_time

    mvar = ds.createVariable('melt', 'f4', ('time',))
    mvar[:] = melt_timeseries

    ds.close()



output_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/Data/Modeling/Downscaled/L3_Scoresby_Sund'
output_file = 'L3_Scoresby_Sund_DJG_Melt_Timeseries.nc'

model_dir = '/Volumes/Helheim/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
            'downscaled_greenland/L3/L3_Scoresby_Sund/results_iceplume_ks22/iceplume'


total_time,melt_timeseries = collect_timeseries_from_iceplume_results(model_dir)

store_melt_timeseries_as_nc(output_dir,output_file,total_time,melt_timeseries)
