
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import interp1d

def collect_timeseries_from_profile_files(model_dir,subsets):

    integration_depths = np.arange(200,450)

    ################################################################################################
    # Get the ECCOv5 timeseries

    ds = nc4.Dataset(os.path.join(model_dir, 'profile_timeseries', 'ECCOv5_THETA_Sound_Entrance_profile_timeseries.nc'))
    time = ds.variables['time'][:]
    depth = np.abs(ds.variables['RC'][:])
    theta = ds.variables['THETA'][:]
    theta_timeseries = np.zeros_like(time)
    for timestep in range(len(time)):
        set_int = interp1d(depth, theta[:, timestep])
        theta_timeseries[timestep] = np.mean(set_int(integration_depths))

    ################################################################################################
    # Get the ECCOv5 extension timeseries

    ds = nc4.Dataset(os.path.join(model_dir, 'profile_timeseries',
                                  'ECCOv5_extension_THETA_Sound_Entrance_profile_timeseries.nc'))
    time_extension = ds.variables['time'][:]
    depth = np.abs(ds.variables['RC'][:])
    theta = ds.variables['THETA'][:]
    theta_timeseries_extension = np.zeros_like(time_extension)
    for timestep in range(len(time_extension)):
        set_int = interp1d(depth, theta[:, timestep])
        theta_timeseries_extension[timestep] = np.mean(set_int(integration_depths))

    ################################################################################################
    # Put the timeseries together and output

    total_time = np.concatenate([time,time_extension])
    theta_timeseries = np.concatenate([theta_timeseries,theta_timeseries_extension])
    theta_timeseries_set = [theta_timeseries]

    return(total_time,theta_timeseries_set)

def store_temperauture_timeseries_as_nc(output_dir,output_file,total_time,subsets,theta_timeseries_set):

    ds = nc4.Dataset(os.path.join(output_dir,output_file),'w')
    ds.createDimension('time',len(total_time))

    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = total_time

    for s in range(len(subsets)):
        grp = ds.createGroup(subsets[s])
        svar = grp.createVariable('THETA','f4',('time',))
        svar[:] = theta_timeseries_set[s]

    ds.close()



ecco_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/Data/Modeling/ECCO'
output_file = 'ECCOv5_Scoresby_Sund_Timeseries.nc'

subsets = ['sound_entrance']

total_time,theta_timeseries_set = collect_timeseries_from_profile_files(ecco_dir,subsets)

store_temperauture_timeseries_as_nc(ecco_dir,output_file,total_time,subsets,theta_timeseries_set)
