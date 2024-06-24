
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import interp1d

def collect_timeseries_from_CTD_results(model_dir,subsets):

    location_dict = {'shelf':1,
                     'sound_entrance': 2,
                     'mid_sound':7,
                     'mid_fjord':8,
                     'near_glacier':10}

    years = np.arange(2000,2008).tolist()

    timeseries_started = False

    n_timesteps = 0
    for year in years:
        if year%4==0 and year!=1992:
            n_timesteps += 366*24
        else:
            n_timesteps += 365*24

    total_time = np.zeros((n_timesteps,))
    theta_timeseries_set = []
    for s in range(len(subsets)):
        theta_timeseries_set.append(np.zeros((n_timesteps,)))

    integration_depths = np.arange(200,401)

    timesteps_counted = 0
    for year in years:
        file_name = 'L3_CE_CTD_profiles_'+str(year)+'.nc'
        if file_name in os.listdir(model_dir):
            print('    - Collecting data from '+file_name)
            ds = nc4.Dataset(os.path.join(model_dir,file_name))
            time = ds.variables['time'][:]
            depth = ds.variables['depth'][:]
            total_time[timesteps_counted:len(time)+timesteps_counted] = time
            for s in range(len(subsets)):
                point_number = location_dict[subsets[s]]
                grp = ds.groups['point_'+'{:02d}'.format(point_number)]
                theta = grp.variables['THETA'][:,:]
                for timestep in range(len(time)):
                    set_int = interp1d(depth,theta[:,timestep])
                    theta_timeseries_set[s][timesteps_counted+timestep] = np.mean(set_int(integration_depths))
                # theta_timeseries_set[s][timesteps_counted:len(time)+timesteps_counted] = theta[18,:]
        timesteps_counted += len(time)

    indices = np.logical_and(theta_timeseries_set[0] != 0, total_time!=0)
    total_time = total_time[indices]
    for s in range(len(subsets)):
        theta_timeseries_set[s] = theta_timeseries_set[s][indices]

    # plt.plot(total_time,theta_timeseries_set[0])
    # plt.show()

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

def collect_profiles_from_CTD_results(model_dir,subsets):

    location_dict = {'shelf':1,
                     'sound_entrance': 2,
                     'mid_sound':7,
                     'mid_fjord':8,
                     'near_glacier':10}

    years = np.arange(2000,2008).tolist()

    timeseries_started = False

    n_timesteps = 0
    for year in years:
        if year%4==0 and year!=1992:
            n_timesteps += 366*24
        else:
            n_timesteps += 365*24

    Nr = 67

    total_time = np.zeros((n_timesteps,))
    theta_profiles_set = []
    for s in range(len(subsets)):
        theta_profiles_set.append(np.zeros((Nr, n_timesteps)))

    integration_depths = np.arange(200,401)

    timesteps_counted = 0
    for year in years:
        file_name = 'L3_CE_CTD_profiles_'+str(year)+'.nc'
        if file_name in os.listdir(model_dir):
            print('    - Collecting data from '+file_name)
            ds = nc4.Dataset(os.path.join(model_dir,file_name))
            time = ds.variables['time'][:]
            depth = ds.variables['depth'][:]
            total_time[timesteps_counted:len(time)+timesteps_counted] = time
            for s in range(len(subsets)):
                point_number = location_dict[subsets[s]]
                grp = ds.groups['point_'+'{:02d}'.format(point_number)]
                theta = grp.variables['THETA'][:,:]
                theta_profiles_set[s][:, timesteps_counted:len(time)+timesteps_counted] = theta
        timesteps_counted += len(time)

    indices = np.logical_and(theta_profiles_set[0][0,:] != 0, total_time!=0)
    total_time = total_time[indices]
    for s in range(len(subsets)):
        theta_profiles_set[s] = theta_profiles_set[s][:,indices]

    # plt.plot(total_time,theta_profiles_set[0])
    # plt.show()

    return(total_time,depth,theta_profiles_set)

def store_temperauture_profiles_as_nc(output_dir,output_file,total_time,depth,subsets,theta_timeseries_set):

    ds = nc4.Dataset(os.path.join(output_dir,output_file),'w')
    ds.createDimension('time',len(total_time))
    ds.createDimension('depth', len(depth))

    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = total_time

    dvar = ds.createVariable('depth', 'f4', ('depth',))
    dvar[:] = depth

    for s in range(len(subsets)):
        grp = ds.createGroup(subsets[s])
        svar = grp.createVariable('THETA','f4',('depth','time'))
        svar[:,:] = theta_timeseries_set[s]

    ds.close()



output_dir = '/Users/michwood/Documents/Research/Projects/Scoresby Sund/Data/Modeling/Downscaled/L3_Scoresby_Sund'


model_dir = '/Volumes/mhwood/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
            'downscaled_greenland/L3/L3_Scoresby_Sund/results_iceplume_ks22/CTD'

subsets = ['shelf','sound_entrance','mid_sound','mid_fjord','near_glacier']
# subsets = ['sound_entrance']

# timeseries first
output_file = 'L3_Scoresby_Sund_dv_CTD_Timeseries.nc'
total_time,theta_timeseries_set = collect_timeseries_from_CTD_results(model_dir,subsets)
store_temperauture_timeseries_as_nc(output_dir,output_file,total_time,subsets,theta_timeseries_set)

# then profiles
output_file = 'L3_Scoresby_Sund_dv_CTD_Profiles.nc'
total_time,depth,theta_profiles_set = collect_profiles_from_CTD_results(model_dir,subsets)
store_temperauture_profiles_as_nc(output_dir,output_file,total_time,depth,subsets,theta_profiles_set)

