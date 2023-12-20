
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4


def collect_timeseries_from_CTD_results(model_dir,point_number):

    years = np.arange(1992,2022).tolist()

    n_timesteps = 0
    for year in years:
        if year%4==0 and year!=1992:
            n_timesteps += 366*24
        else:
            n_timesteps += 365*24

    total_time = np.zeros((n_timesteps,))
    theta_timeseries = np.zeros((50,n_timesteps))

    timesteps_counted = 0
    for year in years:
        file_name = 'L1_CE_CTD_profiles_'+str(year)+'.nc'
        if file_name in os.listdir(model_dir):
            print('    - Collecting data from '+file_name)
            ds = nc4.Dataset(os.path.join(model_dir,file_name))
            time = ds.variables['time'][:]
            depth = ds.variables['depth'][:]
            total_time[timesteps_counted:len(time)+timesteps_counted] = time
            grp = ds.groups['point_'+'{:02d}'.format(point_number)]
            theta = grp.variables['THETA'][:,:]
            theta_timeseries[:, timesteps_counted:timesteps_counted+len(time)] = theta
        timesteps_counted += len(time)

    indices = np.logical_and(theta_timeseries[0,:] != 0, total_time!=0)
    total_time = total_time[indices]
    theta_timeseries = theta_timeseries[:, indices]


    return(total_time, depth, theta_timeseries)


def plot_theta_timeseries(output_file, time, depth, temperature_sets):

    fig = plt.figure(figsize=(8,10))

    plt.subplot(3,1,1)
    C = plt.contourf(time,depth,temperature_sets[0],100,cmap='turbo',vmin=-1.9,vmax=1.5)
    plt.colorbar(C)
    plt.ylabel('Depth (m)')
    plt.gca().set_ylim([450,0])
    plt.title('Sound Entrance')

    plt.subplot(3, 1, 2)
    C = plt.contourf(time, depth, temperature_sets[1], 100, cmap='turbo', vmin=-1.9, vmax=1.5)
    plt.colorbar(C)
    plt.ylabel('Depth (m)')
    plt.gca().set_ylim([450, 0])
    plt.title('Mid-Sound')

    plt.subplot(3, 1, 3)
    C = plt.contourf(time, depth, temperature_sets[2], 100, cmap='turbo', vmin=-1.9, vmax=1.5)
    plt.colorbar(C)
    plt.ylabel('Depth (m)')
    plt.gca().set_ylim([1000, 0])
    plt.title('Near-Glacier')
    # plt.gca().invert_yaxis()

    plt.savefig(output_file)
    plt.close(fig)


project_folder = '/Users/michwood/Documents/Research/Projects/Scoresby Sund'

model_dir = '/Volumes/mhwood/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
            'downscaled_greenland/L1/grid/L1_CE_Greenland/results/CTD'

output_file = os.path.join(project_folder,'Figures','Scoresby_Sund_L1_CTD_Profile_Timeseries.png')

point_numbers = [1,7,10]

temperature_sets = []
for point_number in point_numbers:
    time, depth, theta_timeseries = collect_timeseries_from_CTD_results(model_dir,point_number)
    temperature_sets.append(theta_timeseries)

print()
plot_theta_timeseries(output_file, time, depth, temperature_sets)