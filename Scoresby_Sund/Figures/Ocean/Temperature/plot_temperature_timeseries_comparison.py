
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.gridspec import GridSpec
import scipy.io
from scipy.interpolate import interp1d

def get_hadley_objective_analysis_timeseries(project_folder):
    file_path = os.path.join(project_folder,'Data','In Situ','Hadley','Timeseries','Hadley_IPCC_Timeseries.nc')
    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    temp = ds.variables['Theta'][:]
    ds.close()
    timeseries = np.column_stack([time,temp])
    timeseries[timeseries[:,0]!=0,:]

    return(timeseries)

def get_regional_ECCO_temperature_timeseries(project_folder):
    file_path = os.path.join(project_folder,'Data','Modeling','ECCO','ECCOv5_THETA_IPCC_timeseries.nc')
    ds = nc4.Dataset(file_path)
    time = ds.variables['dec_yr'][:]
    temp = ds.variables['THETA_mean'][:]
    ds.close()
    timeseries = np.column_stack([time,temp])

    timeseries[np.floor(time).astype(int)==2018,1] = np.nan

    return(timeseries)

def get_shelf_ECCO_temperature_timeseries(project_folder):
    file_path = os.path.join(project_folder,'Data','Modeling','ECCO','ECCOv5_Scoresby_Sund_Timeseries.nc')
    ds = nc4.Dataset(file_path)
    grp = ds.groups['sound_entrance']
    time = ds.variables['time'][:]
    temp = grp.variables['THETA'][:]
    ds.close()
    timeseries = np.column_stack([time,temp])

    timeseries[np.floor(time).astype(int)==2018,1] = np.nan

    return(timeseries)

def get_regional_nested_temperature_timeseries(project_folder,model_name):
    file_path = os.path.join(project_folder,'Data','Modeling','Downscaled',model_name,model_name+'_IPCC_timeseries.nc')
    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    temp = ds.variables['Theta'][:]
    ds.close()
    timeseries = np.column_stack([time,temp])
    return(timeseries)

def get_nested_model_temperature_timeseries(project_folder,model_name,locations):

    retrieval_locations = list(locations)
    all_timeseries = []

    file_path = os.path.join(project_folder, 'Data', 'Modeling', 'Downscaled', model_name,
                             model_name + '_THETA_Timeseries.nc')
    print(file_path)
    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    for location_name in retrieval_locations:
        grp_name = '_'.join(location_name.lower().split())
        if grp_name in list(ds.groups.keys()):
            grp = ds.groups[grp_name]
            theta = grp.variables['THETA'][:]
            all_timeseries.append(np.column_stack([time, theta]))
        else:
            all_timeseries.append([])
    ds.close()

    return(time, all_timeseries)

def get_in_situ_temperature_data(project_folder,source,locations):

    all_timeseries = []

    for location_name in locations:

        file_path = os.path.join(project_folder,'Data','In Situ',source,'Data',source+'_'+location_name+'_Temperature_Timeseries.nc')

        if os.path.exists(file_path):

            ds = nc4.Dataset(file_path)
            time = ds.variables['time'][:]
            theta = ds.variables['temperature_mean'][:]
            theta_std = ds.variables['temperature_std'][:]
            ds.close()

            timeseries = np.column_stack([time,theta,theta_std])

        else:
            timeseries = []

        all_timeseries.append(timeseries)

    return(all_timeseries)

def calculate_RMSE(model_timeseries_set, obs_timeseries_set, output_names):
    RMSE_values = []
    for m in range(1,len(model_timeseries_set)):
        set_int = interp1d(np.array(model_timeseries_set[m][:,0]),np.array(model_timeseries_set[m][:,1]))
        errors = set_int(np.array(obs_timeseries_set[m][:,0])) - obs_timeseries_set[m][:,1]
        print(output_names[m],errors)
        RMSE = np.sqrt(np.mean(errors**2))
        # print(output_names[m],RMSE)
        RMSE_values.append(RMSE)
    return(RMSE_values)

def create_timeseries_plot(output_file, plot_titles, ecco_shelf_timeseries,
                           regional_ecco_timeseries, regonal_hadley_timeseries, regional_L1_timeseries,
                           L1_time, L1_temperature_timeseries, L1_locations, L1_RMSE,
                           L2_time, L2_temperature_timeseries, L2_locations, L2_RMSE,
                           omg_temperature_timeseries, omg_locations):

    fig = plt.figure(figsize=(8, 9))

    gs = GridSpec(4, 1, left=0.15, right=0.95, hspace=0.05, top=0.95,bottom=0.05)

    ####################################################################################################
    # regional ocean temperature

    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(regional_L1_timeseries[:,0],regional_L1_timeseries[:,1],'-',color='olivedrab',label='L1')
    ax1.plot(regional_ecco_timeseries[:,0],regional_ecco_timeseries[:,1],'-',color='purple',label='L0')
    # ax1.plot(regonal_hadley_timeseries[:, 0], regonal_hadley_timeseries[:, 1], 'k.', label='EN4.2.2 Objective Analysis')
    ax1.set_ylabel('a) Regional Ocean\nTemperature ($^{\circ}$C)\n(200-500 m)')
    ax1.set_xticklabels([])
    # ax1.set_ylim([0, 1])
    ax1.set_xlim([1991.8, 2022.1])
    ax1.grid(linestyle='--', alpha=0.4, linewidth=0.5)
    ax1.set_title('Atlantic Water Temperature from Model Framework')
    ax1.legend(loc=2,ncol=2)

    ####################################################################################################
    # shelf temperature

    ax2 = fig.add_subplot(gs[1, :])
    ax2.plot(L1_temperature_timeseries[L1_locations.index('Sound Entrance')][:, 0],
             L1_temperature_timeseries[L1_locations.index('Sound Entrance')][:, 1], '-',
             label='L1',color='olivedrab') # (RMSE: '+'{:.2f}'.format(L1_RMSE[0])+'$^{\circ}$C)
    ax2.plot(L2_temperature_timeseries[L2_locations.index('Sound Entrance')][:,0],
             L2_temperature_timeseries[L2_locations.index('Sound Entrance')][:,1], '-',
             label = 'L2',color='steelblue') # (RMSE: '+'{:.2f}'.format(L2_RMSE[0])+'$^{\circ}$C)
    ax2.plot(ecco_shelf_timeseries[:, 0],
             ecco_shelf_timeseries[:, 1], '-',
             label='L0', color='purple')
    # ax2.plot(omg_temperature_timeseries[omg_locations.index('Sound Entrance')][:, 0],
    #          omg_temperature_timeseries[omg_locations.index('Sound Entrance')][:, 1], 'k.', label='OMG CTDs')
    # ax2.set_ylabel('Sound Entrance\nTemperature ($^{\circ}$C)\n(200-450 m)')
    ax2.set_ylabel('b) Sound Entrance\nTemperature ($^{\circ}$C)\n(200-450 m)')
    ax2.legend(loc=2,ncol=2)
    ax2.set_xticklabels([])
    ax2.set_xlim([1991.8,2022.1])
    ax2.grid(linestyle='--',alpha=0.4,linewidth=0.5)

    ####################################################################################################
    # mod-sound temperature

    ax3 = fig.add_subplot(gs[2, :])
    ax3.plot(L1_temperature_timeseries[L1_locations.index('Mid-Sound')][:, 0],
             L1_temperature_timeseries[L1_locations.index('Mid-Sound')][:, 1], '-',
             label='L1',color='olivedrab') # (RMSE: '+'{:.2f}'.format(L1_RMSE[1])+'$^{\circ}$C)
    ax3.plot(L2_temperature_timeseries[L2_locations.index('Mid-Sound')][:, 0],
             L2_temperature_timeseries[L2_locations.index('Mid-Sound')][:, 1], '-',
             label='L2',color='steelblue') # (RMSE: '+'{:.2f}'.format(L2_RMSE[1])+'$^{\circ}$C)
    # ax3.plot(omg_temperature_timeseries[omg_locations.index('Mid-Sound')][:, 0],
    #          omg_temperature_timeseries[omg_locations.index('Mid-Sound')][:, 1], 'k.', label='OMG CTDs')
    ax3.set_ylabel('c) Mid-Sound\nTemperature ($^{\circ}$C)\n(200-400 m)')
    ax3.set_xticklabels([])
    ax3.set_xlim([1991.8, 2022.1])
    ax3.grid(linestyle='--', alpha=0.4, linewidth=0.5)
    ax3.legend(loc=2,ncols=2)

    ####################################################################################################
    # near-glacier temperature

    ax4 = fig.add_subplot(gs[3, :])
    ax4.plot(L2_temperature_timeseries[L2_locations.index('Near-DJG')][:, 0],
             L2_temperature_timeseries[L2_locations.index('Near-DJG')][:, 1], '-',
             label='L2',color='steelblue') # (RMSE: '+'{:.2f}'.format(L2_RMSE[2])+'$^{\circ}$C)
    # ax4.plot(omg_temperature_timeseries[omg_locations.index('Near-DJG')][:, 0],
    #          omg_temperature_timeseries[omg_locations.index('Near-DJG')][:, 1], 'k.', label='OMG CTDs')
    ax4.set_ylabel('d) Near-DJG\nTemperature ($^{\circ}$C)\n(200-450 m)')
    ax4.set_xlim([1991.8, 2022.1])
    ax4.grid(linestyle='--', alpha=0.4, linewidth=0.5)
    ax4.legend(loc=2)


    plt.savefig(output_file)
    plt.close(fig)


project_folder = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund'

plot_titles = ['Regional Ocean', 'Sound Entrance', 'Mid-Sound', 'Near-DJG']
output_names =['Regional_Ocean', 'Sound_Entrance', 'Mid_Sound', 'Near_Glacier']

print('    - Reading in the regional timeseries')
regional_L1_timeseries = get_regional_nested_temperature_timeseries(project_folder,'L1_CE_Greenland')
regonal_hadley_timeseries = get_hadley_objective_analysis_timeseries(project_folder)
regional_ecco_timeseries = get_regional_ECCO_temperature_timeseries(project_folder)
ecco_shelf_timeseries = get_shelf_ECCO_temperature_timeseries(project_folder)

print('    - Reading in OMG timeseries')
omg_locations = plot_titles
OMG_locations = ['']+plot_titles[1:]
OMG_output_names = ['']+output_names[1:]
omg_temperature_timeseries = get_in_situ_temperature_data(project_folder,'OMG',output_names)

print('    - Reading in the L1 timeseries')
L1_model_name = 'L1_CE_Greenland'
L1_locations = ['']+plot_titles[1:]
L1_output_names = ['']+output_names[1:]
L1_time, L1_temperature_timeseries = get_nested_model_temperature_timeseries(project_folder, L1_model_name, L1_output_names)
L1_RMSE = calculate_RMSE(L1_temperature_timeseries, omg_temperature_timeseries, L1_output_names)


print('    - Reading in the L2 timeseries')
L3_model_name = 'L3_Scoresby_Sund'
L2_locations = ['']+plot_titles[1:]
L2_output_names = ['']+output_names[1:]
L2_time, L2_temperature_timeseries = get_nested_model_temperature_timeseries(project_folder, L3_model_name, L2_output_names)
L2_RMSE = calculate_RMSE(L2_temperature_timeseries, omg_temperature_timeseries, L2_output_names)

output_file = os.path.join(project_folder,'Figures','Ocean','Timeseries Comparison',
                                                            'Scoresby_Sund_Temperature_Comparison.png')

print('    - Creating the plot')
create_timeseries_plot(output_file, plot_titles, ecco_shelf_timeseries,
                       regional_ecco_timeseries, regonal_hadley_timeseries,regional_L1_timeseries,
                       L1_time, L1_temperature_timeseries, L1_locations, L1_RMSE,
                       L2_time, L2_temperature_timeseries, L2_locations, L2_RMSE,
                       omg_temperature_timeseries, OMG_locations)