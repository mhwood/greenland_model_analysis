

import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from pyproj import Transformer
from datetime import datetime, timedelta


def get_glacier_subglacial_discharge_timeseries(project_folder,glacier):

    file_path = os.path.join(project_folder,'Data','Modeling','MAR', glacier+' Subglacial Discharge.nc')

    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    sgd = ds.variables['runoff'][:]
    ds.close()

    timeseries = np.column_stack([time,sgd])

    annual_integrated_timeseries = np.copy(timeseries)
    years = np.arange(1992,2020)
    days_counted = 0
    for year in years:
        if year%4==0:
            n_days = 366
        else:
            n_days = 365
        annual_integrated_timeseries[days_counted+1:days_counted+n_days,1] = \
            np.cumsum(timeseries[days_counted+1:days_counted+n_days,1])*24*60*60
        days_counted+=n_days

    return(timeseries, annual_integrated_timeseries)

def get_glacier_melt_timeseries(project_folder,glacier):

    file_path = os.path.join(project_folder, 'Data', 'Modeling', 'Downscaled', 'L3_Scoresby_Sund', 'L3_Scoresby_Sund_DJG_Melt_Timeseries.nc')

    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    melt = ds.variables['melt'][:]
    ds.close()

    # fix an issue with the time steps
    time_step = np.mean(np.diff(time)[:10])
    new_time = np.copy(time)
    for t in range(len(time)):
        if time[t]>2020:
            new_time[t]=time[t-1]+time_step
            time[t] = time[t - 1] + time_step
    time = new_time

    # print(np.diff(time)[:10])
    # plt.plot(np.diff(new_time))
    # plt.show()

    timeseries = np.column_stack([time,melt])

    timeseries = timeseries[melt!=0,:]

    indices = time<2004
    mean_melt = np.mean(timeseries[indices,1])
    # mean_melt = np.mean(timeseries[:, 1])
    anomaly_timeseries = np.copy(timeseries)
    anomaly_timeseries[:,1] -= mean_melt

    timestep = 0.25 # days

    cumulative_melt_rate_anomaly = np.copy(timeseries)
    cumulative_melt_rate_anomaly[:,1] = 0
    for t in range(np.shape(cumulative_melt_rate_anomaly)[0]-1):
        cumulative_melt_rate_anomaly[t+1,1] = cumulative_melt_rate_anomaly[t,1] + anomaly_timeseries[t+1,1]*timestep

    # years = np.arange(1992,2022)
    # cumulative_melt_timeseries = np.zeros((len(years),2))
    # cumulative_melt_timeseries[:,0] = years+0.5
    # for y in range(len(years)):
    #     year = years[y]
    #     indices = np.logical_and(time>=year,time<year+1)
    #     if np.any(indices):
    #         cumulative_melt_timeseries[y,1] = np.sum(timeseries[indices,1])*timestep

    # annual_integrated_timeseries = np.copy(timeseries)
    #
    # days_counted = 0
    # year_days_counted = 0
    # for t in range(np.shape(timeseries)[0]):
    #
    #     if time[t]-int(time[t]) <= 1/(4*365):
    #         year_days_counted = 0
    #         # print(time[t],time[t]-int(time[t]))
    #
    #     if year_days_counted==0:
    #         annual_integrated_timeseries[days_counted,1] = 0
    #     else:
    #         annual_integrated_timeseries[days_counted,1] = annual_integrated_timeseries[days_counted-1,1] + timeseries[days_counted,1]*0.25
    #
    #     days_counted+=1
    #     year_days_counted+=1

    # plt.plot(annual_integrated_timeseries[:,0],annual_integrated_timeseries[:,1])
    # plt.show()

    return(timeseries, cumulative_melt_rate_anomaly)

def get_velocity_timeseries(project_folder,glacier,velocity_source_names):

    file_path = os.path.join(project_folder,'Data','Remote Sensing','Velocity',
                             glacier+' Velocity Timeseries.nc')

    velocity_timeseries = []

    ds = nc4.Dataset(file_path)

    # 'MEaSUREs-ITS_LIVE', 'MEaSUREs-InSAR'

    for source in velocity_source_names:

        if source=='MEaSUREs-ITS_LIVE':
            grp = ds.groups['ITS_LIVE']
        if source=='MEaSUREs-InSAR':
            grp = ds.groups['MEaSUREs']

        # grp = ds.groups[source]
        time = grp.variables['time'][:]
        vel = grp.variables['velocity'][:]
        timeseries = np.column_stack([time, vel])
        velocity_timeseries.append(timeseries)

    ds.close()

    return(velocity_timeseries)

def get_glacier_retreat_timeseries(project_folder,glacier):

    file_path = os.path.join(project_folder,'Data','Remote Sensing','Retreat',
                             glacier+' Retreat Timeseries.nc')

    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    retreat = ds.variables['retreat'][:]
    ds.close()

    timeseries = np.column_stack([time,retreat])

    return(timeseries)

def get_elevation_timeseries(project_folder,glacier,elevation_source_names):

    file_path = os.path.join(project_folder,'Data','Remote Sensing','Elevation',
                             glacier+' Elevation Timeseries.nc')

    elevation_timeseries = []

    ds = nc4.Dataset(file_path)

    for source in elevation_source_names:
        grp = ds.groups[source]
        time = grp.variables['time'][:]
        vel = grp.variables['elevation'][:]
        timeseries = np.column_stack([time, vel])
        timeseries = timeseries[np.logical_and(timeseries[:,1]<500,timeseries[:,1]>220),:]
        elevation_timeseries.append(timeseries)


    ds.close()

    return(elevation_timeseries)


def plot_comparison_figure(output_file, retreat_timeseries,
                           all_velocity_timeseries, velocity_source_names,
                           all_elevation_timeseries, elevation_source_names,
                           melt_timeseries, cumulative_melt_timeseries, for_publication=False):

    dpi=300
    fig = plt.figure(figsize=(8, 10), dpi=dpi)

    ############################################################################################################
    # Melt Timeseries

    plt.subplot(5, 1, 1)

    # plt.plot(melt_timeseries[:,0],melt_timeseries[:,1], linewidth = 0.7, color='b')

    plt.plot(melt_timeseries[:, 0], melt_timeseries[:, 1], color='steelblue')

    plt.text(melt_timeseries[0, 0], melt_timeseries[0, 1],
             'Start of L2 Model $\\rightarrow$ ', ha='right',va='bottom')

    ymin = np.min(melt_timeseries[:, 1])
    ymax = np.max(melt_timeseries[:, 1])
    y_range = ymax - ymin
    ymin -= y_range * 0.1
    ymax += y_range * 0.1
    plt.plot([2005, 2005], [ymin, ymax], 'k--')
    plt.plot([2011, 2011], [ymin, ymax], 'k--')
    plt.gca().set_ylim([ymin, ymax])

    plt.gca().set_xlim([1985, 2025])
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('Melt Rate\n(m/day)')
    plt.gca().set_xticklabels([])
    plt.title('Daugaard-Jensen Modeled Melt vs Observed Glacier Dynamics')
    plt.text(1985.4, np.max(melt_timeseries[:,1]), 'a) ', ha='left', va='top', fontsize=11)

    ############################################################################################################
    # Cumulative Melt Timeseries

    plt.subplot(5, 1, 2)

    # plt.plot(melt_timeseries[:,0],melt_timeseries[:,1], linewidth = 0.7, color='b')

    plt.plot(cumulative_melt_timeseries[:, 0], cumulative_melt_timeseries[:, 1], color='steelblue')

    ymin = np.min(cumulative_melt_timeseries[:, 1])
    ymax = np.max(cumulative_melt_timeseries[:, 1])
    y_range = ymax - ymin
    ymin -= y_range * 0.1
    ymax += y_range * 0.1
    plt.plot([2005, 2005], [ymin, ymax], 'k--')
    plt.plot([2011, 2011], [ymin, ymax], 'k--')
    plt.gca().set_ylim([ymin, ymax])

    # ax2 = plt.gca().twinx()
    # ax2.plot(annual_integrated_melt_timeseries[:,0],annual_integrated_melt_timeseries[:,1],'b-')

    plt.gca().set_xlim([1985, 2025])

    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('Integrated Melt Anomaly\n(m, relative to\n 2000-2004)')
    plt.gca().set_xticklabels([])
    plt.text(1985.4, np.max(cumulative_melt_timeseries[:,1]),
             'b) ', ha='left', va='top', fontsize=11)

    ############################################################################################################
    # Velocity Timeseries

    plt.subplot(5,1,3)
    colors = ['steelblue','olivedrab','k']
    symbols = ['.','s','.']
    for vt in range(len(all_velocity_timeseries)):
        velocity_timeseries = all_velocity_timeseries[vt]
        plt.plot(velocity_timeseries[:, 0]+0.5, velocity_timeseries[:, 1],
                 color=colors[vt], marker=symbols[vt], label=velocity_source_names[vt],linewidth=0)
        for i in range(np.shape(velocity_timeseries)[0]):
            plt.plot([velocity_timeseries[i, 0], velocity_timeseries[i, 0]+1],
                     [velocity_timeseries[i, 1], velocity_timeseries[i, 1]],'-',linewidth=0.5, color=colors[vt])
            plt.plot([velocity_timeseries[i, 0]+0.5, velocity_timeseries[i, 0] + 0.5],
                     [velocity_timeseries[i, 1]-30, velocity_timeseries[i, 1]+30], '-', linewidth=0.5, color=colors[vt])

    total_velocity_timeseries = np.vstack(all_velocity_timeseries)
    ymin = np.min(total_velocity_timeseries[:, 1])
    ymax = np.max(total_velocity_timeseries[:, 1])
    y_range = ymax - ymin
    ymin -= y_range * 0.1
    ymax += y_range * 0.1
    plt.plot([2005, 2005], [ymin, ymax], 'k--')
    plt.plot([2011, 2011], [ymin, ymax], 'k--')
    plt.gca().set_ylim([ymin, ymax])

    plt.text(2005, ymax-0.1*y_range, 'Start of Ice Acceleration $\\rightarrow$ ',
             ha='right', va='top')

    plt.ylabel('Speed (m/yr)')

    plt.gca().set_xlim([1985,2025])
    plt.grid(linestyle='--',alpha=0.5)
    plt.legend(loc=6)
    plt.gca().set_xticklabels([])

    plt.text(1985.4, ymax-y_range*0.1,
             'c) ', ha='left', va='top', fontsize=11)

    ############################################################################################################
    # Elevation Timeseries
    plt.subplot(5, 1, 4)
    colors = ['steelblue', 'k', 'olivedrab','purple','olivedrab','k']
    symbols = ['s', '^', '.','s','.','s']
    markersizes = [8, 5, 10, 8, 10, 8]

    total_elevation_timeseries = np.vstack(all_elevation_timeseries)
    total_elevation_timeseries_subset = total_elevation_timeseries[total_elevation_timeseries[:,0]>2009,:]

    p1 = np.poly1d(np.polyfit(total_elevation_timeseries_subset[:,0],
                               total_elevation_timeseries_subset[:,1], 1))

    threshold = 10

    filtered_elevation_timeseries = []
    for timeseries in all_elevation_timeseries:
        if np.any(timeseries[:,0]>2010):
            diff = np.abs(timeseries[:,1] - p1(timeseries[:,0]))
            timeseries = timeseries[diff<threshold,:]
        filtered_elevation_timeseries.append(timeseries)


    for et in range(len(filtered_elevation_timeseries)):
        elevation_timeseries = filtered_elevation_timeseries[et]
        plt.plot(elevation_timeseries[:, 0], elevation_timeseries[:, 1],
                 color=colors[et], marker=symbols[et], markersize = markersizes[et],
                 label=elevation_source_names[et], linewidth=0)


    ymin = np.min(total_elevation_timeseries[:, 1])
    ymax = np.max(total_elevation_timeseries[:, 1])
    y_range = ymax - ymin
    ymin -= y_range * 0.1
    ymax += y_range * 0.1
    plt.plot([2005, 2005], [ymin, ymax], 'k--')
    plt.plot([2011, 2011], [ymin, ymax], 'k--')
    plt.gca().set_ylim([ymin, ymax])

    plt.gca().set_xlim([1985, 2025])
    plt.grid(linestyle='--', alpha=0.5)
    plt.legend(loc=3, ncol=2)
    plt.ylabel('Ice Elevation (m)')
    plt.gca().set_xticklabels([])

    plt.text(2011, ymax - y_range * 0.1, ' $\\leftarrow$ Start of Ice Thinning\n       and Retreat',
             ha='left', va='top')

    plt.text(1985.4, ymax - y_range * 0.1,
             'd) ', ha='left', va='top', fontsize=11)


    # ############################################################################################################
    # # Subglacial Discharge Timeseries
    # plt.subplot(5, 1, 3)
    # plt.plot(sgd_timeseries[:,0],sgd_timeseries[:,1], linewidth = 0.7, color='k')
    #
    # # ax2 = plt.gca().twinx()
    # # ax2.plot(annual_integrated_sgd_timeseries[:,0],annual_integrated_sgd_timeseries[:,1],'k-')
    #
    #
    # plt.gca().set_xlim([1985, 2022])
    # plt.grid(linestyle='--', alpha=0.5)
    # plt.ylabel('Subglacial\nDischarge (m$^3$/s)')
    # plt.gca().set_xticklabels([])

    ############################################################################################################
    # Retreat Timeseries
    plt.subplot(5, 1, 5)

    y_shift = np.mean(retreat_timeseries[retreat_timeseries[:,0]<2010,1])
    plt.plot(retreat_timeseries[:, 0], (retreat_timeseries[:, 1]-y_shift)/1000, 'k^', label='TermPicks',markersize=5)

    ymin = np.min((retreat_timeseries[:, 1]-y_shift)/1000)
    ymax = np.max((retreat_timeseries[:, 1]-y_shift)/1000)
    y_range = ymax - ymin
    ymin -= y_range * 0.1
    ymax += y_range * 0.1
    plt.plot([2005, 2005], [ymin, ymax], 'k--')
    plt.plot([2011, 2011], [ymin, ymax], 'k--')
    plt.gca().set_ylim([ymin, ymax])

    plt.gca().set_xlim([1985, 2025])
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('Ice Front Position\n(km, relative to\n1985-2010 mean)')
    plt.legend(loc=4, ncol=2)

    plt.text(1985.4, ymax - y_range * 0.1,
             'e) ', ha='left', va='top', fontsize=11)

    if for_publication:
        plt.savefig(output_file[:-3]+'pdf',bbox_inches='tight', dpi=dpi)
    else:
        plt.savefig(output_file, bbox_inches='tight', dpi=dpi)
    plt.close(fig)


glacier = 'Daugaard Jensen'

project_folder = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund'
data_folder = '/Users/mhwood/Documents/Research/Projects/Glacier Variability'

longitude = -28.65249
latitude = 71.90093

velocity_source_names = ['MEaSUREs-ITS_LIVE','MEaSUREs-InSAR']
all_velocity_timeseries = get_velocity_timeseries(project_folder,glacier,velocity_source_names)

elevation_source_names = ['GLISTIN','ArcticDEM','IceBridge','Orthophoto','ICESat2','GIMP']
all_elevation_timeseries = get_elevation_timeseries(project_folder,glacier,elevation_source_names)

retreat_timeseries = get_glacier_retreat_timeseries(project_folder,glacier)

# sgd_timeseries, annual_integrated_sgd_timeseries = get_glacier_subglacial_discharge_timeseries(project_folder,glacier)

melt_timeseries, cumulative_melt_timeseries = get_glacier_melt_timeseries(project_folder,glacier)

output_file = os.path.join(project_folder,'Figures','Glacier','_'.join(glacier.split())+'_Dynamics_Timeseries.png')
plot_comparison_figure(output_file, retreat_timeseries,
                       all_velocity_timeseries, velocity_source_names,
                       all_elevation_timeseries, elevation_source_names,
                       melt_timeseries, cumulative_melt_timeseries, for_publication=True)















