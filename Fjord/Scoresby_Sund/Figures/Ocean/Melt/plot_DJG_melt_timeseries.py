

import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from pyproj import Transformer
from datetime import datetime, timedelta


def get_glacier_melt_timeseries(project_folder,glacier):

    file_path = os.path.join(project_folder, 'Data', 'Modeling', 'Downscaled', 'L3_Scoresby_Sund', 'L3_Scoresby_Sund_DJG_Melt_Timeseries.nc')

    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    melt = ds.variables['melt'][:]
    ds.close()

    timeseries = np.column_stack([time,melt])

    timeseries = timeseries[melt!=0,:]

    mean_melt = np.mean(timeseries[:,1])
    anomaly_timeseries = np.copy(timeseries)
    anomaly_timeseries[:,1] -= mean_melt

    timestep = 0.25 # days

    # cumulative_melt_rate_anomaly = np.copy(timeseries)
    # cumulative_melt_rate_anomaly[:,1] = 0
    # for t in range(np.shape(cumulative_melt_rate_anomaly)[0]-1):
    #     cumulative_melt_rate_anomaly[t+1,1] = cumulative_melt_rate_anomaly[t,1] + anomaly_timeseries[t+1,1]*timestep

    years = np.arange(1992,2022)
    cumulative_melt_timeseries = np.zeros((len(years),2))
    cumulative_melt_timeseries[:,0] = years+0.5
    for y in range(len(years)):
        year = years[y]
        indices = np.logical_and(time>=year,time<year+1)
        if np.any(indices):
            cumulative_melt_timeseries[y,1] = np.sum(timeseries[indices,1])*timestep

    years = np.arange(1992, 2022)
    cumulative_summer_melt_timeseries = np.zeros((len(years), 2))
    cumulative_summer_melt_timeseries[:, 0] = years + 0.5
    for y in range(len(years)):
        year = years[y]
        indices = np.logical_and(time >= year, time < year + 1)
        if np.any(indices):
            baseline = np.min(timeseries[indices, 1])
            print(baseline)
            cumulative_summer_melt_timeseries[y, 1] = np.sum(timeseries[indices, 1]-baseline) * timestep

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

    return(timeseries, cumulative_melt_timeseries, cumulative_summer_melt_timeseries)


def plot_comparison_figure(output_file, melt_timeseries, cumulative_melt_timeseries, cumulative_summer_melt_timeseries):

    fig = plt.figure(figsize=(8, 10))


    ############################################################################################################
    # Velocity Timeseries

    plt.subplot(3,1,1)
    plt.plot(melt_timeseries[:,0],melt_timeseries[:,1], linewidth = 0.7)

    plt.ylabel('Maximal Melt Rate (m/day)')
    plt.title('Daugaard-Jensen Glacier Ice Front Melt')

    plt.gca().set_xlim([1999,2006])
    plt.grid(linestyle='--',alpha=0.5)
    plt.gca().set_xticklabels([])

    ############################################################################################################
    # Elevation Timeseries
    plt.subplot(3, 1, 2)
    plt.bar(cumulative_melt_timeseries[:, 0], cumulative_melt_timeseries[:, 1])

    plt.gca().set_xlim([1999,2006])
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('Cumulative Annual Melt (m)')
    plt.gca().set_xticklabels([])

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
    plt.subplot(3, 1, 3)

    # plt.plot(melt_timeseries[:,0],melt_timeseries[:,1], linewidth = 0.7, color='b')

    plt.bar(cumulative_summer_melt_timeseries[:, 0], cumulative_summer_melt_timeseries[:, 1])

    # ax2 = plt.gca().twinx()
    # ax2.plot(annual_integrated_melt_timeseries[:,0],annual_integrated_melt_timeseries[:,1],'b-')

    plt.gca().set_xlim([1999,2006])
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('Cumulative\nPlume-Induced Melt\n(m)')

    plt.savefig(output_file,bbox_inches='tight')
    plt.close(fig)


glacier = 'Daugaard Jensen'

project_folder = '/Users/michwood/Documents/Research/Projects/Scoresby Sund'

melt_timeseries, cumulative_melt_timeseries, cumulative_summer_melt_timeseries = get_glacier_melt_timeseries(project_folder,glacier)

output_file = os.path.join(project_folder,'Figures','Glacier','_'.join(glacier.split())+'_Melt_Timeseries.png')
plot_comparison_figure(output_file, melt_timeseries, cumulative_melt_timeseries, cumulative_summer_melt_timeseries)















