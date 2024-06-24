



import os
import argparse
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import cmocean.cm as cm
from matplotlib.gridspec import GridSpec
from scipy.interpolate import interp1d
import datetime


def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)


def calculate_isopycnal_timeseries(depth,density,ref_density,temperature, temp_density):

    isopycnal_timeseries = np.zeros((np.shape(density)[0], np.shape(density)[2]))
    temp_timeseries = np.zeros((np.shape(density)[0], np.shape(density)[2]))

    for t in range(np.shape(isopycnal_timeseries)[0]):
        if t%10==0:
            print('    - Calculating step '+str(t)+' of '+str(np.shape(isopycnal_timeseries)[0]))
        for i in range(np.shape(isopycnal_timeseries)[1]):
            density_profile = density[t,:,i]
            temperature_profile = temperature[t, :, i]
            indices = density_profile>1005 # filter for the bathy
            if np.min(density_profile)<ref_density and np.max(density_profile)>ref_density:
                set_int = interp1d(density_profile[indices],depth[indices])
                isopycnal_timeseries[t, i] = set_int(ref_density)
            else:
                isopycnal_timeseries[t,i] = np.NaN
            if np.min(density_profile)<temp_density and np.max(density_profile)>temp_density:
                # set_int = interp1d(density_profile[indices], depth[indices])
                temp_set_int = interp1d(density_profile[indices],temperature_profile[indices])
                temp_timeseries[t, i] = temp_set_int(temp_density)
            else:
                temp_timeseries[t,i] = np.NaN

    return(isopycnal_timeseries, temp_timeseries)


def plot_transect_timeseries(project_dir,output_dir,for_publication=False):

    model_level = 'L3'
    model_name = 'L3_Scoresby_Sund'
    testing = False
    plot_anomaly = True

    transect_file = os.path.join(project_dir,'Data','Modeling','Downscaled',
                               'L3_Scoresby_Sund','Scoresby_Sund_transects.nc')

    ds = nc4.Dataset(transect_file)
    years = ds.variables['years'][:]
    months = ds.variables['months'][:]
    depth = ds.variables['depth'][:]
    distance = ds.variables['transect_distance'][:]
    density = ds.variables['transect_density'][:, :, :]
    temperature = ds.variables['transect_temp'][:, :, :]
    ds.close()

    distance *= 1e-3

    time = np.zeros(np.shape(density)[0])
    for t in range(len(time)):
        time[t] = YMD_to_DecYr(int(years[t]), int(months[t]), 15)

    if testing:
        time=time[:24]
        density = density[:24,:,:]
        temperature = temperature[:24, : ,:]

    ref_density = 1028.5
    temp_density = 1029
    isopycnal_timeseries, temp_timeseries = calculate_isopycnal_timeseries(depth, density, ref_density, temperature, temp_density)

    if for_publication:
        fig = plt.figure(figsize=(10,10),dpi=300)
    else:
        fig = plt.figure(figsize=(10, 10))

    gs = GridSpec(21, 11, left=0.1, right=0.9, bottom = 0.05, top=0.95)

    #####################################################################
    # Isopycnals

    ax1 = fig.add_subplot(gs[:4, :])
    for t in range(len(time)):
        ax1.plot(distance, isopycnal_timeseries[t, :], '-', color='silver', linewidth=0.8)
    ax1.plot(distance, np.mean(isopycnal_timeseries,axis=0),'k-')
    ax1.set_xlabel('Distance Along Transect')
    ax1.invert_yaxis()
    ax1.set_title('Isopycnal Range for $\\rho = $' + '{:.2f}'.format(ref_density) + ' g/cm$^3$ on Along-Fjord Transect')

    #####################################################################
    # Isopycnals

    ax2 = fig.add_subplot(gs[6:-2, 6:])
    if not plot_anomaly:
        C = ax2.pcolormesh(distance,time,isopycnal_timeseries, shading='nearest')
        cbar = plt.colorbar(C)
        cbar.ax.invert_yaxis()
    else:
        isopycnal_timeseries_anomaly = np.copy(isopycnal_timeseries)
        for t in range(np.shape(isopycnal_timeseries)[0]):
            isopycnal_timeseries_anomaly[t,:] -= np.mean(isopycnal_timeseries, axis=0)
        anomaly_max = np.max(np.abs(isopycnal_timeseries_anomaly[~np.isnan(isopycnal_timeseries_anomaly)]))
        C = ax2.pcolormesh(distance, time, isopycnal_timeseries_anomaly,
                           shading='nearest', cmap='seismic',
                           vmin=-1*anomaly_max,vmax=anomaly_max)
        # cbar = plt.colorbar(C)
        # cbar.ax.invert_yaxis()
        # cbar.set_label('Depth Anomaly (m)\n($\leftarrow$ Anomalously Deep   to   Anomalously Shallow $\\rightarrow$)')
    ax2.set_yticks(np.arange(int(np.min(time)) + 1, int(np.max(time)) + 1))
    ax2.set_xlabel('Distance Into Fjord')
    ax2.set_title('Isopycnal Depth Anomaly\n($\\rho = $' + '{:.2f}'.format(ref_density) + ' g/cm$^3$)')

    ax2c = fig.add_subplot(gs[-1, 6:])
    x = np.linspace(-anomaly_max, anomaly_max, 100)
    y = np.arange(0, 1.1, 0.1)
    X, Y = np.meshgrid(x, y)
    ax2c.pcolormesh(X, Y, X, cmap='seismic')
    ax2c.set_yticklabels([])
    ax2c.set_xlabel('Isopycnal Depth Anomaly (m)')

    #####################################################################
    # Temperature

    ax3 = fig.add_subplot(gs[6:-2, :5])
    plot_depth = 300
    depth_index = np.argmin(np.abs(depth - plot_depth))
    if not plot_anomaly:
        C = ax3.pcolormesh(distance, time, temp_timeseries, shading='nearest')
        cbar = plt.colorbar(C)
        cbar.ax.invert_yaxis()
    else:
        temperature_anomaly = np.copy(temp_timeseries)
        for t in range(np.shape(temperature_anomaly)[0]):
            temperature_anomaly[t, :] -= np.mean(temp_timeseries, axis=0)
        temp_anomaly_max = np.max(np.abs(temperature_anomaly[~np.isnan(temperature_anomaly)]))
        C = ax3.pcolormesh(distance, time, temperature_anomaly,
                           shading='nearest', cmap='seismic',
                           vmin=-1 * temp_anomaly_max, vmax = temp_anomaly_max)
        # cbar = plt.colorbar(C, orientation='horizontal')
        # cbar.set_label('Temperature Anomaly ($^{\circ}$C)')
    ax3.set_yticks(np.arange(int(np.min(time))+1, int(np.max(time))+1))
    ax3.set_xlabel('Distance Into Fjord')
    ax3.set_title('Temperature Anomaly\n($\\rho = $' + '{:.2f}'.format(temp_density) + ' g/cm$^3$)')

    ax3c = fig.add_subplot(gs[-1, :5])
    x = np.linspace(-temp_anomaly_max,temp_anomaly_max,100)
    y = np.arange(0,1.1,0.1)
    X, Y = np.meshgrid(x,y)
    ax3c.pcolormesh(X,Y,X,cmap='seismic')
    ax3c.set_yticklabels([])
    ax3c.set_xlabel('Temperature Anomaly ($^{\circ}$C)')


    if for_publication:
        output_file = os.path.join(output_dir,'L3_isopycnal_transect_timeseries.pdf')
        plt.savefig(output_file, dpi=300)
    else:
        output_file = os.path.join(output_dir, 'L3_isopycnal_transect_timeseries.png')
        plt.savefig(output_file)
    plt.close(fig)





config_dir = '/Volumes/helheim/Ocean_Modeling/Projects/Downscale_Greenland/' \
             'MITgcm/configurations/downscaled_greenland'
project_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund'
output_dir = project_dir+'/Figures/Ocean'

plot_transect_timeseries(project_dir, output_dir, for_publication = True)

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#
#     parser.add_argument("-d", "--config_dir", action="store",
#                         help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
#                         type=str, required=True)
#
#     args = parser.parse_args()
#     config_dir = args.config_dir





