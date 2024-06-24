import os
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import datetime

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def time_to_dec_yr(time):

    dec_yr = np.zeros((len(time),))
    for t in range(len(time)):
        date = datetime.datetime(1950,1,1)+datetime.timedelta(days=int(time[t]))
        dec_yr[t] = YMD_to_DecYr(date.year, date.month, date.day)

    return(dec_yr)

def get_Qsg_timeseries(project_folder):
    file_path = os.path.join(project_folder,'Data','Modeling','MAR','Glacier Discharge Timeseries','Glacier_118_discharge.nc')
    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    land_discharge = ds.variables['land_discharge'][:]
    ice_discharge = ds.variables['ice_discharge'][:]
    ds.close()
    discharge = ice_discharge + land_discharge

    dec_yr = time_to_dec_yr(time)

    return(dec_yr,discharge)

def get_dv_theta_timeseries(model_folder,years):

    N = 0
    for year in years:
        if year % 4 == 0:
            N += 366 * 4
        else:
            N += 365 * 4
    Nr = 67

    full_time = np.zeros((N,))
    theta_timeseries = np.zeros((Nr,N))

    counter = 0
    for year in years:
        ds = nc4.Dataset(os.path.join(model_folder,'results_iceplume','CTD','L3_CE_CTD_profiles_' + str(year) + '.nc'))
        time = ds.variables['time'][:]
        grp = ds.groups['point_10']
        theta = grp.variables['THETA'][:, :]
        ds.close()

        full_time[counter:counter + len(time)] = time
        theta_timeseries[:, counter:counter + len(time)] = theta
        counter += len(time)

    indices = theta_timeseries[20,:]!=0

    full_time = full_time[indices]
    theta_timeseries = theta_timeseries[:,indices]

    return(full_time, theta_timeseries)

def get_dv_melt_rate_timeseries(model_folder,years):

    N = 0
    for year in years:
        if year % 4 == 0:
            N += 366 * 4
        else:
            N += 365 * 4
    Nr = 67

    full_time = np.zeros((N,))
    plume_timeseries = np.zeros((Nr,N))
    ambient_timeseries = np.zeros((Nr,N))

    counter = 0
    for year in years:
        ds = nc4.Dataset(os.path.join(model_folder,'results_iceplume','iceplume','L3_CE_iceplume_profiles_' + str(year) + '.nc'))
        time = ds.variables['time'][:]
        plume_melt = ds.variables['ICEFRNTM'][:, :, :]
        ambient_melt = ds.variables['ICEFRNTA'][:, :, :]
        depth = ds.variables['depth'][:]
        # C = plt.imshow(plume_melt[:,:,84],cmap='turbo')
        # plt.show()
        plume_melt = plume_melt[:, :, 84]
        ambient_melt = ambient_melt[:, :, 84]
        # plt.plot(time,plume_melt)
        # plt.show()
        ds.close()

        full_time[counter:counter + len(time)] = time
        plume_timeseries[:, counter:counter + len(time)] = plume_melt
        ambient_timeseries[:, counter:counter + len(time)] = ambient_melt
        counter += len(time)

    indices = ambient_timeseries[20,:]!=0

    full_time = full_time[indices]
    plume_timeseries = plume_timeseries[:,indices]
    ambient_timeseries = ambient_timeseries[:,indices]

    return(full_time, depth, plume_timeseries, ambient_timeseries)

def plot_melt_rate_components(project_folder, depth,
                              discharge_time, discharge,
                              theta_time, theta_timeseries,
                              melt_time, plume_melt_timeseries, ambient_melt_timeseries):

    fig = plt.figure(figsize=(10,12))

    gs = GridSpec(4, 11, left=0.1, right=0.90, hspace=0.05)

    ax1 = fig.add_subplot(gs[0, :-1])
    ax1.plot(discharge_time,discharge)
    ax1.set_xlim([2000,2002])
    ax1.set_ylim([0,600])
    plt.ylabel('Subglacial Discharge\n(m$^3$/s)')
    ax1.set_xticks(np.arange(2000, 2002, 1 / 12))
    plt.gca().set_xticklabels([])
    plt.grid(linewidth=0.5,alpha=0.5,linestyle='--')
    plt.title('Melt Rate Components for Daugaard-Jensen Glacier')


    ax2 = fig.add_subplot(gs[1, :-1])
    cax2 = fig.add_subplot(gs[1, -1])
    C = ax2.pcolormesh(theta_time, depth, theta_timeseries, cmap='RdYlBu_r', shading='nearest')
    ax2.grid(linewidth=0.5, alpha=0.5, linestyle='--')
    cbar = plt.colorbar(C,ax=cax2)
    cbar.ax.set_ylabel('Potential Temperature\n($^{\circ}$C)')
    cbar.ax.yaxis.set_label_position('left')
    cax2.axis('off')
    ax2.set_ylim([530, 0])
    ax2.set_xlim([2000, 2002])
    ax2.set_ylabel('Depth (m)')
    ax2.set_xticks(np.arange(2000, 2002, 1 / 12))
    ax2.set_xticklabels([])

    ax3 = fig.add_subplot(gs[2, :-1])
    cax3 = fig.add_subplot(gs[2, -1])
    C = ax3.pcolormesh(melt_time, depth, plume_melt_timeseries, cmap='turbo', shading='nearest')
    ax3.grid(linewidth=0.5, alpha=0.5, linestyle='--')
    cbar = plt.colorbar(C,ax=cax3)
    cbar.ax.set_ylabel('Plume Melt Rate\n(m/day)')
    cbar.ax.yaxis.set_label_position('left')
    cax3.axis('off')
    ax3.set_ylim([530,0])
    ax3.set_xlim([2000, 2002])
    ax3.set_ylabel('Depth (m)')
    ax3.set_xticks(np.arange(2000, 2002, 1 / 12))
    ax3.set_xticklabels([])

    ax4 = fig.add_subplot(gs[3, :-1])
    cax4 = fig.add_subplot(gs[3, -1])
    C = ax4.pcolormesh(melt_time, depth, ambient_melt_timeseries, cmap='turbo', shading='nearest')
    ax4.grid(linewidth=0.5, alpha=0.5, linestyle='--')
    cbar = plt.colorbar(C, ax=cax4)
    cbar.ax.set_ylabel('Ambient Melt Rate\n(m/day)')
    cbar.ax.yaxis.set_label_position('left')
    cax4.axis('off')
    ax4.set_ylim([530, 0])
    ax4.set_xlim([2000, 2002])
    ax4.set_ylabel('Depth (m)')
    ax4.set_xticks(np.arange(2000, 2002, 1 / 12))
    ax4.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'])

    plt.savefig(os.path.join(project_folder,'Figures','Ocean','DJG Melt Rate Components.png'),bbox_inches='tight')
    plt.close(fig)




project_folder = '/Users/michwood/Documents/Research/Projects/Scoresby Sund'

model_folder = '/Volumes/mhwood/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
            'downscaled_greenland/L3/L3_Scoresby_Sund'

years = [2000, 2001]

discharge_time, discharge = get_Qsg_timeseries(project_folder)

theta_time, theta_timeseries = get_dv_theta_timeseries(model_folder,years)

melt_time, depth, plume_melt_timeseries, ambient_melt_timeseries = get_dv_melt_rate_timeseries(model_folder,years)

plot_melt_rate_components(project_folder, depth,
                          discharge_time, discharge,
                          theta_time, theta_timeseries,
                          melt_time, plume_melt_timeseries, ambient_melt_timeseries)