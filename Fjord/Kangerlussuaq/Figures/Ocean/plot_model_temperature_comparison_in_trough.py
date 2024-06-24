

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime
from matplotlib.gridspec import GridSpec
from scipy.interpolate import griddata

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)


def iter_number_to_dec_yr(iter_number,seconds_per_iter=60):

    total_seconds = iter_number*seconds_per_iter
    date = datetime.datetime(1992,1,1) + datetime.timedelta(seconds=total_seconds)

    dec_yr = YMD_to_DecYr(date.year,date.month,date.day)
    # print(date)
    return(dec_yr)


def read_L1_model_results_form_nc(project_dir):

    output_file = os.path.join(project_dir, 'Data', 'Modeling', 'L1_state_at_L2_CTDs.nc')
    ds = nc4.Dataset(output_file)

    depths = ds.variables['depths'][:]
    iterations = ds.variables['iterations'][:]
    longitudes = ds.variables['longitude'][:]
    latitudes = ds.variables['latitude'][:]
    Theta = ds.variables['Theta'][:, :, :]

    ds.close()

    dec_yrs = np.zeros((len(iterations),))
    for i in range(len(iterations)):
        dec_yrs[i] = iter_number_to_dec_yr(iterations[i],seconds_per_iter=300)

    return(depths, dec_yrs, longitudes, latitudes, Theta)


def plot_model_temperature_comparison(project_dir, depths, dec_yrs, longitudes, latitudes, Theta):

    min_year = 2015
    max_year = 2022

    fig = plt.figure(figsize=(8, 10))

    location_number = 8
    # for i in range(len(longitudes)):
    #     print(i+1, longitudes[i], latitudes[i])

    gs = GridSpec(3, 1)

    ax1 = fig.add_subplot(gs[0, 0])
    C = ax1.pcolormesh(dec_yrs, depths, Theta[location_number-1,:,:].T, cmap='turbo', shading='nearest')
    plt.colorbar(C)
    ax1.set_xlim([min_year, max_year])
    max_depth_index = np.sum(Theta[location_number-1,10,:]!=0)
    max_depth = depths[max_depth_index]
    ax1.set_ylim([max_depth, 0])
    ax1.set_ylabel('L1 Model\nDepth(m')
    ax1.set_title('Model Temperature Comparison (Lon: '+
                  '{:.4f}'.format(longitudes[location_number-1])+', Lat: '+
                  '{:.4f}'.format(latitudes[location_number-1])+')')

    ax2 = fig.add_subplot(gs[1, 0])
    # C = ax1.pcolormesh(dec_yrs, depths, Theta[location_number - 1, :, :].T, cmap='turbo', shading='nearest')
    # plt.colorbar(C)
    ax2.set_xlim([min_year, max_year])
    ax2.set_ylim([max_depth, 0])
    ax2.set_ylabel('L2 Model\nDepth(m')

    aw_depth = 350
    Dec_yrs, Depth = np.meshgrid(dec_yrs,depths)
    points = np.column_stack([Dec_yrs.ravel(), Depth.ravel()])
    values = np.ravel(Theta[location_number-1,:,:].T)
    L1_AW_theta = griddata(points, values, (dec_yrs,aw_depth))


    ax3 = fig.add_subplot(gs[2, 0])
    ax3.plot(dec_yrs[dec_yrs>=min_year], L1_AW_theta[dec_yrs>=min_year])
    ax3.set_xlim([min_year, max_year])


    output_file = os.path.join(project_dir, 'Figures', 'Modeling', 'Kangerlussuaq Trough Model Temperature.png')
    plt.savefig(output_file)
    plt.close(fig)

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'


depths, dec_yrs, longitudes, latitudes, Theta = read_L1_model_results_form_nc(project_dir)

plot_model_temperature_comparison(project_dir, depths, dec_yrs, longitudes, latitudes, Theta)


