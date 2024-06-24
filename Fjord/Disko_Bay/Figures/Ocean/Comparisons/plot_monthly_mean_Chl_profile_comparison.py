


import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
import datetime

def iter_number_to_date(iter_number):
    seconds_per_iter = 30

    total_seconds = iter_number*seconds_per_iter
    date = datetime.datetime(1992,1,1) + datetime.timedelta(seconds=total_seconds)
    return(date)

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_timeseries_from_nc(project_dir, var_name, year):

    output_file = os.path.join(project_dir,'Data','Modeling','L2',var_name+'_'+str(year)+'_mean_profile_timeseries.nc')
    ds = nc4.Dataset(output_file)

    grid = ds.variables[var_name][:, :]
    depth = ds.variables['depth'][:]
    iterations = ds.variables['iterations'][:]

    ds.close()

    dec_yrs = np.zeros((len(iterations),))
    for d in range(len(dec_yrs)):
        date = iter_number_to_date(iterations[d])
        dec_yrs[d] = YMD_to_DecYr(date.year, date.month, date.day, date.hour, date.minute)

    return(depth, dec_yrs, grid)

def compute_ticks_and_labels(year):
    month_labels = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    x_ticks = []
    x_grid_locations = []
    for month in range(1, 13):
        x_ticks.append(YMD_to_DecYr(year, month, 15))
        if month in [1, 3, 5, 7, 8, 10, 12]:
            x_grid_locations.append(YMD_to_DecYr(year, month, 31))
        elif month in [4, 6, 9, 11]:
            x_grid_locations.append(YMD_to_DecYr(year, month, 30))
        else:
            if year % 4 == 0:
                x_grid_locations.append(YMD_to_DecYr(year, month, 29))
            else:
                x_grid_locations.append(YMD_to_DecYr(year, month, 28))
    return(x_ticks, month_labels, x_grid_locations)

def plot_grids_with_differences(project_dir, var_name, depth,
                                dec_yrs_2008, timeseries_2008,
                                dec_yrs_2012, timeseries_2012,
                                dec_yrs_2017, timeseries_2017,
                                dec_yrs_2019, timeseries_2019):

    fig = plt.figure(figsize=(14, 10))

    image_rows = 6
    image_cols = 6

    max_depth = 100

    gs1 = GridSpec(image_rows * 4 + 3, image_cols * 2 + 1,
                   left=0.1, right=0.95, bottom=0.05, top=0.95)

    ########################################################################
    # 2008 timeseries

    x_ticks_2008, x_tick_labels_2008, x_grid_locations_2008 = compute_ticks_and_labels(year=2008)

    plot_row = 0
    ax21m = fig.add_subplot(gs1[plot_row * image_rows + 1:(plot_row + 1) * image_rows + 1, :image_cols])
    # ax21m.imshow(np.flipud(landsat_image),
    #              extent=[np.min(landsat_X) / 1000, np.max(landsat_X) / 1000, np.min(landsat_Y) / 1000,
    #                      np.max(landsat_Y) / 1000],
    #              alpha=0.8)
    ax21m.pcolormesh(dec_yrs_2008, depth, timeseries_2008, cmap='turbo', vmin=0, vmax=5, shading='nearest')
    # ax21m.set_xlim([np.min(polygon[:, 0] / 1000), np.max(polygon[:, 0] / 1000)])
    # ax21m.set_ylim([np.min(polygon[:, 1] / 1000), np.max(polygon[:, 1] / 1000)])
    # ax21m.set_yticks([])
    ax21m.set_ylim([max_depth, 0])
    ax21m.set_xticks(x_ticks_2008)
    ax21m.set_xticklabels(x_tick_labels_2008)
    for loc in range(len(x_grid_locations_2008)):
        ax21m.plot([x_grid_locations_2008[loc], x_grid_locations_2008[loc]], [max_depth, 0], '-', linewidth=0.5,
                   color='silver')

    ########################################################################
    # 2012 timeseries

    x_ticks_2012, x_tick_labels_2012, x_grid_locations_2012 = compute_ticks_and_labels(year=2012)

    plot_row = 1
    ax21m = fig.add_subplot(gs1[plot_row * image_rows + 1:(plot_row + 1) * image_rows + 1, :image_cols])
    # ax21m.imshow(np.flipud(landsat_image),
    #              extent=[np.min(landsat_X) / 1000, np.max(landsat_X) / 1000, np.min(landsat_Y) / 1000,
    #                      np.max(landsat_Y) / 1000],
    #              alpha=0.8)
    ax21m.pcolormesh(dec_yrs_2012, depth, timeseries_2012, cmap='turbo', vmin=0, vmax=5, shading='nearest')
    # ax21m.set_xlim([np.min(polygon[:, 0] / 1000), np.max(polygon[:, 0] / 1000)])
    # ax21m.set_ylim([np.min(polygon[:, 1] / 1000), np.max(polygon[:, 1] / 1000)])
    # ax21m.set_yticks([])
    ax21m.set_ylim([max_depth, 0])
    ax21m.set_xticks(x_ticks_2012)
    ax21m.set_xticklabels(x_tick_labels_2012)
    for loc in range(len(x_grid_locations_2012)):
        ax21m.plot([x_grid_locations_2012[loc], x_grid_locations_2012[loc]], [max_depth, 0], '-', linewidth=0.5,
                   color='silver')

    ########################################################################
    # 2017 timeseries

    x_ticks_2017, x_tick_labels_2017, x_grid_locations_2017 = compute_ticks_and_labels(year=2017)

    plot_row = 2
    ax21m = fig.add_subplot(gs1[plot_row * image_rows + 1:(plot_row+1) * image_rows + 1, :image_cols])
    # ax21m.imshow(np.flipud(landsat_image),
    #              extent=[np.min(landsat_X) / 1000, np.max(landsat_X) / 1000, np.min(landsat_Y) / 1000,
    #                      np.max(landsat_Y) / 1000],
    #              alpha=0.8)
    ax21m.pcolormesh(dec_yrs_2017, depth, timeseries_2017, cmap='turbo', vmin=0, vmax=5, shading='nearest')
    # ax21m.set_xlim([np.min(polygon[:, 0] / 1000), np.max(polygon[:, 0] / 1000)])
    # ax21m.set_ylim([np.min(polygon[:, 1] / 1000), np.max(polygon[:, 1] / 1000)])
    # ax21m.set_yticks([])
    ax21m.set_ylim([max_depth,0])
    ax21m.set_xticks(x_ticks_2017)
    ax21m.set_xticklabels(x_tick_labels_2017)
    for loc in range(len(x_grid_locations_2017)):
        ax21m.plot([x_grid_locations_2017[loc], x_grid_locations_2017[loc]], [max_depth, 0], '-', linewidth=0.5, color='silver')

    # ax21m.text(np.max(polygon[:, 0] / 1000) - 1, np.max(polygon[:, 1] / 1000) - 1, '2016/03/23', ha='right', va='top',
    #            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    ########################################################################
    # 2017 timeseries

    x_ticks_2019, x_tick_labels_2019, x_grid_locations_2019 = compute_ticks_and_labels(year=2019)

    plot_row = 3
    ax21m = fig.add_subplot(gs1[plot_row * image_rows + 1:(plot_row + 1) * image_rows + 1, :image_cols])
    # ax21m.imshow(np.flipud(landsat_image),
    #              extent=[np.min(landsat_X) / 1000, np.max(landsat_X) / 1000, np.min(landsat_Y) / 1000,
    #                      np.max(landsat_Y) / 1000],
    #              alpha=0.8)
    ax21m.pcolormesh(dec_yrs_2019, depth, timeseries_2019, cmap='turbo', vmin=0, vmax=5, shading='nearest')
    # ax21m.set_xlim([np.min(polygon[:, 0] / 1000), np.max(polygon[:, 0] / 1000)])
    # ax21m.set_ylim([np.min(polygon[:, 1] / 1000), np.max(polygon[:, 1] / 1000)])
    # ax21m.set_yticks([])
    ax21m.set_ylim([max_depth, 0])
    ax21m.set_xticks(x_ticks_2019)
    ax21m.set_xticklabels(x_tick_labels_2019)
    for loc in range(len(x_grid_locations_2019)):
        ax21m.plot([x_grid_locations_2019[loc], x_grid_locations_2019[loc]], [max_depth, 0], '-', linewidth=0.5,
                   color='silver')

        # ax21m.text(np.max(polygon[:, 0] / 1000) - 1, np.max(polygon[:, 1] / 1000) - 1, '2016/03/23', ha='right', va='top',
        #            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    output_file = os.path.join(project_dir,'Figures','Ocean','Modeling','L2','L2_Disko_Bay_'+var_name+'_Profile_Comparison.png')
    plt.savefig(output_file)
    plt.close()


config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/' \
             'configurations/downscale_darwin'

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Disko Bay'

model_name = 'L2_Disko_Bay'
var_name = 'Total_Chl'

depth, dec_yrs_2008, timeseries_2008 = read_timeseries_from_nc(project_dir, var_name, year=2008)
depth, dec_yrs_2012, timeseries_2012 = read_timeseries_from_nc(project_dir, var_name, year=2012)
_, dec_yrs_2017, timeseries_2017 = read_timeseries_from_nc(project_dir, var_name, year=2017)
_, dec_yrs_2019, timeseries_2019 = read_timeseries_from_nc(project_dir, var_name, year=2019)


plot_grids_with_differences(project_dir, var_name, depth,
                            dec_yrs_2008, timeseries_2008,
                            dec_yrs_2012, timeseries_2012,
                            dec_yrs_2017, timeseries_2017,
                            dec_yrs_2019, timeseries_2019)


