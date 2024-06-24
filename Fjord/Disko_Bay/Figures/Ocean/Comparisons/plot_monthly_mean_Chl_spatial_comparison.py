


import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec


def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    ds.close()
    return(XC, YC, Depth)


def read_grid_from_nc(project_dir, var_name, year, month):

    output_file = os.path.join(project_dir,'Data','Modeling','L2',var_name+'_'+str(year)+'_monthly_mean_profile_fields.nc')
    ds = nc4.Dataset(output_file)

    grid = ds.variables[var_name][:, :, :]
    grid = grid[month-1, : , :]

    ds.close()

    return(grid)


def plot_grids_with_differences(project_dir, var_name, XC, YC, grid_2008, grid_2012, grid_2017):

    fig = plt.figure(figsize=(10, 10))

    image_rows = 6
    image_cols = 6
    gs1 = GridSpec(image_rows * 4 + 3, image_cols * 3 + 2,
                   left=0.1, right=0.95, bottom=0.05, top=0.95)

    ########################################################################
    # 2008

    plot_row = 0
    ax21m = fig.add_subplot(gs1[plot_row * image_rows + 1:(plot_row + 1) * image_rows + 1, :image_cols])
    # ax21m.imshow(np.flipud(landsat_image),
    #              extent=[np.min(landsat_X) / 1000, np.max(landsat_X) / 1000, np.min(landsat_Y) / 1000,
    #                      np.max(landsat_Y) / 1000],
    #              alpha=0.8)
    ax21m.pcolormesh(grid_2008, cmap='turbo', vmin=0, vmax=5)
    # ax21m.set_xlim([np.min(polygon[:, 0] / 1000), np.max(polygon[:, 0] / 1000)])
    # ax21m.set_ylim([np.min(polygon[:, 1] / 1000), np.max(polygon[:, 1] / 1000)])
    ax21m.set_xticks([])
    ax21m.set_yticks([])
    # ax21m.text(np.max(polygon[:, 0] / 1000) - 1, np.max(polygon[:, 1] / 1000) - 1, '2016/03/23', ha='right', va='top',
    #            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    ########################################################################
    # 2012

    plot_row = 1
    ax21m = fig.add_subplot(gs1[plot_row * image_rows + 1:(plot_row + 1) * image_rows + 1, :image_cols])
    # ax21m.imshow(np.flipud(landsat_image),
    #              extent=[np.min(landsat_X) / 1000, np.max(landsat_X) / 1000, np.min(landsat_Y) / 1000,
    #                      np.max(landsat_Y) / 1000],
    #              alpha=0.8)
    ax21m.pcolormesh(grid_2012, cmap='turbo', vmin=0, vmax=5)
    # ax21m.set_xlim([np.min(polygon[:, 0] / 1000), np.max(polygon[:, 0] / 1000)])
    # ax21m.set_ylim([np.min(polygon[:, 1] / 1000), np.max(polygon[:, 1] / 1000)])
    ax21m.set_xticks([])
    ax21m.set_yticks([])
    # ax21m.text(np.max(polygon[:, 0] / 1000) - 1, np.max(polygon[:, 1] / 1000) - 1, '2016/03/23', ha='right', va='top',
    #            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    ########################################################################
    # 2017

    plot_row = 2
    ax21m = fig.add_subplot(gs1[plot_row * image_rows + 1:(plot_row+1) * image_rows + 1, :image_cols])
    # ax21m.imshow(np.flipud(landsat_image),
    #              extent=[np.min(landsat_X) / 1000, np.max(landsat_X) / 1000, np.min(landsat_Y) / 1000,
    #                      np.max(landsat_Y) / 1000],
    #              alpha=0.8)
    ax21m.pcolormesh(grid_2017, cmap='turbo', vmin=0, vmax=5)
    # ax21m.set_xlim([np.min(polygon[:, 0] / 1000), np.max(polygon[:, 0] / 1000)])
    # ax21m.set_ylim([np.min(polygon[:, 1] / 1000), np.max(polygon[:, 1] / 1000)])
    ax21m.set_xticks([])
    ax21m.set_yticks([])
    # ax21m.text(np.max(polygon[:, 0] / 1000) - 1, np.max(polygon[:, 1] / 1000) - 1, '2016/03/23', ha='right', va='top',
    #            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    output_file = os.path.join(project_dir,'Figures','Ocean','Modeling','L2','L2_Disko_Bay_'+var_name+'_Spatial_Comparison.png')
    plt.savefig(output_file)
    plt.close()

config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/' \
             'configurations/downscale_darwin'

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Disko Bay'

model_name = 'L2_Disko_Bay'
var_name = 'Total_Chl'
month = 8

XC, YC,  Depth = read_grid_geometry_from_nc(config_dir, model_name)

grid_2008 = read_grid_from_nc(project_dir, var_name, year=2008, month=month)
grid_2012 = read_grid_from_nc(project_dir, var_name, year=2012, month=month)
grid_2017 = read_grid_from_nc(project_dir, var_name, year=2017, month=month)

plot_grids_with_differences(project_dir, var_name, XC, YC, grid_2008, grid_2012, grid_2017)


