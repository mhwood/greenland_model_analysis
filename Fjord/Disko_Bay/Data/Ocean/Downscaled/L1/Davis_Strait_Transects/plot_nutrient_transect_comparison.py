
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon

def read_data_from_nc(folder, file_name):
    ds = nc4.Dataset(os.path.join(folder, 'Data', 'Nitrate', file_name))
    lon = ds.variables['longitude'][:]
    lat = ds.variables['latitude'][:]
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    years = ds.variables['year'][:]
    months = ds.variables['month'][:]
    depth = ds.variables['depth'][:]
    distance = ds.variables['distance'][:]
    bathymetry = ds.variables['bathymetry'][:]
    var_grid = ds.variables['Nitrate'][:, :, :]
    ds.close()
    return(lon, lat, x, y, years, months, depth, distance, bathymetry, var_grid)

def plot_annual_comparison(folder, depth, distance, bathymetry, insitu_var_grid, model_var_grid, years, months):

    vmax = np.max([np.max(model_var_grid[model_var_grid<20]), np.max(insitu_var_grid[insitu_var_grid<20])])

    for yy in range(np.shape(model_var_grid)[0]):
        fig = plt.figure(figsize=(12,5))

        gs1 = GridSpec(1, 3, left=0.05, right=0.95, bottom=0.05, top=0.95)

        ax1 = fig.add_subplot(gs1[0, 0])
        C = ax1.pcolormesh(distance, depth, insitu_var_grid[yy, :, :], vmin=0, vmax=vmax, cmap='turbo')
        ax1.set_ylim([np.max(bathymetry),0])
        cbar = plt.colorbar(C, orientation='horizontal')
        cbar.set_label('Nitrate Conc. ($\mu$M)')
        bathy = np.vstack([np.column_stack([distance,bathymetry]),
                          np.flipud(np.column_stack([distance,np.ones_like(bathymetry)*(np.max(bathymetry)+10)]))])
        poly = Polygon(bathy, facecolor='silver', edgecolor='k')
        ax1.add_patch(poly)
        ax1.set_title('In Situ Data')

        ax2 = fig.add_subplot(gs1[0, 1])
        C = ax2.pcolormesh(distance, depth, model_var_grid[yy, :, :], vmin=0, vmax=vmax, cmap='turbo')
        ax2.set_ylim([np.max(bathymetry), 0])
        cbar = plt.colorbar(C, orientation='horizontal')
        cbar.set_label('Nitrate Conc. ($\mu$M)')
        poly = Polygon(bathy, facecolor='silver', edgecolor='k')
        ax2.add_patch(poly)
        ax2.set_title('Model Data')

        ax3 = fig.add_subplot(gs1[0, 2])
        diff = model_var_grid[yy, :, :]-insitu_var_grid[yy, :, :]
        dvmax = np.max(np.abs(diff[diff<20]))
        dvmin = -1*dvmax
        C = ax3.pcolormesh(distance, depth, diff, vmin=dvmin, vmax=dvmax, cmap='seismic')
        ax3.set_ylim([np.max(bathymetry), 0])
        cbar = plt.colorbar(C, orientation='horizontal')
        cbar.set_label('Nitrate Conc. ($\mu$M)')
        poly = Polygon(bathy, facecolor='silver', edgecolor='k')
        ax3.add_patch(poly)
        ax3.set_title('Model - In Situ')

        date_str = str(years[yy])+'{:02d}'.format(months[yy])
        plt.savefig(os.path.join(folder,'Figures','Nutrients','Davis_Strait_Nitrate_Profile_'+date_str+'.png'))
        plt.close(fig)
    a=1

folder = '/Users/mike/Documents/Research/Projects/Disko Bay'
transect_number = 3
insitu_file_name = 'Transect_' + str(transect_number) + '_Nitrate_Profile.nc'
model_file_name = 'L1_Transect_' + str(transect_number) + '_Nitrate_Profile.nc'

lon, lat, x, y, years, months, depth, distance, bathymetry, insitu_var_grid = read_data_from_nc(folder, insitu_file_name)
years_subset = []
months_subset = []
for yy in range(len(years)):
    if years[yy] in [2010,2011]:
        years_subset.append(years[yy])
        months_subset.append(months[yy])
years = years_subset
months = months_subset

_, _, _, _, _, _, _, _, _, model_var_grid = read_data_from_nc(folder, model_file_name)

plot_annual_comparison(folder, depth, distance, bathymetry, insitu_var_grid, model_var_grid, years, months)





