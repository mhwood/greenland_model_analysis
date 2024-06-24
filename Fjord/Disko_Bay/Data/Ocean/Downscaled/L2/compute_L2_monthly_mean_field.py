

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime


def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    drF = ds.variables['drF'][:]
    hFaC = ds.variables['HFacC'][:, :, :]
    rA = ds.variables['rA'][:, :]
    ds.close()
    Z_bottom = np.cumsum(drF)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2
    return(XC, YC, drF, Z, Depth, hFaC, rA)


def compute_mean_field(results_dir, var_name, year, month, hFaC, rA, drF, max_depth_index):

    if var_name=='Total_Chl':
        file_path = os.path.join(results_dir, 'Chl01', 'Chl01_' + str(year) + '{:02d}'.format(month) + '.nc')
        ds = nc4.Dataset(file_path)
        grid = ds.variables['Chl01'][:, :, :, :]
        # iterations = ds.variables['iterations'][:]
        ds.close()
        for chl_number in [2,3,4,5]:
            tmp_name = 'Chl'+'{:02d}'.format(chl_number)
            file_path = os.path.join(results_dir, tmp_name, tmp_name+'_' + str(year) + '{:02d}'.format(month) + '.nc')
            ds = nc4.Dataset(file_path)
            grid += ds.variables[tmp_name][:, :, :, :]
            ds.close()
    else:
        file_path = os.path.join(results_dir, var_name, var_name + '_' + str(year) + '{:02d}'.format(month) + '.nc')
        ds = nc4.Dataset(file_path)
        grid = ds.variables[var_name][:, :, :, :]
        # iterations = ds.variables['iterations'][:]
        ds.close()

    grid = np.mean(grid, axis=0)

    total_grid = np.zeros((np.shape(grid)[-2],np.shape(grid)[-1]))
    column_volume = np.zeros((np.shape(grid)[-2],np.shape(grid)[-1]))

    for d in range(max_depth_index+1):
        level_volume = rA*drF[d]
        column_volume += level_volume
        total_grid += grid[d,:,:] * level_volume

    total_grid = total_grid/column_volume

    return(total_grid)


def write_timeseries_to_nc(project_dir, year, var_name, XC, YC, annual_grid):

    output_file = os.path.join(project_dir,'Data','Modeling','L2',var_name+'_'+str(year)+'_monthly_mean_profile_fields.nc')
    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('longitude',np.shape(XC)[1])
    ds.createDimension('latitude',np.shape(XC)[0])
    ds.createDimension('month', 12)

    ivar = ds.createVariable('longitude','f4',('latitude','longitude'))
    ivar[:, :] = XC

    dvar = ds.createVariable('latitude', 'f4', ('latitude','longitude'))
    dvar[:, :] = YC

    mvar = ds.createVariable('months', 'f4', ('month', ))
    mvar[:] = np.arange(1,13)

    vvar = ds.createVariable(var_name, 'f4', ('month','latitude','longitude'))
    vvar[:,:,:] = annual_grid

    ds.close()


config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/' \
             'configurations/downscale_darwin'

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Disko Bay'

year = 2019
model_name = 'L2_Disko_Bay'
var_name = 'Total_Chl'
# var_name = 'Chl02'

max_depth = 100

results_dir = os.path.join(config_dir,'L2','L2_Disko_Bay','results_30_'+str(year))

XC, YC, drF, Z, Depth, hFaC, rA = read_grid_geometry_from_nc(config_dir, model_name)

max_depth_index = np.sum(np.cumsum(drF)<max_depth)

if year%4==0:
    timeseries = np.zeros((53, 365))
    iterations = np.zeros((365,))
else:
    timeseries = np.zeros((53, 364))
    iterations = np.zeros((364,))

for month in range(1,13):
    print('     - Reading in month '+str(month))
    grid = compute_mean_field(results_dir, var_name, year, month, hFaC, rA, drF, max_depth_index)

    if month==1:
        annual_grid = np.zeros((12,np.shape(grid)[-2], np.shape(grid)[-1]))

    annual_grid[month-1] = grid

    # C = plt.imshow(grid, origin='lower')
    # plt.colorbar(C)
    # plt.show()

write_timeseries_to_nc(project_dir, year, var_name, XC, YC, annual_grid)



