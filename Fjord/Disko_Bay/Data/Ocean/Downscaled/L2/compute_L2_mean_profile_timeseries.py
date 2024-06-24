

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
    ds.close()
    Z_bottom = np.cumsum(drF)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2
    return(XC, YC, Z, Depth, hFaC)


def compute_mean_profile_timeseries(results_dir, var_name, year, month, hFaC):

    if var_name=='Total_Chl':
        file_path = os.path.join(results_dir, 'Chl01', 'Chl01_' + str(year) + '{:02d}'.format(month) + '.nc')
        ds = nc4.Dataset(file_path)
        grid = ds.variables['Chl01'][:, :, :, :]
        iterations = ds.variables['iterations'][:]
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
        iterations = ds.variables['iterations'][:]
        ds.close()

    timeseries = np.zeros((np.shape(grid)[1],np.shape(grid)[0]))

    for t in range(np.shape(grid)[0]):
        for d in range(np.shape(grid)[1]):
            level_set = grid[t,d,:,:]
            if np.any(hFaC[d,:,:]!=0):
                timeseries[d,t] = np.mean(level_set[hFaC[d,:,:]!=0])

    return(iterations,timeseries)


def write_timeseries_to_nc(project_dir, year, var_name, iterations, depth, timeseries):

    output_file = os.path.join(project_dir,'Data','Modeling','L2',var_name+'_'+str(year)+'_mean_profile_timeseries.nc')
    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('depth',len(depth))
    ds.createDimension('iterations',len(iterations))

    ivar = ds.createVariable('iterations','f4',('iterations',))
    ivar[:] = iterations

    dvar = ds.createVariable('depth', 'f4', ('depth',))
    dvar[:] = depth

    vvar = ds.createVariable(var_name, 'f4', ('depth','iterations'))
    vvar[:,:] = timeseries

    ds.close()


config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/' \
             'configurations/downscale_darwin'

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Disko Bay'

year = 2019
model_name = 'L2_Disko_Bay'
var_name = 'Total_Chl'

results_dir = os.path.join(config_dir,'L2','L2_Disko_Bay','results_30_'+str(year))

XC, YC, Z, Depth, hFaC = read_grid_geometry_from_nc(config_dir, model_name)

if year%4==0:
    timeseries = np.zeros((53, 366))
    iterations = np.zeros((366,))
else:
    timeseries = np.zeros((53, 365))
    iterations = np.zeros((365,))

counter = 0
for month in range(1,13):
    print('     - Reading in month '+str(month))
    month_iterations, month_timeseries = compute_mean_profile_timeseries(results_dir, var_name, year, month, hFaC)
    # plt.imshow(month_timeseries)
    # plt.show()

    timeseries[:,counter:counter+np.shape(month_timeseries)[1]] = month_timeseries
    iterations[counter:counter+np.shape(month_timeseries)[1]] = month_iterations

    counter += np.shape(month_timeseries)[1]

C = plt.pcolormesh(iterations, Z, timeseries, shading='nearest')
plt.colorbar(C)
plt.show()

write_timeseries_to_nc(project_dir, year, var_name, iterations, Z, timeseries)





