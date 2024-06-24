

import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, interp1d


def read_insitu_data_from_nc(folder, transect_number):
    ds = nc4.Dataset(os.path.join(folder, 'Data', 'Nitrate', 'Transect_' + str(transect_number) + '_Nitrate_Profile.nc'))
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

def interpolate_model_profiles_to_nc(results_dir, lon, lat, years, months, depth, distance, bathymetry, var_grid):

    output_var_grid = np.zeros((len(years),len(depth), len(distance)))

    for yy in range(len(years)):
        datestr = str(years[yy])+'{:02d}'.format(months[yy])
        file_to_use = ''
        for file_name in os.listdir(results_dir):
            if file_name.split('.')[1]==datestr:
                file_to_use=file_name

        print('        - Reading from file '+str(file_to_use))
        ds=nc4.Dataset(os.path.join(results_dir,file_to_use))
        model_lon = np.array(ds.variables['longitude'][:, :])
        model_lat = np.array(ds.variables['latitude'][:, :])
        model_var_grid = ds.variables['NO3'][:, :, :, :]
        model_depth = np.array(ds.variables['depths'][:])
        ds.close()
        model_var_grid = model_var_grid[0,:,:,:]

        min_dist = (model_lon - np.min(lon))**2 + (model_lat - np.min(lat))**2
        min_row, min_col = np.where(min_dist==np.min(min_dist))
        min_row = min_row[0]
        min_col = min_col[0]

        max_dist = (model_lon - np.max(lon)) ** 2 + (model_lat - np.max(lat)) ** 2
        max_row, max_col = np.where(max_dist == np.min(max_dist))
        max_row = max_row[0]
        max_col = max_col[0]

        model_lon = model_lon[min_row:max_row+1, min_col:max_col]
        model_lat = model_lat[min_row:max_row + 1, min_col:max_col]
        model_var_grid = model_var_grid[:, min_row:max_row + 1, min_col:max_col]

        points = np.column_stack([model_lon.ravel(), model_lat.ravel()])

        # interpolate onto the transect first
        print('        - Interpolating field to transect')
        intermediate_var_grid = np.zeros((len(model_depth), len(distance)))
        for d in range(len(model_depth)):
            indices = model_var_grid[d,:,:].ravel()!=0
            if np.sum(indices)>0:
                points_subset = points[indices,:]
                values = model_var_grid[d, :, :].ravel()[indices]
                interp_layer = griddata(points_subset, values, (lon, lat), fill_value=0)
                intermediate_var_grid[d,:] = interp_layer

        # plt.subplot(2,1,1)
        # C = plt.pcolormesh(distance, model_depth, intermediate_var_grid)
        # plt.colorbar(C)
        # plt.gca().set_ylim([1100,0])
        # plt.subplot(2,1, 2)
        # C = plt.pcolormesh(distance, depth, var_grid[yy, :, :])
        # plt.colorbar(C)
        # plt.gca().set_ylim([1100, 0])
        # plt.show()

        # then interpolate vertically
        var_grid_interpolated = np.zeros((len(depth), len(distance)))
        for d in range(len(distance)):
            non_zero_indices = intermediate_var_grid[:,d]!=0
            if np.count_nonzero(non_zero_indices)>1:
                k_max = np.sum(non_zero_indices)
                set_int = interp1d(model_depth[non_zero_indices], intermediate_var_grid[non_zero_indices, d],
                                   bounds_error=False,
                                   fill_value=(intermediate_var_grid[0,d],  intermediate_var_grid[k_max-1, d]))
                var_grid_interpolated[:,d] = set_int(np.array(depth))

                # mask out the dry depths
                depth_indices = depth >= bathymetry[d]+1
                var_grid_interpolated[depth_indices, d]=0

        # plt.subplot(2,1,1)
        # C = plt.pcolormesh(distance, depth, var_grid_interpolated)
        # plt.colorbar(C)
        # plt.gca().set_ylim([1100,0])
        # plt.subplot(2,1, 2)
        # C = plt.pcolormesh(distance, depth, var_grid[yy, :, :])
        # plt.colorbar(C)
        # plt.gca().set_ylim([1100, 0])
        # plt.show()

        output_var_grid[yy, :, :] = var_grid_interpolated
    return(output_var_grid)

def outputting_transects_to_nc(folder, transect_number, years, months,
                               depth, distance, bathymetry, x, y, lon, lat, interpolated_var_grid):

    ds = nc4.Dataset(os.path.join(folder,'Data','Nitrate','L1_Transect_'+str(transect_number)+'_Nitrate_Profile.nc'),'w')

    ds.createDimension('time',len(years))
    ds.createDimension('depth',len(depth))
    ds.createDimension('distance', len(distance))

    var = ds.createVariable('year','i4',('time',))
    var[:] = years

    var = ds.createVariable('month', 'i4', ('time',))
    var[:] = months

    var = ds.createVariable('depth', 'f4', ('depth',))
    var[:] = depth

    var = ds.createVariable('distance', 'f4', ('distance',))
    var[:] = distance

    var = ds.createVariable('bathymetry', 'f4', ('distance',))
    var[:] = bathymetry

    var = ds.createVariable('x', 'f4', ('distance',))
    var[:] = x

    var = ds.createVariable('y', 'f4', ('distance',))
    var[:] = y

    var = ds.createVariable('longitude', 'f4', ('distance',))
    var[:] = lon

    var = ds.createVariable('latitude', 'f4', ('distance',))
    var[:] = lat

    var = ds.createVariable('Nitrate', 'f4', ('time','depth','distance'))
    var[:,:,:] = interpolated_var_grid

    ds.close()



folder = '/Users/mhwood/Documents/Research/Projects/Disko Bay'
config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/' \
             'configurations/downscale_darwin'
model_name = 'L1_W_Greenland'
results_dir = os.path.join(config_dir,'L1',model_name,'results')
transect_number = 3

# step 1: read in the nutrient profiles
print('  - Reading in the in situ data')
lon, lat, x, y, years, months, depth, distance, bathymetry, var_grid = read_insitu_data_from_nc(folder, transect_number)
years_subset = []
months_subset = []
for yy in range(len(years)):
    if years[yy] in [2010,2011]:
        years_subset.append(years[yy])
        months_subset.append(months[yy])
years = years_subset
months = months_subset
print(years, months)

# step 2: read in the model profiles on the grid
print('  - Interpolating the model files')
interpolated_var_grid = interpolate_model_profiles_to_nc(results_dir, lon, lat, years, months, depth, distance, bathymetry, var_grid)

# step 3: store the model profiles to nc
outputting_transects_to_nc(folder, transect_number, years, months,
                           depth, distance, bathymetry, x, y, lon, lat, interpolated_var_grid)
