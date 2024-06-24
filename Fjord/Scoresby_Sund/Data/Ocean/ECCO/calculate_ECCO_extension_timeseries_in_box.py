
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
import matplotlib.path as mplPath
from pyproj import Transformer
from scipy.interpolate import interp1d
import cmocean.cm as cm
import shutil
import argparse
import sys
import shapefile

def read_sample_area_from_shapefile(project_dir,sample_area_name):
    if sample_area_name=='Scoresby_Shelf':
        shp_path = os.path.join(project_dir,'Map','Shapefiles','Ocean_Sample_Areas','Scoresby_Shelf_Ocean_Sample_Area')
    else:
        shp_path = os.path.join(project_dir, 'Map', 'Shapefiles', 'Ocean_Sample_Areas', sample_area_name)
    sf = shapefile.Reader(shp_path)

    shapes = sf.shapes()
    records = sf.records()
    # points = []
    # for s in range(len(records)):
    #     if records[s][0]==sample_area_name:
    points = shapes[0].points
    points = np.array(points)

    return(points)

def reproject_points(points,inputCRS,outputCRS,x_column=0,y_column=1):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(points[:, x_column], points[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
        run_test = False
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon = np.copy(points)
    output_polygon[:, x_column] = x2
    output_polygon[:, y_column] = y2
    return output_polygon

def read_ecco_geometry(ecco_dir):

    XC_tiles = {}
    YC_tiles = {}

    for tile in [3,4,7,11]:#range(1,14):
        file_path = os.path.join(ecco_dir, 'GRID', 'GRID.' + '{:04d}'.format(tile) + '.nc')
        # print('Reading from grid GRID.' + '{:04d}'.format(tile+1) + '.nc')
        ds = nc4.Dataset(file_path)
        XC = ds.variables['XC'][:, :]
        YC = ds.variables['YC'][:, :]
        RC = ds.variables['RC'][:]
        ds.close()

        XC_tiles[tile] = XC
        YC_tiles[tile] = YC

    return(XC_tiles, YC_tiles, RC)

def find_sample_points_in_ecco_tiles(sample_polygon, ecco_XC_tiles, ecco_YC_tiles):

    p = mplPath.Path(sample_polygon)

    sample_tiles = []
    sample_rows = []
    sample_cols = []
    for tile in list(ecco_XC_tiles.keys()):
        tile_points = np.column_stack([ecco_XC_tiles[tile].ravel(), ecco_YC_tiles[tile].ravel()])
        inside = p.contains_points(tile_points)
        inside = np.reshape(inside,np.shape(ecco_XC_tiles[tile]))

        if np.any(inside==1):
            # plt.plot(sample_polygon[:,0],sample_polygon[:,1], 'g-')
            # plt.plot(tile_points[inside,0], tile_points[inside,1], 'k.')
            # plt.title('Tile: '+str(tile))
            # plt.show()
            rows, cols = np.where(inside)
            for ri in range(len(rows)):
                sample_tiles.append(tile)
                sample_rows.append(rows[ri])
                sample_cols.append(cols[ri])

    return(sample_tiles, sample_rows, sample_cols)

def generate_ecco_field_timeseries(ecco_dir, field_name, sample_tiles, sample_rows, sample_cols,
                                   min_depth, max_depth, ecco_RC, is_extension):

    if not is_extension:
        years = np.arange(1992, 2018)
    else:
        years = np.arange(2017, 2022)

    if is_extension:
        timeseries = np.zeros((12*(len(years)-1)+6,3))
    else:
        timeseries = np.zeros((12 * len(years), 3))

    timestep_counter = 0

    interp_depth = np.arange(min_depth,max_depth)

    for year in years:
        print('   - Reading file in year '+str(year))
        if is_extension:
            file_path = os.path.join(ecco_dir,field_name+'_Extension',field_name+'_'+str(year)+'.nc')
        else:
            file_path = os.path.join(ecco_dir, field_name, field_name + '_' + str(year) + '.nc')
        ds = nc4.Dataset(file_path)
        grid = ds.variables[field_name][:,:,:,:,:]
        for timestep in range(np.shape(grid)[0]):
            # print(timestep, timestep_counter + timestep)

            all_profiles = np.zeros((len(interp_depth),len(sample_tiles)))

            for ri in range(len(sample_tiles)):
                profile = grid[timestep,:,sample_tiles[ri]-1,sample_rows[ri],sample_cols[ri]]

                if np.any(profile[np.abs(ecco_RC)<max_depth]==0):
                    # print('Note: Skipping some profiles which dont cover the entire requested depth range')
                    interp_profile = np.zeros_like(interp_depth)
                    # sub_ecco_RC = ecco_RC[profile!=0]
                    #
                    # if np.any(np.abs(sub_ecco_RC)>min_depth):
                    #     print(ri)
                    #     sub_profile = profile[profile!=0]
                    #     set_int = interp1d(np.abs(sub_ecco_RC), sub_profile)
                    #
                    #     interp_depth_sub = interp_depth[interp_depth<np.max(np.abs(sub_ecco_RC))]
                    #     interp_profile_sub = set_int(interp_depth_sub)
                    #     interp_profile[:len(interp_profile_sub)]=interp_profile_sub
                    #
                    #     # plt.plot(sub_profile, np.abs(sub_ecco_RC),'k.')
                    #     # plt.plot(interp_profile_sub, interp_depth_sub)
                    #     # plt.show()
                else:
                    set_int = interp1d(np.abs(ecco_RC),profile)
                    interp_profile = set_int(interp_depth)

                # plt.plot(profile, np.abs(ecco_RC), '-.',color='silver')
                # plt.plot(interp_profile, interp_depth)
                # plt.show()

                all_profiles[:,ri] = interp_profile

            # X, D = np.meshgrid(np.arange(len(sample_tiles)),interp_depth)
            # C = plt.pcolormesh(X,D,all_profiles,cmap=cm.thermal)
            # plt.colorbar(C)
            # plt.show()

            timeseries[timestep_counter + timestep, 1] = np.mean(all_profiles[all_profiles!=0])
            timeseries[timestep_counter + timestep, 2] = np.std(all_profiles[all_profiles!=0])

        ds.close()


        if is_extension and year == 2017:
            year_time = year+np.arange(1/2+1/12+1/24, 1+1/24, 1/12)
            first_index = 0
            last_index = 6
        else:
            year_time = year + np.arange(1 / 24, 1 + 1 / 24, 1 / 12)
            first_index = 12 * (year - years[0])
            last_index = 12 * (year - years[0] + 1)
            if is_extension:
                first_index -= 6
                last_index -= 6

        # print(first_index, last_index)
        # year_time = year_time[:np.shape(grid)[0]]
        timeseries[first_index:last_index,0] = year_time

        if is_extension and year==2017:
            timestep_counter += 6
        else:
            timestep_counter += 12

    # plt.plot(timeseries[:,0],timeseries[:,1],'g-')
    # plt.plot(timeseries[:, 0], timeseries[:, 1]+timeseries[:, 2],'-',alpha=0.5)
    # plt.plot(timeseries[:, 0], timeseries[:, 1] - timeseries[:, 2], '-', alpha=0.5)
    # plt.show()

    return (timeseries)

def write_timeseries_to_nc(output_file,timeseries,field_name):
    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('time',np.shape(timeseries)[0])

    tvar = ds.createVariable('dec_yr','f4',('time',))
    tvar[:] = timeseries[:,0]
    mvar = ds.createVariable(field_name+'_mean', 'f4', ('time',))
    mvar[:] = timeseries[:, 1]
    svar = ds.createVariable(field_name+'_std', 'f4', ('time',))
    svar[:] = timeseries[:, 2]

    ds.close()

project_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund'

ecco_dir = '/Volumes/petermann/Research/Projects/Ocean_Modeling/ECCO'
field_name = 'THETA'

sample_area_name = 'Scoresby_Shelf'
# sample_area_name = 'IPCC Sample Area Off Shelf'

if sample_area_name=='IPCC Sample Area Off Shelf':
    min_depth = 200
    max_depth = 500
if sample_area_name == 'Scoresby_Shelf':
    min_depth = 200
    max_depth = 450


is_extension = True

sample_polygon = read_sample_area_from_shapefile(project_dir,sample_area_name)
# sample_polygon = reproject_points(sample_polygon, inputCRS=3413, outputCRS=4326)

ecco_XC_tiles, ecco_YC_tiles, ecco_RC = read_ecco_geometry(ecco_dir)

sample_tiles, sample_rows, sample_cols = find_sample_points_in_ecco_tiles(sample_polygon,ecco_XC_tiles, ecco_YC_tiles)
print(len(sample_tiles))
print('    - First 5 sample points')
for i in range(5):
    print('         - tile '+str(sample_tiles[i])+', row '+str(sample_rows[i])+', col '+str(sample_cols[i]))

timeseries = generate_ecco_field_timeseries(ecco_dir, field_name, sample_tiles, sample_rows, sample_cols,
                               min_depth, max_depth, ecco_RC, is_extension)

if is_extension:
    timeseries = timeseries[:-1,:]
plt.plot(timeseries[:,0], timeseries[:,1])
plt.show()

sample_area_name_output = sample_area_name.split()[0]
if is_extension:
    output_file = os.path.join(project_dir,'Data','Modeling','ECCO','ECCOv5_'+field_name+'_'+sample_area_name_output+'_extension_timeseries.nc')
else:
    output_file = os.path.join(project_dir, 'Data', 'Modeling', 'ECCO', 'ECCOv5_' + field_name + '_' + sample_area_name_output + '_Alpha_timeseries.nc')
write_timeseries_to_nc(output_file,timeseries,field_name)



