
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.path as mplPath
import matplotlib.pyplot as plt
from pyproj import Transformer
from scipy.interpolate import interp1d
import shapefile

def read_sample_area_from_shapefile(project_dir,shapefile_name):
    shp_path = os.path.join(project_dir,'Map','Shapefiles','Ocean_Sample_Areas',shapefile_name)
    sf = shapefile.Reader(shp_path)

    shapes = sf.shapes()
    records = sf.records()

    sample_area = np.array(shapes[0].points)

    # points = []
    # for s in range(len(records)):
    #     if records[s][0]==sample_area_name:
    #         points = shapes[s].points
    #         points = np.array(points)

    return(sample_area)

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
            rows, cols = np.where(inside)
            for ri in range(len(rows)):
                sample_tiles.append(tile)
                sample_rows.append(rows[ri])
                sample_cols.append(cols[ri])

    return(sample_tiles, sample_rows, sample_cols)

def generate_ecco_field_timeseries_by_profiles(ecco_dir, field_name, sample_tiles, sample_rows, sample_cols,
                                   min_depth, max_depth, ecco_RC):

    years = np.arange(1992,1995)#2022)
    # years = np.arange(2017,2020)
    # years= [1992]

    timeseries = np.zeros((12*len(years),3))

    timestep_counter = 0

    interp_depth = np.arange(min_depth,max_depth)

    for year in years:
        print('   - Reading file in year '+str(year))
        file_path = os.path.join(ecco_dir,field_name,field_name+'_'+str(year)+'.nc')
        ds = nc4.Dataset(file_path)
        grid = ds.variables[field_name][:,:,:,:,:]
        for timestep in range(12):

            sum_profile = np.zeros((50,))
            count_profile = np.zeros((50,))

            for ri in range(len(sample_tiles)):
                profile = grid[timestep,:,sample_tiles[ri]-1,sample_rows[ri],sample_cols[ri]]

                sum_profile[profile!=0] += profile[profile!=0]
                count_profile[profile!=0] += 1

            sum_profile[count_profile!=0] = sum_profile[count_profile!=0]/count_profile[count_profile!=0]

            set_int = interp1d(np.abs(ecco_RC),sum_profile)
            interp_profile = set_int(interp_depth)

            timeseries[timestep_counter + timestep, 1] = np.mean(interp_profile)
            timeseries[timestep_counter + timestep, 2] = np.std(interp_profile)

        ds.close()

        first_index = 12*(year-years[0])
        last_index = 12*(year-years[0]+1)
        timeseries[first_index:last_index,0] = year+np.arange(1/24, 1+1/24, 1/12)

        timestep_counter+=12

    plt.plot(timeseries[:,0],timeseries[:,1],'g-')
    # plt.plot(timeseries[:, 0], timeseries[:, 1]+timeseries[:, 2],'-',alpha=0.5)
    # plt.plot(timeseries[:, 0], timeseries[:, 1] - timeseries[:, 2], '-', alpha=0.5)
    plt.show()

    return (timeseries)

def generate_ecco_field_timeseries(ecco_dir, field_name, sample_tiles, sample_rows, sample_cols,
                                   min_depth, max_depth, ecco_RC):

    years = np.arange(1992,2022)
    # years = np.arange(2017,2020)
    # years= [1992]

    timeseries = np.zeros((12*len(years),3))

    timestep_counter = 0

    interp_depth = np.arange(min_depth,max_depth)

    for year in years:
        print('   - Reading file in year '+str(year))
        file_path = os.path.join(ecco_dir,field_name,field_name+'_'+str(year)+'.nc')
        ds = nc4.Dataset(file_path)
        grid = ds.variables[field_name][:,:,:,:,:]
        for timestep in range(12):

            sum_profile = np.zeros((50,))
            count_profile = np.zeros((50,))

            # print('sample_tiles',sample_tiles)

            for ri in range(len(sample_tiles)):
                profile = grid[timestep,:,sample_tiles[ri]-1,sample_rows[ri],sample_cols[ri]]

                sum_profile[profile!=0] += profile[profile!=0]
                count_profile[profile!=0] += 1

            sum_profile[count_profile!=0] = sum_profile[count_profile!=0]/count_profile[count_profile!=0]

            set_int = interp1d(np.abs(ecco_RC),sum_profile)
            interp_profile = set_int(interp_depth)

            timeseries[timestep_counter + timestep, 1] = np.mean(interp_profile)
            timeseries[timestep_counter + timestep, 2] = np.std(interp_profile)

        ds.close()

        first_index = 12*(year-years[0])
        last_index = 12*(year-years[0]+1)
        timeseries[first_index:last_index,0] = year+np.arange(1/24, 1+1/24, 1/12)

        timestep_counter+=12

    # plt.plot(timeseries[:,0],timeseries[:,1],'g-')
    # # plt.plot(timeseries[:, 0], timeseries[:, 1]+timeseries[:, 2],'-',alpha=0.5)
    # # plt.plot(timeseries[:, 0], timeseries[:, 1] - timeseries[:, 2], '-', alpha=0.5)
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

ecco_dir = '/Volumes/helheim/Data_Repository/ECCO'
field_name = 'THETA'

# get sample area from shapefile
shapefile_name = 'IPCC Sample Area Off Shelf'

sample_area_name = 'IPCC'
min_depth = 200
max_depth = 500

sample_polygon = read_sample_area_from_shapefile(project_dir,shapefile_name)
# sample_polygon = reproject_points(sample_polygon, inputCRS=3413, outputCRS=4326)

ecco_XC_tiles, ecco_YC_tiles, ecco_RC = read_ecco_geometry(ecco_dir)

sample_tiles, sample_rows, sample_cols = find_sample_points_in_ecco_tiles(sample_polygon,ecco_XC_tiles, ecco_YC_tiles)

timeseries = generate_ecco_field_timeseries(ecco_dir, field_name, sample_tiles, sample_rows, sample_cols,
                               min_depth, max_depth, ecco_RC)

output_file = os.path.join(project_dir,'Data','Modeling','ECCO','ECCOv5_'+field_name+'_'+sample_area_name+'_timeseries.nc')
write_timeseries_to_nc(output_file,timeseries,field_name)



