
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.path as mplPath
import matplotlib.pyplot as plt
from pyproj import Transformer
from scipy.interpolate import interp1d, griddata
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

def read_interpolated_ECCO_grid(project_folder, ecco_dir, ecco_XC_tiles, ecco_YC_tiles, year):

    output_file_name = 'L0_CE_Greenland_THETA_'+str(year)+'.nc'
    if output_file_name not in os.listdir(os.path.join(project_folder,'Data','Modeling','ECCO','interpolated_THETA')):

        print('    - '+output_file_name+' not yet created - creating it now')
        full_XC = np.concatenate([ecco_XC_tiles[3], np.rot90(ecco_XC_tiles[7], k=3)], axis=-2)
        full_YC = np.concatenate([ecco_YC_tiles[3], np.rot90(ecco_YC_tiles[7], k=3)], axis=-2)

        print('        - Sample polygon lon range: ' + str(np.min(sample_polygon[:, 0])) + ' to ' + str(
            np.max(sample_polygon[:, 0])))
        print('        - Sample polygon lat range: ' + str(np.min(sample_polygon[:, 1])) + ' to ' + str(
            np.max(sample_polygon[:, 1])))

        min_lon = np.min(np.floor(sample_polygon[:, 0]))
        max_lon = np.max(np.ceil(sample_polygon[:, 0]))
        min_lat = np.min(np.floor(sample_polygon[:, 1]))
        max_lat = np.max(np.ceil(sample_polygon[:, 1]))

        print('        - Interp grid lon range: ' + str(min_lon) + ' to ' + str(max_lon))
        print('        - Interp grid lat range: ' + str(min_lat) + ' to ' + str(max_lat))

        # plt.subplot(1, 2, 1)
        # C = plt.imshow(full_XC)
        # plt.colorbar(C)
        # plt.subplot(1,2,2)
        # C = plt.imshow(full_YC)
        # plt.colorbar(C)
        # plt.show()

        lon = np.arange(min_lon, max_lon + 1)
        lat = np.arange(min_lat, max_lat + 1)
        Lon, Lat = np.meshgrid(lon, lat)

        file_path = os.path.join(ecco_dir, field_name, field_name + '_' + str(year) + '.nc')
        ds = nc4.Dataset(file_path)
        grid = ds.variables[field_name][:, :, :, :, :]

        full_grid = np.concatenate([grid[:, :, 2, :, :],
                                    np.rot90(grid[:, :, 6, :, :], axes=(2, 3), k=3)], axis=-2)

        # we dont need the bottom?
        full_grid = full_grid[:,:29,:,:]

        # plt.imshow(full_grid[0,0,:,:],origin='lower')
        # plt.show()

        interp_grid = np.zeros((np.shape(full_grid)[0], np.shape(full_grid)[1],
                                np.shape(Lon)[0], np.shape(Lon)[1]))

        print('        - Interpolating the ECCO grid onto a regular grid')
        for t in range(np.shape(full_grid)[0]):
            print('               - Working on timestep '+str(t))
            for d in range(np.shape(full_grid)[1]):
                if d%10==0:
                    print('                    - Working on depth index '+str(d))
                interp_grid_slice = griddata(np.column_stack([full_XC.ravel(), full_YC.ravel()]),
                                             full_grid[t, d, :, :].ravel(), (Lon, Lat))

                # C = plt.pcolormesh(Lon, Lat, interp_grid_slice)
                # plt.plot(sample_polygon[:,0],sample_polygon[:,1],'k-')
                # plt.colorbar(C)
                # plt.show()

                interp_grid[t, d, :, :] = interp_grid_slice

        ds.close()

        ds = nc4.Dataset(os.path.join(project_folder,'Data','Modeling','ECCO','interpolated_THETA',output_file_name),'w')
        ds.createDimension('time',np.shape(interp_grid)[0])
        ds.createDimension('depth',np.shape(interp_grid)[1])
        ds.createDimension('latitude', np.shape(interp_grid)[2])
        ds.createDimension('longitude', np.shape(interp_grid)[3])

        lonvar = ds.createVariable('Longitude','f4',('latitude','longitude'))
        lonvar[:, :] = Lon

        latvar = ds.createVariable('Latitude', 'f4', ('latitude', 'longitude'))
        latvar[:, :] = Lat

        thetavar = ds.createVariable('THETA', 'f4', ('time', 'depth', 'latitude', 'longitude'))
        thetavar[:, :, :, :] = interp_grid

        ds.close()

    else:
        ds = nc4.Dataset(os.path.join(project_folder,'Data','Modeling','ECCO','interpolated_THETA',output_file_name))

        Lon = ds.variables['Longitude'][:, :]
        Lat = ds.variables['Latitude'][:, :]
        interp_grid = ds.variables['THETA'][:, :, :, :]

        ds.close()

    return(Lon, Lat, interp_grid)

def generate_ecco_field_timeseries(project_folder, ecco_dir, field_name, min_depth, max_depth, sample_polygon, ecco_XC_tiles, ecco_YC_tiles, ecco_RC):

    years = np.arange(1992,1993)#2022)
    # years = np.arange(2017,2020)
    # years= [1992]

    timeseries = np.zeros((12*len(years),3))

    timestep_counter = 0

    interp_depth = np.arange(min_depth,max_depth)
    p = mplPath.Path(sample_polygon)

    for year in years:
        print('   - Reading file in year '+str(year))

        Lon, Lat, interp_grid = read_interpolated_ECCO_grid(project_folder,ecco_dir, ecco_XC_tiles, ecco_YC_tiles, year)

        inside = p.contains_points(np.column_stack([Lon.ravel(), Lat.ravel()]))
        inside = np.reshape(inside, np.shape(Lon))

        for t in range(np.shape(interp_grid)[0]):

            mean_profile = np.zeros((len(ecco_RC), ))
            std_profile = np.zeros((len(ecco_RC), ))

            for d in range(np.shape(interp_grid)[1]):
                # if depth_indices[d]:
                depth_slice = interp_grid[t, d, :, :]
                inside_points = depth_slice[inside]
                if np.any(inside_points !=0):
                    inside_points = inside_points[inside_points !=0]
                    mean_profile[d] = np.mean(inside_points)
                    std_profile[d] = np.std(inside_points)

            # plt.plot(mean_profile,ecco_RC)
            # plt.gca().set_ylim([-900,0])
            # plt.show()

            print(np.shape(np.abs(ecco_RC[:np.shape(interp_grid)[1]])))
            print(np.shape(mean_profile[:np.shape(interp_grid)[1]]))

            set_int = interp1d(np.abs(ecco_RC[:np.shape(interp_grid)[1]]), mean_profile[:np.shape(interp_grid)[1]])
            interp_profile = set_int(interp_depth)

            timeseries[timestep_counter + t, 1] = np.mean(interp_profile)
            timeseries[timestep_counter + t, 2] = np.std(interp_profile)

        first_index = 12*(year-years[0])
        last_index = 12*(year-years[0]+1)
        timeseries[first_index:last_index,0] = year+np.arange(1/24, 1+1/24, 1/12)

        timestep_counter+=12

    plt.plot(timeseries[:,0],timeseries[:,1],'g-')
    # plt.plot(timeseries[:, 0], timeseries[:, 1]+timeseries[:, 2],'-',alpha=0.5)
    # plt.plot(timeseries[:, 0], timeseries[:, 1] - timeseries[:, 2], '-', alpha=0.5)
    plt.show()

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

# get sample area from shapefile
shapefile_name = 'IPCC Sample Area Off Shelf'

sample_area_name = 'IPCC'
min_depth = 200
max_depth = 500

is_extension = True

sample_polygon = read_sample_area_from_shapefile(project_dir,shapefile_name)
# sample_polygon = reproject_points(sample_polygon, inputCRS=3413, outputCRS=4326)

ecco_XC_tiles, ecco_YC_tiles, ecco_RC = read_ecco_geometry(ecco_dir)

print(len(ecco_RC[ecco_RC>-1000]))

# sample_tiles, sample_rows, sample_cols = find_sample_points_in_ecco_tiles(sample_polygon,ecco_XC_tiles, ecco_YC_tiles)

timeseries = generate_ecco_field_timeseries(project_dir, ecco_dir, field_name, min_depth, max_depth, sample_polygon, ecco_XC_tiles, ecco_YC_tiles, ecco_RC)


# output_file = os.path.join(project_dir,'Data','Modeling','ECCO','ECCOv5_'+field_name+'_'+sample_area_name+'_timeseries.nc')
# write_timeseries_to_nc(output_file,timeseries,field_name)



