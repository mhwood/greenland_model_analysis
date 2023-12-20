
import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt



plot_titles = ['Near-Glacier','Mid-Fjord','Fjord Entrance','Shelf Break']
output_names =['near_glacier','mid_fjord','fjord_entrance','shelf_break']

locations = [(-28.3818,71.9304),(-24.9112,71.0983),(-21.3842,70.1023),(-17.515,71.976)]

location_name = 'Sound_Entrance'
location = [(-21.379,70.091)]

year_months = [[2017,10],[2016,9]]


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

def great_circle_distance(lon_ref, lat_ref, Lon, Lat):
    earth_radius = 6371000
    lon_ref_radians = np.radians(lon_ref)
    lat_ref_radians = np.radians(lat_ref)
    lons_radians = np.radians(Lon)
    lats_radians = np.radians(Lat)
    lat_diff = lats_radians - lat_ref_radians
    lon_diff = lons_radians - lon_ref_radians
    d = np.sin(lat_diff * 0.5) ** 2 + np.cos(lat_ref_radians) * np.cos(lats_radians) * np.sin(lon_diff * 0.5) ** 2
    h = 2 * earth_radius * np.arcsin(np.sqrt(d))
    return(h)

def find_sample_point_in_ecco_tiles(lon ,lat, ecco_XC_tiles, ecco_YC_tiles):

    min_dist = 1e22
    for tile in list(ecco_XC_tiles.keys()):
        tile_XC = ecco_XC_tiles[tile]
        tile_YC = ecco_YC_tiles[tile]

        dist = great_circle_distance(lon, lat, tile_XC.ravel(), tile_YC.ravel())
        if np.min(dist)<min_dist:
            min_dist = np.min(dist)
            dist = np.reshape(dist,(np.shape(tile_XC)))
            row, col = np.where(dist==np.min(dist))
            row = row[0]
            col = col[0]
            tile_to_use = tile

    x = ecco_XC_tiles[tile_to_use][row, col]
    y = ecco_YC_tiles[tile_to_use][row, col]

    return(row,col,tile_to_use,min_dist,x,y)

def read_profiles_from_ECCO_extension_files(ecco_dir,field_name, row, col, tile, RC):

    years = np.arange(2019,2022)
    # years= [2019]

    timestep_counter = 0

    time = np.zeros((12*len(years),))
    timeseries = np.zeros((50, 12 * len(years)))

    # time = np.zeros((2,))
    # timeseries = np.zeros((50, 2))

    for year in years:
        print('   - Reading file in year ' + str(year))
        file_path = os.path.join(ecco_dir, field_name, field_name + '_' + str(year) + '.nc')
        ds = nc4.Dataset(file_path)
        grid = ds.variables[field_name][:, :, :, :, :]
        grid = grid[:,:,tile-1,row,col]
        ds.close()
        year_time = year+np.arange(1/24, 1+1/24, 1/12)
        time[timestep_counter:timestep_counter+12] = year_time
        for t in range(12):
            timeseries[:, timestep_counter+t] = grid[t,:]

        timestep_counter += 12

    # plt.plot(timeseries[:,0],timeseries[:,1],'g-')
    # plt.plot(timeseries[:, 0], timeseries[:, 1]+timeseries[:, 2],'-',alpha=0.5)
    # plt.plot(timeseries[:, 0], timeseries[:, 1] - timeseries[:, 2], '-', alpha=0.5)
    # plt.show()

    RC = RC[timeseries[:,0]!=0]
    timeseries = timeseries[timeseries[:,0]!=0,:]

    C = plt.imshow(timeseries)
    plt.colorbar(C)
    plt.show()

    return (time, RC, timeseries)

def write_profiles_to_nc(output_file, field_name, x, y, time, RC, timeseries):

    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('time',len(time))
    ds.createDimension('RC',len(RC))

    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = time

    rvar = ds.createVariable('RC', 'f4', ('RC',))
    rvar[:] = RC

    var = ds.createVariable(field_name, 'f4', ('RC','time'))
    var[:] = timeseries

    ds.longitude = x
    ds.latitude = y

    ds.close()


project_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund'

ecco_dir =  '/Volumes/helheim/Data_Repository/ECCO'
field_name = 'THETA'

location_name = 'Sound_Entrance'
location = (-21.379,70.091)

# location_name = 'Scoresby_Sund_Fjord'
# location = (-24.9012,71.1012)

# location_name = 'EGC1'
# location = (-11.1343,74.8293)

ecco_XC_tiles, ecco_YC_tiles, RC = read_ecco_geometry(ecco_dir)

row, col, tile, dist_from_point, x, y = \
    find_sample_point_in_ecco_tiles(location[0], location[1], ecco_XC_tiles, ecco_YC_tiles)

time, RC, timeseries = read_profiles_from_ECCO_extension_files(ecco_dir,field_name, row, col, tile,RC)

output_file = os.path.join(project_dir,'Data','Modeling','ECCO','ECCOv5_extension_'+field_name+'_'+location_name+'_profile_timeseries.nc')
write_profiles_to_nc(output_file, field_name, x, y, time, RC, timeseries)

