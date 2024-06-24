import datetime
import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import shapefile
from pyproj import Transformer
from scipy.interpolate import interp1d, griddata

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1,run_test = True):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(polygon_array[:, x_column], polygon_array[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
        run_test = False
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon

def series_to_N_points(series,N):
    #find the total length of the series
    totalDistance=0
    for s in range(len(series[:,0])-1):
        totalDistance+=((series[s,0]-series[s+1,0])**2+(series[s,1]-series[s+1,1])**2)**0.5
    intervalDistance=totalDistance/(N-1)

    #make the list of points
    newSeries=series[0,:]
    currentS = 0
    currentPoint1=series[currentS,:]
    currentPoint2=series[currentS+1,:]
    for p in range(N-2):
        distanceAccrued = 0
        while distanceAccrued<intervalDistance:
            currentLineDistance=((currentPoint1[0]-currentPoint2[0])**2+(currentPoint1[1]-currentPoint2[1])**2)**0.5
            if currentLineDistance<intervalDistance-distanceAccrued:
                distanceAccrued+=currentLineDistance
                currentS+=1
                currentPoint1 = series[currentS, :]
                currentPoint2 = series[currentS + 1, :]
            else:
                distance=intervalDistance-distanceAccrued
                newX=currentPoint1[0]+(distance/currentLineDistance)*(currentPoint2[0]-currentPoint1[0])
                newY = currentPoint1[1] + (distance / currentLineDistance) * (currentPoint2[1] - currentPoint1[1])
                distanceAccrued=intervalDistance+1
                newSeries=np.vstack([newSeries,np.array([newX,newY])])
                currentPoint1=np.array([newX,newY])
    newSeries = np.vstack([newSeries, series[-1,:]])
    return(newSeries)

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

def read_transect_from_shapefile(folder, transect_number):
    r = shapefile.Reader(os.path.join(folder,'Map','Shapefiles','Davis_Strait_Transect_'+str(transect_number)))
    points = np.array(r.shapes()[0].points).astype(float)
    r.close()
    points = series_to_N_points(points, N=1000)
    return(points)

def compute_diastance_and_bathymetry_profile(folder,transect):

    transect_4326 = reproject_polygon(transect,3413,4326)

    distance = np.zeros((np.shape(transect)[0],))
    for i in range(len(distance)-1):
        distance[i+1] = distance[i] + great_circle_distance(transect_4326[i,0], transect_4326[i,1],
                                                            transect_4326[i+1,0], transect_4326[i+1,1])

    ds = nc4.Dataset(os.path.join(folder,'Data','Bathymetry','Davis_Strait_Bathymetry.nc'))
    bathy_lon = ds.variables['lon'][:]
    bathy_lat = ds.variables['lat'][:]
    bathy_grid = -1*np.array(ds.variables['elevation'][:,:])
    ds.close()

    skip = 20
    bathy_lon = bathy_lon[::skip]
    bathy_lat = bathy_lat[::skip]
    bathy_grid = bathy_grid[::skip,::skip]

    Lon, Lat = np.meshgrid(bathy_lon,bathy_lat)
    points = np.column_stack([Lon.ravel(), Lat.ravel()])

    bathymetry = griddata(points, bathy_grid.ravel(),(transect_4326[:,0], transect_4326[:,1]))

    return(distance, bathymetry, transect_4326)

def find_profile_files_on_transect(folder,transect,year):

    file_names = []
    lons = []
    lats = []
    file_dates = []

    f = open(os.path.join(folder,'Map','Shapefiles','Nitrate_OSD_locations.csv'))
    lines = f.read()
    f.close()
    lines = lines.split('\n')
    lines.pop(0)

    for line in lines:
        line = line.split(',')
        if line[3][:4]==str(year):
            file_names.append(line[0])
            lons.append(line[1])
            lats.append(line[2])
            file_dates.append(line[3])

    points = np.column_stack([np.array(lons),np.array(lats)])
    points_3413 = reproject_polygon(points,4326,3413)
    points_3413 = points_3413.astype(float)

    output_file_names = []
    output_x = []
    output_y = []
    output_dates = []

    distance_threshold = 5000
    for p in range(len(file_names)):
        dist = ((transect[:,0] - points_3413[p,0])**2 + (transect[:,1] - points_3413[p,1])**2)**0.5
        if np.min(dist)<distance_threshold:
            output_file_names.append(file_names[p])
            output_x.append(points_3413[p,0])
            output_y.append(points_3413[p,1])
            output_dates.append(file_dates[p])

    return(output_file_names, output_x, output_y, output_dates)

def read_profiles_from_nc(folder, file_names, common_depth, variable_name='Nitrate', subset='OSD'):

    profiles = []

    for file_name in file_names:
        file_path = os.path.join(folder,'Data',variable_name,subset,file_name)
        ds = nc4.Dataset(file_path)
        depth = ds.variables['z'][:]
        var_points = ds.variables[variable_name][:]
        ds.close()

        indices = var_points>=0
        depth = depth[indices]
        var_points = var_points[indices]

        if len(depth)<=1:
            profile = np.nan*np.ones_like(common_depth)
        else:
            set_int = interp1d(depth,var_points,fill_value=np.nan,bounds_error=False)
            profile = set_int(common_depth)

        # plt.plot(profile,common_depth,'k-')
        # plt.plot(var_points,depth,'g.')
        # plt.show()

        profiles.append(profile)

    return(profiles)

def interpolated_profiles_to_common_grid(transect, distance, common_depth, profile_x, profile_y, profiles):

    distance_indices = []
    for p in range(len(profiles)):
        dist = ((transect[:,0]-profile_x[p])**2 + (transect[:,1]-profile_y[p])**2)**0.5
        dist_index = np.argmin(dist)
        distance_indices.append(dist_index)

    sorted_indices = sorted(range(len(distance_indices)), key=lambda k: distance_indices[k])

    d = []
    V_started = False
    for index in sorted_indices:
        d.append(distance[distance_indices[index]])
        if not V_started:
            V_started = True
            V = np.reshape(profiles[index], (len(profiles[index]),1))
        else:
            V = np.hstack([V,np.reshape(profiles[index], (len(profiles[index]),1))])

    D, Z = np.meshgrid(np.array(d), common_depth)
    points = np.column_stack([D.ravel(), Z.ravel()])

    new_D, Z = np.meshgrid(distance, common_depth)

    interpolated_profile_grid = griddata(points,V.ravel(),(new_D, Z))

    return(interpolated_profile_grid)

def outputting_transects_to_nc(folder, transect_number, output_years, output_months, output_days,
                               common_depth, distance, bathymetry, transect, transect_4326, interpolation_grid):

    ds = nc4.Dataset(os.path.join(folder,'Data','Nitrate','Transect_'+str(transect_number)+'_Nitrate_Profile.nc'),'w')

    ds.createDimension('time',len(output_years))
    ds.createDimension('depth',len(common_depth))
    ds.createDimension('distance', len(distance))

    var = ds.createVariable('year','i4',('time',))
    var[:] = output_years

    var = ds.createVariable('month', 'i4', ('time',))
    var[:] = output_months

    var = ds.createVariable('day', 'i4', ('time',))
    var[:] = output_days

    var = ds.createVariable('depth', 'f4', ('depth',))
    var[:] = common_depth

    var = ds.createVariable('distance', 'f4', ('distance',))
    var[:] = distance

    var = ds.createVariable('bathymetry', 'f4', ('distance',))
    var[:] = bathymetry

    var = ds.createVariable('x', 'f4', ('distance',))
    var[:] = transect[:, 0]

    var = ds.createVariable('y', 'f4', ('distance',))
    var[:] = transect[:, 1]

    var = ds.createVariable('longitude', 'f4', ('distance',))
    var[:] = transect_4326[:, 0]

    var = ds.createVariable('latitude', 'f4', ('distance',))
    var[:] = transect_4326[:, 1]

    var = ds.createVariable('Nitrate', 'f4', ('time','depth','distance'))
    var[:,:,:] = interpolation_grid

    ds.close()

def create_transect_plot(distance, bathymetry, common_depth, interpolated_profile_grid):

    plt.pcolormesh(distance,common_depth,interpolated_profile_grid, cmap='jet')
    plt.plot(distance, bathymetry, 'k-')
    plt.gca().invert_yaxis()
    plt.show()
    a=1


folder = '/Users/mike/Documents/Research/Projects/Disko Bay'
transect_number = 3

print('  - Retrieving the transect')
# step 1: get transect
transect = read_transect_from_shapefile(folder, transect_number)
if transect[0,0]>transect[-1,0]:
    transect = np.flipud(transect)

print('  - Computing the distance and bathymetry')
distance, bathymetry, transect_4326 = compute_diastance_and_bathymetry_profile(folder,transect)
common_depth = np.arange(int(np.max(bathymetry))+10)

print('  - Computing interpolated transects')
# years = [2009, 2010]
years = np.arange(1992,2024)

interpolation_grid_started = False

output_years = []
output_months = []
output_days = []
for year in years:
    # print('      - Working on year '+str(year))

    # step 2: get list of files
    file_list, x, y, dates = find_profile_files_on_transect(folder,transect,year)
    # print('      - Found '+str(len(file_list))+' profile files')

    if len(file_list)>1:
        print('        - Adding '+str(len(file_list))+' from year '+str(year))
        output_years.append(year)

        all_dates = []
        for date in dates:
            all_dates.append(datetime.datetime(year,int(date[4:6]), int(date[6:8])))
        ref_date = datetime.datetime(1900, 1, 1)
        mean_date =  ref_date + sum([date - ref_date for date in all_dates], datetime.timedelta()) / len(all_dates)
        # mean_date = np.mean(all_dates)
        output_months.append(mean_date.month)
        output_days.append(mean_date.day)

        # step 3: read in the files to a common depth
        profiles = read_profiles_from_nc(folder, file_list, common_depth)

        # step 4: make interpolated profile files
        interpolated_profile_grid = interpolated_profiles_to_common_grid(transect, distance, common_depth, x, y, profiles)

        if not interpolation_grid_started:
            interpolation_grid = np.reshape(interpolated_profile_grid,
                                            (1,np.shape(interpolated_profile_grid)[0],np.shape(interpolated_profile_grid)[1]))
            interpolation_grid_started=True
        else:
            interpolation_grid = np.concatenate([interpolation_grid,
                                           np.reshape(interpolated_profile_grid,
                                                      (1, np.shape(interpolated_profile_grid)[0],
                                                       np.shape(interpolated_profile_grid)[1]))], axis=0)
        # print(np.shape(interpolation_grid))
    else:
        pass
        # print('         - Skipping this year with too few profiles')

# create_transect_plot(distance, bathymetry, common_depth, interpolated_profile_grid)

outputting_transects_to_nc(folder, transect_number, output_years, output_months, output_days,
                           common_depth, distance, bathymetry, transect, transect_4326, interpolation_grid)


