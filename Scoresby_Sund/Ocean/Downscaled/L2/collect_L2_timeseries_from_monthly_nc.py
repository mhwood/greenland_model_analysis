
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime
from scipy.interpolate import interp1d
import shapefile as sf

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_transect_file(config_dir):

    model_level = 'L3'
    model_name = 'L3_Scoresby_Sund'
    transect_file = os.path.join(config_dir,model_level,model_name,'results_iceplume_ks22',
                               'Scoresby_Sund_transects','Scoresby_Sund_transects.nc')

    ds = nc4.Dataset(transect_file)
    years = ds.variables['years'][:]
    months = ds.variables['months'][:]
    depth = ds.variables['depth'][:]
    distance = ds.variables['transect_distance'][:]
    density = ds.variables['transect_density'][:, :, :]
    temperature = ds.variables['transect_temp'][:, :, :]
    ds.close()

    distance *= 1e-3

    time = np.zeros(np.shape(density)[0])
    for t in range(len(time)):
        time[t] = YMD_to_DecYr(int(years[t]), int(months[t]), 15)

    return(time, depth, density, temperature)

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

def read_transect_points_from_shp(file_path):
    r = sf.Reader(file_path)
    shapes = r.shapes()
    points = np.array(shapes[0].points)

    points = series_to_N_points(points,1000)

    XC = points[:,0]
    YC = points[:,1]

    d = np.zeros_like(XC)
    for i in range(len(d)-1):
        d[i+1] = d[i] + great_circle_distance(XC[i],YC[i],XC[i+1],YC[i+1])

    return(XC, YC, d)

def find_closest_transect_point(lon_ref, lat_ref, XC, YC):

    dist = great_circle_distance(lon_ref, lat_ref, XC, YC)
    index = np.argmin(dist)
    return(index)

def average_temp_transect_in_density_range(point_index, density, temperature, min_density, max_density):
    temp_timeseries = np.zeros((np.shape(density)[0],))

    interp_densities = np.linspace(min_density,max_density,100)

    for t in range(np.shape(density)[0]):
        # if t % 10 == 0:
        #     print('    - Calculating step ' + str(t) + ' of ' + str(np.shape(density)[0]))
        density_profile = density[t, :, point_index]
        temperature_profile = temperature[t, :, point_index]

        if np.min(density_profile) < max_density and np.max(density_profile) > min_density:
            # set_int = interp1d(density_profile[indices], depth[indices])
            temp_set_int = interp1d(density_profile, temperature_profile)
            temp_timeseries[t] = np.mean(temp_set_int(interp_densities))
        else:
            temp_timeseries[t] = np.NaN

    return(temp_timeseries)

def plot_temperature_timeseries():
    a=1
    # plt.subplot(2,1,1)
    # C = plt.pcolormesh(time,depth,temperature[:,:,point_index].T, cmap='turbo')
    # plt.contour(time,depth,density[:,:,point_index].T,levels=[1029,1030.1],colors='k')
    # plt.colorbar(C)
    # plt.gca().invert_yaxis()
    # plt.gca().set_ylim([1000,0])
    #
    # plt.subplot(2,1,2)
    # plt.plot(time, temp_timeseries)
    #
    # plt.show()

def store_temperature_timeseries_as_nc(output_dir,output_file,total_time,subsets,theta_timeseries_set,densities):

    ds = nc4.Dataset(os.path.join(output_dir,output_file),'w')
    ds.createDimension('time',len(total_time))

    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = total_time

    for s in range(len(subsets)):
        grp = ds.createGroup(subsets[s])
        svar = grp.createVariable('THETA','f4',('time',))
        svar[:] = theta_timeseries_set[s]
        grp.min_density = densities[s][0]
        grp.max_density = densities[s][1]
    ds.close()


config_dir = '/Volumes/helheim/Ocean_Modeling/Projects/Downscale_Greenland/' \
             'MITgcm/configurations/downscaled_greenland'

output_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/Data/Modeling/Downscaled/L3_Scoresby_Sund'
output_file = 'L3_Scoresby_Sund_THETA_Timeseries.nc'

point_fjord_entrance = [-21.379,70.091]
point_mid_fjord = [-24.888,71.094]
point_DJG = [-28.35017,71.93849]

densities = [[1029,1030.1],[1029,1030.1],[1029,1029.65]]
subsets = ['near_glacier','sound_entrance','mid_sound']
points = [point_DJG,point_fjord_entrance,point_mid_fjord]

time, depth, density, temperature = read_transect_file(config_dir)

theta_timeseries_set = []

project_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/'
shp_file_name = 'Scoresby_Sund_Transect'
shp_file_path = os.path.join(project_dir,'Map','Shapefiles','Analysis_Model_Sample_Areas',shp_file_name)
XC, YC, transect_distance = read_transect_points_from_shp(shp_file_path)

for s in range(len(subsets)):
    lon_ref = points[s][0]
    lat_ref = points[s][1]

    min_density = densities[s][0]
    max_density = densities[s][1]

    point_index = find_closest_transect_point(lon_ref, lat_ref, XC, YC)
    print(point_index)
    temp_timeseries = average_temp_transect_in_density_range(point_index, density, temperature, min_density, max_density)

    theta_timeseries_set.append(temp_timeseries)

store_temperature_timeseries_as_nc(output_dir,output_file,time,subsets,theta_timeseries_set,densities)



