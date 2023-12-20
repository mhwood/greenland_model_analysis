
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import datetime
from scipy.interpolate import interp1d

def read_observation_metadata(project_folder):
    file_path = os.path.join(project_folder,'Data','In Situ','Hadley','Metadata','Hadley_CTD_Locations.csv')
    f = open(file_path)
    lines = f.read()
    f.close()
    lines = lines.split('\n')
    lines.pop(0)
    file_names = []
    lons = []
    lats = []
    for line in lines:
        line = line.split(',')
        if len(line)>3:
            file_names.append(line[0])
            lons.append(float(line[7]))
            lats.append(float(line[8]))

    lons = np.array(lons)
    lats = np.array(lats)

    return(file_names, lons, lats)

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

def subset_file_list_close_to_point(file_names,lons,lats,ref_lon, ref_lat, threshold):
    dist = great_circle_distance(ref_lon, ref_lat, lons,lats)
    indices = np.where(dist<threshold)[0]
    file_names_subset = []
    for index in indices:
        file_names_subset.append(file_names[index])
    return(file_names_subset)

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_temperature_timeseries_from_profiles(project_dir, file_names_subset, min_depth, max_depth):

    profile_dir = os.path.join(project_dir, 'Data', 'In Situ', 'Hadley', 'Data')

    dec_yrs = []
    temperature_means = []

    interp_depths = np.arange(min_depth,max_depth)

    for file_name in file_names_subset:
        year = file_name[:4]
        ds = nc4.Dataset(os.path.join(profile_dir,year,'CTD_'+file_name+'.nc'))
        depth = ds.variables['depth'][:]
        temp = ds.variables['potential_temperature'][:]
        year = ds.year
        month = ds.month
        day = ds.day
        hour = ds.hour
        minute = int(str(ds.minute)[:2])
        # print(year,month,day,hour,minute)
        if day==32:
            day = 31
        dec_yr = YMD_to_DecYr(year, month, day, hour, minute)
        ds.close()
        # print(file_name,np.min(depth),np.max(depth))
        if np.min(depth)<=min_depth and np.max(depth)>=max_depth:
            print(file_name,np.min(depth),np.max(depth))
            dec_yrs.append(dec_yr)
            set_int = interp1d(depth,temp)
            interp_temp = set_int(interp_depths)
            temperature_means.append(np.mean(interp_temp))

    timeseries = np.column_stack([np.array(dec_yrs),
                                  np.array(temperature_means)])

    return(timeseries)

def output_timeseries_to_nc(output_file,timeseries,min_depth,max_depth):

    ds = nc4.Dataset(output_file,'w')
    ds.createDimension('time', np.shape(timeseries)[0])

    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = timeseries[:,0]
    tvar = ds.createVariable('temperature_mean', 'f4', ('time',))
    tvar[:] = timeseries[:, 1]
    # tvar = ds.createVariable('temperature_std', 'f4', ('time',))
    # tvar[:] = timeseries[:, 2]

    ds.min_depth = min_depth
    ds.max_depth = max_depth

    ds.close()



project_folder = '/Users/michwood/Documents/Research/Projects/Scoresby Sund'

# step 1: read in the data list as file names, lon, lat
file_names, lons, lats = read_observation_metadata(project_folder)

# step 2: get the profiles within a threshold of this point
lon_ref = -21.070
lat_ref = 69.369
threshold = 40000
file_names_subset = subset_file_list_close_to_point(file_names,lons,lats,lon_ref, lat_ref, threshold)

min_depth = 200
max_depth = 400

# step 3: store as a timeseries
timeseries = read_temperature_timeseries_from_profiles(project_folder, file_names_subset, min_depth, max_depth)

output_file = os.path.join(project_folder,'Data','In Situ','Hadley','Data','Hadley_Shelf_Temperature_Timeseries.nc')
output_timeseries_to_nc(output_file,timeseries,min_depth,max_depth)

