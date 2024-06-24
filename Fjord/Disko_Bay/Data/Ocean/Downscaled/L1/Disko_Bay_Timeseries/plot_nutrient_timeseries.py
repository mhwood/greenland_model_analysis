

import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from scipy.interpolate import interp1d
import datetime


def find_profile_files_in_box(folder, min_lon, min_lat, max_lon, max_lat):

    file_names = []
    file_dates = []

    f = open(os.path.join(folder,'Map','Shapefiles','Nitrate_OSD_locations.csv'))
    lines = f.read()
    f.close()
    lines = lines.split('\n')
    lines.pop(0)

    bbox = np.array([[min_lon, min_lat],
                     [max_lon, min_lat],
                     [max_lon, max_lat],
                     [min_lon, max_lat],
                     [min_lon, min_lat]])

    p = mplPath.Path(bbox)

    for line in lines:
        line = line.split(',')
        if p.contains_point([float(line[1]), float(line[2])]):
            file_names.append(line[0])
            file_dates.append(line[3])

    return(file_names, file_dates)


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
        # plt.title(file_name)
        # plt.show()

        profiles.append(profile)

    return(profiles)


def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)


def compute_averaged_profiles(common_depth, dates, profiles, min_year, max_year):

    date_indices = []
    for p in range(len(dates)):
        if int(dates[p][:4])>=min_year and int(dates[p][:4])<=max_year:
            date_indices.append(p)

    profile_mean = []
    profile_std = []
    for d in range(len(common_depth)):
        depth_points = []
        for index in date_indices:
            if ~np.isnan(profiles[index][d]):
                depth_points.append(profiles[index][d])
        if len(depth_points)>1:
            profile_mean.append(np.mean(depth_points))
            profile_std.append(np.std(depth_points))
        else:
            profile_mean.append(np.nan)
            profile_std.append(np.nan)

    profile_mean = np.array(profile_mean)
    profile_std = np.array(profile_std)

    # plt.plot(common_depth,profile_mean,'k-')
    # plt.plot(common_depth, profile_mean+profile_std,'k--')
    # plt.show()

    return(profile_mean, profile_std)


def create_timeseries_plot(folder, common_depth, dates, profiles):

    fig = plt.figure(figsize=(10,5))

    day_spacing = 45

    for d in range(len(dates)):
        dec_yr = YMD_to_DecYr(int(dates[d][:4]),int(dates[d][4:6]),int(dates[d][6:8]))
        t = np.array([dec_yr-day_spacing/365.25, dec_yr, dec_yr+day_spacing/365.25])
        T, D = np.meshgrid(t, common_depth)
        V = np.column_stack([profiles[d], profiles[d], profiles[d]])
        plt.pcolormesh(T,D,V, vmin=0, vmax=15, cmap='jet')

    plt.gca().set_xlim(1980,2015)
    plt.gca().set_ylim([np.max(common_depth),0])

    plt.savefig(os.path.join(folder,'Figures','Nutrients','Disko_Bay_Nitrate_Timeseries.png'))
    plt.close(fig)


def create_stacked_profile_plot(folder, common_depth, dates, profiles):

    fig = plt.figure(figsize=(10,5))

    day_spacing = 45
    min_yr = 1980
    max_yr = 2015

    sorted_indices = sorted(range(len(dates)), key=lambda k: dates[k])

    for d in sorted_indices:
        dec_yr = YMD_to_DecYr(int(dates[d][:4]),int(dates[d][4:6]),int(dates[d][6:8]))
        rgb = plt.cm.jet((np.clip(dec_yr, min_yr, max_yr) - min_yr) / (max_yr-min_yr))
        plt.plot(profiles[d],common_depth,color=rgb)

    plt.gca().set_xlim(0, 15)
    plt.gca().set_ylim([np.max(common_depth),0])

    plt.savefig(os.path.join(folder,'Figures','Nutrients','Disko_Bay_Nitrate_Profiles.png'))
    plt.close(fig)


def create_mean_profile_plot(folder, common_depth,
                             profile_1980_2000_mean, profile_1980_2000_std,
                             profile_2000_2020_mean, profile_2000_2020_std):

    fig = plt.figure(figsize=(10,5))

    plt.plot(profile_1980_2000_mean, common_depth, 'b-')
    plt.plot(profile_1980_2000_mean-profile_1980_2000_std, common_depth, 'b--')
    plt.plot(profile_1980_2000_mean+profile_1980_2000_std, common_depth, 'b--')

    plt.plot(profile_2000_2020_mean, common_depth, 'r-')
    plt.plot(profile_2000_2020_mean-profile_2000_2020_std, common_depth, 'r--')
    plt.plot(profile_2000_2020_mean+profile_2000_2020_std, common_depth, 'r--')

    plt.gca().set_xlim(0, 15)
    plt.gca().set_ylim([np.max(common_depth),0])

    plt.savefig(os.path.join(folder,'Figures','Nutrients','Disko_Bay_Nitrate_Mean_Profile.png'))
    plt.close(fig)




folder = '/Users/mike/Documents/Research/Projects/Disko Bay'

# just the trough in Disko
min_lat = 68.8750
min_lon = -53.3744
max_lat = 68.9840
max_lon = -53.1117

# all of Disko
# min_lat = 68.5252
# min_lon = -53.8670
# max_lat = 69.8216
# max_lon = -50.5690

common_depth = np.arange(700)

file_list, dates = find_profile_files_in_box(folder, min_lon, min_lat, max_lon, max_lat)

profiles = read_profiles_from_nc(folder, file_list, common_depth)

profile_1980_2000_mean, profile_1980_2000_std = compute_averaged_profiles(common_depth, dates, profiles, 1980, 2000)
profile_2000_2020_mean, profile_2000_2020_std = compute_averaged_profiles(common_depth, dates, profiles, 2000, 2020)

# create_timeseries_plot(folder, common_depth, dates, profiles)

# create_stacked_profile_plot(folder, common_depth, dates, profiles)

create_mean_profile_plot(folder, common_depth,
                         profile_1980_2000_mean, profile_1980_2000_std,
                         profile_2000_2020_mean, profile_2000_2020_std)
