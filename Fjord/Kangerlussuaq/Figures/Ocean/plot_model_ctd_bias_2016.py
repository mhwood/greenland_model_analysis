

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime
from matplotlib.gridspec import GridSpec
from scipy.interpolate import interp1d

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)


def iter_number_to_dec_yr(iter_number,seconds_per_iter=60):

    total_seconds = iter_number*seconds_per_iter
    date = datetime.datetime(1992,1,1) + datetime.timedelta(seconds=total_seconds)

    dec_yr = YMD_to_DecYr(date.year,date.month,date.day)
    # print(date)
    return(dec_yr)


def read_ctd_files_dict_from_list(project_folder):

    output_file = os.path.join(project_folder, 'Data', 'Ocean', 'AXCTDs', 'AXCTD List at Model Points.csv')

    file_dict = {}

    f = open(output_file)
    lines = f.read()
    f.close()
    lines = lines.split('\n')
    lines.pop(0)

    for line in lines:
        line = line.split(',')
        if len(line)>2:
            file_name = line[0]
            point = int(line[1])

            if point not in list(file_dict.keys()):
                file_dict[point] = [file_name]
            else:
                file_dict[point].append(file_name)

    return(file_dict)


def read_profile_from_repository(axctd_folder, file_dict, year):

    theta_profiles = []
    salt_profiles = []

    longitudes = []
    latitudes = []

    #points = list(file_dict.keys())
    points = [3, 5, 8, 9, 10, 11, 12]

    folder = os.path.join(axctd_folder, 'Processed', str(year))

    for point in points:
        # print(point)
        file_list = file_dict[point]
        if point==11:
            file_list=['20161006_141711']
        # print(file_list[0])
        for file_name in file_list:
            if file_name[:4] == str(year):
                print('    - Reading from file '+file_name +' (point '+str(point)+')')

                if 'CTD_'+file_name+'.nc' in os.listdir(folder):
                    ds = nc4.Dataset(os.path.join(folder,'CTD_'+file_name+'.nc'))
                    depth = ds.variables['depth'][:]
                    temp = ds.variables['potential_temperature'][:]
                    salt = ds.variables['practical_salinity'][:]
                    longitude = ds.longitude
                    latitude = ds.latitude
                    # print(longitude,latitude)
                    ds.close()

                    theta_profiles.append(np.column_stack([depth, temp]))
                    salt_profiles.append(np.column_stack([depth, salt]))
                    longitudes.append(longitude)
                    latitudes.append(latitude)

    return(points, theta_profiles, salt_profiles, longitudes, latitudes)


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


def read_profiles_from_L2_model(project_folder, year, lon_ctd, lat_ctd):
    output_file = os.path.join(project_folder, 'Data','Ocean', 'Modeling', 'L2_state_at_L2_CTDs.nc')
    ds = nc4.Dataset(output_file)
    depths = ds.variables['depths'][:]
    iterations = ds.variables['iterations'][:]
    longitudes = ds.variables['longitude'][:]
    latitudes = ds.variables['latitude'][:]
    Theta = ds.variables['Theta'][:, :, :]
    Salt = ds.variables['Salt'][:, :, :]
    ds.close()

    dec_yrs = np.zeros((len(iterations),))
    for i in range(len(iterations)):
        dec_yrs[i] = iter_number_to_dec_yr(iterations[i], seconds_per_iter=60)

    dist = great_circle_distance(lon_ctd, lat_ctd, longitudes, latitudes)
    dist_index = np.where(dist==np.min(dist))[0][0]

    time_index = np.argmin(np.abs(dec_yrs-(year+213/365.25)))

    theta_profile = np.column_stack([depths.ravel(), Theta[dist_index, time_index, :].ravel()])
    theta_profile = theta_profile[theta_profile[:,1]!=0,:]

    salt_profile = np.column_stack([depths.ravel(), Salt[dist_index, time_index, :].ravel()])
    salt_profile = salt_profile[salt_profile[:, 1] != 0, :]

    return(theta_profile, salt_profile)



def plot_model_bias(project_folder, var_name, profiles_ctd, profiles_model, points):

    min_y = 0
    max_y = 600

    min_T = -1.9
    max_T = 6

    fig = plt.figure(figsize=(8, 10))

    plot_height = 1

    gs2 = GridSpec(plot_height * 4, 4, left=0.09, right=0.98, bottom=0.05, top=0.95)

    for p in range(len(profiles_ctd)):

        # plot the ctd and model profiles

        if p<3:
            ax = fig.add_subplot(gs2[(p + 1) * plot_height:(p + 2) * plot_height, 0])
        else:
            ax = fig.add_subplot(gs2[(p-3) * plot_height:(p -2) * plot_height, 2])

        profile_ctd = profiles_ctd[p]
        ax.plot(profile_ctd[:, 1], profile_ctd[:, 0], '-', color='k', linewidth=2)

        profile_model = profiles_model[p]
        ax.plot(profile_model[:, 1], profile_model[:, 0], '-', color='g', linewidth=2)

        ax.set_ylim([max_y, min_y])
        ax.set_xlim([min_T, max_T])

        # ax.text(min_T, max_y, str(i + 2016), ha='left', va='bottom')

        ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)

    output_file = os.path.join(project_folder, 'Figures', 'Ocean', 'Model_Bias_Comparison_'+str(var_name)+'_2016.png')
    plt.savefig(output_file)
    plt.close()


project_folder = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

axctd_folder = '/Users/mhwood/Documents/Research/Data Repository/Greenland/Ocean Properties/OMG CTDs'

file_dict = read_ctd_files_dict_from_list(project_folder)

points, theta_profiles_ctd, salt_profiles_ctd, lon_ctds, lat_ctds =\
    read_profile_from_repository(axctd_folder, file_dict, year=2016)

theta_profiles_model = []
salt_profiles_model = []

for ll in range(len(lon_ctds)):
    theta, salt = read_profiles_from_L2_model(project_folder, 2016, lon_ctds[ll], lat_ctds[ll])
    theta_profiles_model.append(theta)
    salt_profiles_model.append(salt)

print(theta_profiles_model)

plot_model_bias(project_folder, 'Theta', theta_profiles_ctd, theta_profiles_model, points)


