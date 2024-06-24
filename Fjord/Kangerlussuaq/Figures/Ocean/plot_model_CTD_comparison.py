
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

    points = list(file_dict.keys())
    points = [12]

    folder = os.path.join(axctd_folder, 'Processed', str(year))

    for point in points:
        file_list = file_dict[point]
        for file_name in file_list:
            if file_name[:4] == str(year):
                print('    - Reading from file '+file_name)

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

    return(theta_profiles, salt_profiles, longitude, latitude)


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


def read_profiles_from_L1_model(project_folder, year, lon_ctd, lat_ctd):
    output_file = os.path.join(project_folder, 'Data','Ocean', 'Modeling', 'L1_state_at_L2_CTDs.nc')
    ds = nc4.Dataset(output_file)
    depths = ds.variables['depths'][:]
    iterations = ds.variables['iterations'][:]
    longitudes = ds.variables['longitude'][:]
    latitudes = ds.variables['latitude'][:]
    Theta = ds.variables['Theta'][:, :, :]
    ds.close()

    dec_yrs = np.zeros((len(iterations),))
    for i in range(len(iterations)):
        dec_yrs[i] = iter_number_to_dec_yr(iterations[i], seconds_per_iter=300)

    dist = great_circle_distance(lon_ctd, lat_ctd, longitudes, latitudes)
    dist_index = np.where(dist==np.min(dist))[0][0]

    time_index = np.argmin(np.abs(dec_yrs-(year+243/365.25)))

    profile = np.column_stack([depths.ravel(), Theta[dist_index, time_index, :].ravel()])
    profile = profile[profile[:, 1] != 0, :]

    return(profile)


def read_profiles_from_L2_model(project_folder, year, lon_ctd, lat_ctd):
    output_file = os.path.join(project_folder, 'Data','Ocean', 'Modeling', 'L2_state_at_L2_CTDs.nc')
    ds = nc4.Dataset(output_file)
    depths = ds.variables['depths'][:]
    iterations = ds.variables['iterations'][:]
    longitudes = ds.variables['longitude'][:]
    latitudes = ds.variables['latitude'][:]
    Theta = ds.variables['Theta'][:, :, :]
    ds.close()

    dec_yrs = np.zeros((len(iterations),))
    for i in range(len(iterations)):
        dec_yrs[i] = iter_number_to_dec_yr(iterations[i], seconds_per_iter=60)

    dist = great_circle_distance(lon_ctd, lat_ctd, longitudes, latitudes)
    dist_index = np.where(dist==np.min(dist))[0][0]

    time_index = np.argmin(np.abs(dec_yrs-(year+213/365.25)))

    profile = np.column_stack([depths.ravel(), Theta[dist_index, time_index, :].ravel()])
    profile = profile[profile[:,1]!=0,:]

    return(profile)



def plot_theta_profiles(project_folder, theta_profile_ctd_set, theta_profile_L1_set, theta_profile_L2_set):

    min_y = 0
    max_y = 800

    min_T = -1.9
    max_T = 6

    fig = plt.figure(figsize=(8,10))

    plot_height = 1

    gs2 = GridSpec(plot_height*6, 3, left=0.09, right=0.98, bottom=0.05, top=0.95)

    ##################################################################################

    plot_depth = 300

    for i in range(len(theta_profile_ctd_set)):
        ax = fig.add_subplot(gs2[i*plot_height:(i+1)*plot_height, 0])

        for profile in theta_profile_ctd_set[i]:
            plt.plot(profile[:,1], profile[:,0], '-', color='silver', linewidth=1)
            set_int = interp1d(profile[:,0], profile[:,1])
            plot_temp = set_int(plot_depth)
            plt.plot(plot_temp, plot_depth, 'ko', markersize=6)
            plt.text(plot_temp, plot_depth, '{:.2f}'.format(plot_temp)+'$^{\circ}$C')

        ax.set_ylim([max_y, min_y])
        ax.set_xlim([min_T, max_T])

        ax.text(min_T, max_y, str(i+2016), ha='left', va='bottom')

        ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)

        if i==0:
            ax.set_title('AXCTDs')
        if i==2:
            ax.set_ylabel('Depth (m)')
        if i<5:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('Potential Temperature ($^{\circ}$C)')

    ##################################################################################

    for i in range(len(theta_profile_L1_set)):
        ax = fig.add_subplot(gs2[i * plot_height:(i + 1) * plot_height, 1])

        for profile in theta_profile_L1_set[i]:
            # print(np.shape(profile))
            plt.plot(profile[:, 1], profile[:, 0], '-', color='silver', linewidth=1)
            set_int = interp1d(profile[:, 0], profile[:, 1])
            plot_temp = set_int(plot_depth)
            plt.plot(plot_temp, plot_depth, 'ko', markersize=6)
            plt.text(plot_temp, plot_depth, '{:.2f}'.format(plot_temp) + '$^{\circ}$C')

        ax.set_ylim([max_y, min_y])
        ax.set_xlim([min_T, max_T])

        ax.text(min_T, max_y, str(i + 2016), ha='left', va='bottom')

        ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)

        if i == 0:
            ax.set_title('L1')
        if i == 2:
            ax.set_ylabel('Depth (m)')
        if i < 5:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('Potential Temperature ($^{\circ}$C)')

    ##################################################################################

    for i in range(len(theta_profile_L2_set)):
        ax = fig.add_subplot(gs2[i * plot_height:(i + 1) * plot_height, 2])

        for profile in theta_profile_L2_set[i]:
            # print(np.shape(profile))
            plt.plot(profile[:, 1], profile[:, 0], '-', color='silver', linewidth=1)
            set_int = interp1d(profile[:, 0], profile[:, 1])
            plot_temp = set_int(plot_depth)
            plt.plot(plot_temp, plot_depth, 'ko', markersize=6)
            plt.text(plot_temp, plot_depth, '{:.2f}'.format(plot_temp) + '$^{\circ}$C')

        ax.set_ylim([max_y, min_y])
        ax.set_xlim([min_T, max_T])

        ax.text(min_T, max_y, str(i + 2016), ha='left', va='bottom')

        ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)

        if i == 0:
            ax.set_title('L2')
        if i == 2:
            ax.set_ylabel('Depth (m)')
        if i < 5:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('Potential Temperature ($^{\circ}$C)')


    output_file = os.path.join(project_folder,'Figures','Ocean','AXCTD_Model_Comparison.png')
    plt.savefig(output_file)
    plt.close()






project_folder = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

axctd_folder = '/Users/mhwood/Documents/Research/Data Repository/Greenland/Ocean Properties/OMG CTDs'

file_dict = read_ctd_files_dict_from_list(project_folder)

theta_profile_ctd_set = []
theta_profile_L1_set = []
theta_profile_L2_set = []
for year in range(2016,2022):
    theta_profiles_ctd, salt_profiles_ctd, lon_ctd, lat_ctd = read_profile_from_repository(axctd_folder, file_dict, year)
    theta_profile_ctd_set.append(theta_profiles_ctd)

    theta_profile_model = read_profiles_from_L1_model(project_folder, year, lon_ctd, lat_ctd)
    theta_profile_L1_set.append([theta_profile_model])

    if year<2017:
        theta_profile_model = read_profiles_from_L2_model(project_folder, year, lon_ctd, lat_ctd)
        theta_profile_L2_set.append([theta_profile_model])


plot_theta_profiles(project_folder, theta_profile_ctd_set, theta_profile_L1_set, theta_profile_L2_set)











