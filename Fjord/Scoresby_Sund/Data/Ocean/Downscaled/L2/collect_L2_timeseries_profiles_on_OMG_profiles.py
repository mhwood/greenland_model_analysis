
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime
from scipy.interpolate import interp1d
import shapefile as sf


def get_nominal_dates(subset,profile_subset):
    if profile_subset=='OMG':
        if subset=='near_glacier':
            date_list = ['20190819_134458','20180827_143507','20171016_121232']#'20210820_131357','20200904_183011',
        if subset=='sound_entrance':
            date_list = ['20190819_134458','20180825_140502','20171015_123129']#'20210825_125131','20200904_155014',
    else:
        date_list = []
        for year in range(2001,2019):
            for month in range(1,13):
                date_str = str(year)+'{:02d}'.format(month)+'15_120000'
                date_list.append(date_str)
    return(date_list)

def get_lists_of_src_files(data_folder, date_strs):

    src_list = []

    for date_str in date_strs:
        year = int(date_str[:4])
        month = int(date_str[4:6])
        day = int(date_str[6:8])

        date = datetime.datetime(year, month, day)
        if month==12:
            date_str_1 = str(year)+'{:02d}'.format(month)
            date_str_2 = str(year+1) + '{:02d}'.format(1)
        else:
            date_str_1 = str(year) + '{:02d}'.format(month)
            date_str_2 = str(year) + '{:02d}'.format(month + 1)
        # print(date_str, date_str_1, date_str_2)

        file_name_1 = ''
        file_name_2 = ''

        for file_name in os.listdir(data_folder):
            if '.'+date_str_1+'.' in file_name and file_name[0]!='.':
                file_name_1 = file_name
            if '.'+date_str_2+'.' in file_name and file_name[0]!='.':
                file_name_2 = file_name

        # print(date_str,file_name_1,file_name_2)

        if month != 12:
            date_1 = datetime.datetime(year, month, 1)
            date_2 = datetime.datetime(year, month+1, 1)
        else:
            date_1 = datetime.datetime(year+1, month, 1)
            date_2 = datetime.datetime(year, 1, 1)

        scalar = (date-date_1)/(date_2-date_1)

        if file_name_1!='' and file_name_2!='':
            src_list.append([file_name_1,file_name_2,scalar])

    return(src_list)

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

def interpolate_profile_from_output(data_folder, point, src_list):

    file_name_1 = src_list[0]
    file_name_2 = src_list[1]
    scalar = src_list[2]

    ds = nc4.Dataset(os.path.join(data_folder,file_name_1))
    theta_1 = ds.variables['Theta'][:,:,:,:]
    salt_1 = ds.variables['Salt'][:,:,:,:]
    lon = ds.variables['longitude'][:,:]
    lat = ds.variables['latitude'][:, :]
    depth = ds.variables['depths'][:]
    ds.close()

    ds = nc4.Dataset(os.path.join(data_folder,file_name_2))
    theta_2 = ds.variables['Theta'][:, :, :, :]
    salt_2 = ds.variables['Salt'][:, :, :, :]
    ds.close()

    theta = (1 - scalar) * theta_1 + scalar * theta_2
    salt = (1 - scalar) * salt_1 + scalar * salt_2

    dist = great_circle_distance(point[0], point[1], lon, lat)
    row, col = np.where(dist==np.min(dist))
    row = row[0]
    col = col[0]

    theta_profile = np.column_stack([depth, theta[0,:,row,col]])
    salt_profile = np.column_stack([depth, salt[0, :, row, col]])

    return(depth, theta_profile,salt_profile)



def store_profiles_to_nc(output_dir,output_file,points,subsets,depth,
                         all_times, all_theta_profiles, all_salt_profiles):

    ds = nc4.Dataset(os.path.join(output_dir,output_file),'w')
    ds.createDimension('depth',len(all_theta_profiles[0][0][:, 0]))

    dvar = ds.createVariable('depth', 'f4', ('depth',))
    dvar[:] = all_theta_profiles[0][0][:, 0]

    for s in range(len(subsets)):
        grp = ds.createGroup(subsets[s])

        if s == 0:
            dvar = grp.createVariable('depth', 'f4', ('depth',))
            dvar[:] = all_theta_profiles[s][0][:, 0]

        grp.lon = points[0]
        grp.lat = points[1]

        for d in range(len(all_theta_profiles[s])):
            sub_grp = grp.createGroup(str(all_times[s][d]))
            tvar = sub_grp.createVariable('Theta', 'f4', ('depth',))
            tvar[:] = all_theta_profiles[s][d][:, 1]
            svar = sub_grp.createVariable('Salt', 'f4', ('depth',))
            svar[:] = all_salt_profiles[s][d][:, 1]


    ds.close()


config_dir = '/Volumes/helheim/Ocean_Modeling/Projects/Downscale_Greenland/' \
             'MITgcm/configurations/downscaled_greenland'
data_folder = os.path.join(config_dir,'L3','L3_Scoresby_Sund','results_iceplume_ks22','state_3D_mon_mean')

output_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/Data/Modeling/Downscaled/L3_Scoresby_Sund'

profile_subset = 'OMG'
profile_subset = 'All'

output_file = 'L3_Scoresby_Sund_'+profile_subset+'_Profiles.nc'

point_fjord_entrance = [-21.379,70.091]
point_mid_fjord = [-24.888,71.094]
point_DJG = [-28.35017,71.93849]

subsets = ['near_glacier','sound_entrance']
points = [point_DJG,point_fjord_entrance]
all_times = []
all_theta_profiles = []
all_salt_profiles = []

project_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/'

for s in range(len(subsets)):
    subset = subsets[s]

    print('  - Collecting the profiles in the '+subset+' subset')

    print('       - Getting a list of dates')
    date_strs = get_nominal_dates(subset,profile_subset)

    print('       - Getting the files corresponding to these dates')
    src_lists = get_lists_of_src_files(data_folder, date_strs)

    print('       - Interpolating profile from output')
    point = points[s]
    time = []
    theta_profiles = []
    salt_profiles = []
    for s in range(len(src_lists)):
        src_list = src_lists[s]
        depth, theta_profile, salt_profile = interpolate_profile_from_output(data_folder, point, src_list)
        time.append(int(date_strs[s][:4]))
        theta_profiles.append(theta_profile)
        salt_profiles.append(salt_profile)

    all_times.append(time)
    all_theta_profiles.append(theta_profiles)
    all_salt_profiles.append(salt_profiles)

    # plt.subplot(1,2,1)
    # plt.plot(theta_profile[:,1],theta_profile[:,0])
    #
    # plt.subplot(1, 2, 2)
    # plt.plot(salt_profile[:, 1], salt_profile[:, 0])
    #
    # plt.show()

print('   - Outputting profiles to nc')
store_profiles_to_nc(output_dir, output_file, points,subsets, depth,
                     all_times, all_theta_profiles, all_salt_profiles)









