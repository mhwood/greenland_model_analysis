
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime
from matplotlib.gridspec import GridSpec
from scipy.signal import hann
from scipy.interpolate import interp1d
import shutil

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


def find_ctd_points(target_points, config_dir):
    target_points = np.array(target_points)

    ds = nc4.Dataset(os.path.join(config_dir, 'L2', 'L2_Kanger', 'results', 'dv', 'L2_Kanger_CTD.nc'))
    longitude = ds.variables['longitude'][:]
    latitude = ds.variables['latitude'][:]
    ds.close()

    # print(longitude,latitude)

    ctd_points = []

    for p in range(len(target_points)):
        dist = great_circle_distance(target_points[p, 0], target_points[p, 1], longitude, latitude)
        ctd_points.append(np.argmin(dist))

    return(ctd_points)


def read_mean_profile_from_nc(config_dir, results_dir, var_name, ctd_point):

    ds = nc4.Dataset(os.path.join(config_dir,'L2','L2_Kanger',results_dir,'dv','L2_Kanger_CTD.nc'))
    longitude = ds.variables['longitude'][:]
    latitude = ds.variables['latitude'][:]
    iterations = ds.variables['iterations'][:]
    depth = ds.variables['depth'][:]
    var_grids = ds.variables[var_name][:, :, :]
    ds.close()

    # print(np.max(depth[var_grids[0,:,ctd_point]!=0]))

    longitude = longitude[ctd_point]
    latitude = latitude[ctd_point]
    var_grid = var_grids[:, :, ctd_point]
    var_grid = var_grid[var_grid[:,0]!=0,:]
    profile = np.mean(var_grid, axis = 0)
    profile = np.column_stack([depth,profile])
    profile = profile[profile[:,1]!=0,:]

    return(profile, longitude, latitude)


def read_and_integrate_axctds(project_folder, axctd_folder, target_points, min_depth, max_depth):

    ctd_list_path = os.path.join(project_folder,'Data','Ocean','AXCTDs','AXCTD List at Model Points.csv')
    f = open(ctd_list_path)
    lines = f.read()
    f.close()
    lines = lines.split('\n')
    lines.pop(0)

    output_timeseries = []

    for t in range(len(target_points)):
        print('Searching for CTDs for point '+str(t+1))
        target_point = target_points[t]
        temperature_timeseries = []
        for ll in range(len(lines)):
            line = lines[ll].split(',')
            # print(line)
            lon = float(line[2])
            lat = float(line[3])
            dist = great_circle_distance(lon, lat, target_point[0], target_point[1])
            if dist<10000:
                file_name = 'CTD_'+line[0]+'.nc'
                if file_name not in ['CTD_20171015_104006.nc']:

                    dec_yr = YMD_to_DecYr(int(line[0][:4]), int(line[0][4:6]), int(line[0][6:8]))

                    ds = nc4.Dataset(os.path.join(axctd_folder,line[0][:4],file_name))
                    longitude = ds.longitude
                    latitude = ds.latitude
                    print(' -Reading from ' + file_name, longitude, latitude)
                    depth = ds.variables['depth'][:]
                    temp = ds.variables['potential_temperature'][:]
                    ds.close()

                    mean_temp = np.mean(temp[np.logical_and(depth>=min_depth,depth<=max_depth)])
                    min_temp = np.min(temp[np.logical_and(depth >= min_depth, depth <= max_depth)])
                    max_temp = np.max(temp[np.logical_and(depth >= min_depth, depth <= max_depth)])
                    # plt.plot(temp,depth)
                    # plt.title(str(mean_temp))
                    # plt.show()
                    temperature_timeseries.append([dec_yr, mean_temp, min_temp, max_temp])

        output_timeseries.append(np.array(temperature_timeseries))

    return(output_timeseries)


def plot_profile_comparison(project_dir, var_names, control_profiles, melange_profiles,
                               melange_plume_profiles, plume_profiles, longitudes, latitudes):


    fig = plt.figure(figsize=(8,8))

    gs = GridSpec(2, len(longitudes), top=0.92, bottom=0.07, left = 0.1, right = 0.95)

    max_depth = 650

    # temp
    for ll in range(len(longitudes)):
        ax = fig.add_subplot(gs[0, ll])
        ax.plot(control_profiles[ll][:,1], control_profiles[ll][:, 0], 'black', label='control')
        ax.plot(melange_profiles[ll][:, 1], melange_profiles[ll][:, 0], 'green', label='melange only')
        ax.plot(plume_profiles[ll][:, 1], plume_profiles[ll][:, 0], 'purple', label='plume only')
        ax.plot(melange_plume_profiles[ll][:, 1], melange_plume_profiles[ll][:, 0], 'darkorange', label='melange + plume')

        ax.set_ylim([650,0])

        min_temp = np.min([np.min(control_profiles[ll][:, 1]),
                           np.min(melange_profiles[ll][:, 1]),
                           np.min(plume_profiles[ll][:, 1]),
                           np.min(melange_plume_profiles[ll][:, 1])])
        max_temp = np.max([np.max(control_profiles[ll][:, 1]),
                           np.max(melange_profiles[ll][:, 1]),
                           np.max(plume_profiles[ll][:, 1]),
                           np.max(melange_plume_profiles[ll][:, 1])])
        # ax.plot(axctd_profile[ll][:,0],axctd_profile[ll][:,1],'k.',markersize=10, label='CTD')
        # ax.plot([axctd_profile[ll][:, 0],axctd_profile[ll][:, 0]],
        #         [axctd_profile[ll][:, 2],axctd_profile[ll][:, 3]], 'k-')
        min_temp = 1.25
        max_temp = 3
        ax.set_xlim([min_temp, max_temp])
        ax.grid(linewidth=0.5, alpha = 0.5, linestyle='--')

        if ll==0:
            plt.title('Continental Shelf\n('+'{:.2f}'.format(-1*longitudes[ll])+'$^{\circ}$W, '+
                      '{:.2f}'.format(latitudes[ll])+'$^{\circ}$N)')
            ax.text(min_temp+0.09, max_depth-30, 'a)', ha='left', va='bottom')#,
                     #  bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
            ax.set_ylabel('Depth (m)')


        if ll==1:
            ax.set_xlabel('Temperature ($^{\circ}C$)')
            ax.set_yticklabels([])
            ax.set_title('Fjord Mouth\n('+'{:.2f}'.format(-1*longitudes[ll])+'$^{\circ}$W, '+
                          '{:.2f}'.format(latitudes[ll])+'$^{\circ}$N)')
            ax.text(min_temp+0.09, max_depth-30, 'b)', ha='left', va='bottom')#,
                  #  bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

        if ll==2:
            ax.set_title('Near Glacier\n('+'{:.2f}'.format(-1*longitudes[ll])+'$^{\circ}$W, '+
                          '{:.2f}'.format(latitudes[ll])+'$^{\circ}$N)')
            ax.set_yticklabels([])
            ax.text(min_temp+0.09, max_depth-30, 'c)', ha='left', va='bottom')#,
                 #   bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    for ll in range(len(longitudes)):
        ax = fig.add_subplot(gs[1, ll])

        ax.plot(control_profiles[ll+len(longitudes)][:,1], control_profiles[ll+len(longitudes)][:, 0], 'black', label='control')
        ax.plot(melange_profiles[ll+len(longitudes)][:, 1], melange_profiles[ll+len(longitudes)][:, 0], 'green', label='melange')
        ax.plot(plume_profiles[ll+len(longitudes)][:, 1], plume_profiles[ll+len(longitudes)][:, 0], 'purple', label='plume')
        ax.plot(melange_plume_profiles[ll+len(longitudes)][:, 1], melange_plume_profiles[ll+len(longitudes)][:, 0], 'darkorange', label='melange\n+ plume')

        ax.set_ylim([max_depth,0])

        min_salt = np.min([np.min(control_profiles[ll+len(longitudes)][:, 1]),
                           np.min(melange_profiles[ll+len(longitudes)][:, 1]),
                           np.min(plume_profiles[ll+len(longitudes)][:, 1]),
                           np.min(melange_plume_profiles[ll+len(longitudes)][:, 1])])
        max_salt = np.max([np.max(control_profiles[ll+len(longitudes)][:, 1]),
                           np.max(melange_profiles[ll+len(longitudes)][:, 1]),
                           np.max(plume_profiles[ll+len(longitudes)][:, 1]),
                           np.max(melange_plume_profiles[ll+len(longitudes)][:, 1])])
        min_salt = 33.5
        max_salt = 35.5
        # ax.plot(axctd_profile[ll][:,0],axctd_profile[ll][:,1],'k.',markersize=10, label='CTD')
        # ax.plot([axctd_profile[ll][:, 0],axctd_profile[ll][:, 0]],
        #         [axctd_profile[ll][:, 2],axctd_profile[ll][:, 3]], 'k-')

        ax.set_xlim([min_salt, max_salt])
        ax.grid(linewidth=0.5, alpha = 0.5, linestyle='--')

        if ll==0:
            ax.text(min_salt+0.09, max_depth-30, 'd)', ha='left', va='bottom')#,
                       #bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
            ax.set_ylabel('Depth (m)')
            ax.legend(loc=6)


        if ll==1:
            ax.set_xlabel('Salinity (psu)')
            ax.set_yticklabels([])
            ax.text(min_salt+0.09, max_depth-30, 'e)', ha='left', va='bottom')#,
                  #  bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

        if ll==2:
            ax.set_yticklabels([])
            ax.text(min_salt+0.09, max_depth-30, 'f)', ha='left', va='bottom')#,
                   # bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    output_file = os.path.join(project_dir, 'Figures', 'Ocean', 'Modeling', 'Temperature Profile Comparison.png')
    plt.savefig(output_file)
    plt.close(fig)


project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

axctd_folder = '/Users/mike/Documents/Research/Data Repository/Greenland/Ocean Properties/OMG_CTDs/Processed'

config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
                 'configurations/downscaled_greenland'


var_names = ['THETA','SALT']

# shelf, fjord mouth, near-glacier
target_points = [[-31.0496,67.2902],[-31.8221, 68.0413],[-32.5201, 68.4994]]

# get the axctds
#axctd_profile = read_and_integrate_axctds(project_dir, axctd_folder, target_points, min_depth, max_depth)

ctd_points = find_ctd_points(target_points, config_dir)
# print(ctd_points)

control_profiles = []
melange_profiles = []
melange_plume_profiles = []
plume_profiles = []

longitudes = []
latitudes = []

for var_name in var_names:
    for ctd_point in ctd_points:

        # read in the control profile
        profile, longitude, latitude = \
            read_mean_profile_from_nc(config_dir, 'results', var_name, ctd_point)
        control_profiles.append(profile)

        # add the lat/lon to the lists
        if var_name=='THETA':
            longitudes.append(longitude)
            latitudes.append(latitude)

        profile, _, _ = \
            read_mean_profile_from_nc(config_dir, 'results_melange', var_name, ctd_point)
        melange_profiles.append(profile)

        profile, _, _ = \
            read_mean_profile_from_nc(config_dir, 'results_melange_plume', var_name, ctd_point)
        melange_plume_profiles.append(profile)

        profile, _, _ = \
            read_mean_profile_from_nc(config_dir, 'results_plume', var_name, ctd_point)
        plume_profiles.append(profile)

# plt.plot(iterations,var_grid)
# plt.show()

print(len(control_profiles))

plot_profile_comparison(project_dir, var_names, control_profiles, melange_profiles,
                           melange_plume_profiles, plume_profiles, longitudes, latitudes)


shutil.copyfile(os.path.join(project_dir, 'Figures', 'Ocean', 'Modeling', 'Temperature Profile Comparison.png'),
                 os.path.join(project_dir, 'Manuscript','Draft 1', 'Figure_4.png'))