
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


def read_and_integrate_timeseries_from_nc(config_dir, results_dir, var_name, ctd_point, min_depth, max_depth):

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

    integrated_vargrid = np.zeros((np.shape(var_grid)[0],))
    interp_depths = np.arange(min_depth,max_depth+1)
    for v in range(len(integrated_vargrid)):
        set_int = interp1d(depth,var_grid[v,:])
        integrated_vargrid[v] = np.mean(set_int(interp_depths))

    # compute the dec yrs
    dec_yrs = np.zeros_like(iterations)
    for d in range(len(iterations)):
        dec_yrs[d] = iter_number_to_dec_yr(iterations[d])

    timeseries = np.column_stack([dec_yrs, integrated_vargrid])

    timeseries = timeseries[np.logical_and(timeseries[:,1]!=0,~np.isnan(timeseries[:,1])),:]

    win = hann(24)
    win = win / np.sum(win)
    timeseries[:,1] = np.convolve(timeseries[:,1], win, mode='same')

    timeseries = timeseries[:-7,:]

    return(timeseries, longitude, latitude)


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


def plot_timeseries_comparison(project_dir, control_timeseries, melange_timeseries,
                               melange_plume_timeseries, plume_timeseries, longitudes, latitudes,
                               axctd_timeseries):


    fig = plt.figure(figsize=(8,8))

    gs = GridSpec(len(longitudes), 1, top=0.92, bottom=0.05, left = 0.1, right = 0.95)

    for ll in range(len(longitudes)):
        ax = fig.add_subplot(gs[ll, 0])
        ax.plot(control_timeseries[ll][:,0], control_timeseries[ll][:, 1], 'black', label='control')
        ax.plot(melange_timeseries[ll][:, 0], melange_timeseries[ll][:, 1], 'green', label='melange only')
        ax.plot(plume_timeseries[ll][:, 0], plume_timeseries[ll][:, 1], 'purple', label='plume only')
        ax.plot(melange_plume_timeseries[ll][:, 0], melange_plume_timeseries[ll][:, 1], 'darkorange', label='melange + plume')

        min_temp = np.min([np.min(control_timeseries[ll][:, 1]),
                           np.min(melange_timeseries[ll][:, 1]),
                           np.min(plume_timeseries[ll][:, 1]),
                           np.min(melange_plume_timeseries[ll][:, 1])])
        max_temp = np.max([np.max(control_timeseries[ll][:, 1]),
                           np.max(melange_timeseries[ll][:, 1]),
                           np.max(plume_timeseries[ll][:, 1]),
                           np.max(melange_plume_timeseries[ll][:, 1])])
        # ax.plot(axctd_timeseries[ll][:,0],axctd_timeseries[ll][:,1],'k.',markersize=10, label='CTD')
        # ax.plot([axctd_timeseries[ll][:, 0],axctd_timeseries[ll][:, 0]],
        #         [axctd_timeseries[ll][:, 2],axctd_timeseries[ll][:, 3]], 'k-')

        ax.set_xlim([2015,2022])
        ax.grid(linewidth=0.5, alpha = 0.5, linestyle='--')

        if ll==0:
            plt.title('Mean Temperature Timeseries (200-500 m)\n'
                      'Location: Continental Shelf ('+'{:.2f}'.format(-1*longitudes[ll])+'$^{\circ}$W, '+
                      '{:.2f}'.format(latitudes[ll])+'$^{\circ}$N)')
            plt.legend(loc=3, ncol=2)
            ax.set_xticklabels([])
            ax.text(2015.1,max_temp, 'a)', ha='left', va='top')#,
                   #    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

        if ll==1:
            ax.set_ylabel('Temperature ($^{\circ}C$)')
            ax.set_xticklabels([])
            ax.set_title('Location: Fjord Mouth ('+'{:.2f}'.format(-1*longitudes[ll])+'$^{\circ}$W, '+
                          '{:.2f}'.format(latitudes[ll])+'$^{\circ}$N)')
            ax.text(2015.1, max_temp, 'b)', ha='left', va='top')#,
                #    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

        if ll==2:
            ax.set_title('Location: Near Kangerlussuaq Glacier ('+'{:.2f}'.format(-1*longitudes[ll])+'$^{\circ}$W, '+
                          '{:.2f}'.format(latitudes[ll])+'$^{\circ}$N)')
            ax.text(2015.1, max_temp, 'c)', ha='left', va='top')#,
                    #bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    output_file = os.path.join(project_dir, 'Figures', 'Ocean', 'Modeling', 'Temperature Timeseries Comparison.png')
    plt.savefig(output_file)
    plt.close(fig)


project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

axctd_folder = '/Users/mike/Documents/Research/Data Repository/Greenland/Ocean Properties/OMG_CTDs/Processed'

config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
                 'configurations/downscaled_greenland'


var_name = 'THETA'
min_depth = 200
max_depth = 500

# shelf, fjord mouth, near-glacier
target_points = [[-31.0496,67.2902],[-31.8221, 68.0413],[-32.5201, 68.4994]]

# get the axctds
axctd_timeseries = read_and_integrate_axctds(project_dir, axctd_folder, target_points, min_depth, max_depth)

ctd_points = find_ctd_points(target_points, config_dir)
# print(ctd_points)

control_timeseries = []
melange_timeseries = []
melange_plume_timeseries = []
plume_timeseries = []

longitudes = []
latitudes = []

for ctd_point in ctd_points:

    # read in the control timeseries
    timeseries, longitude, latitude = \
        read_and_integrate_timeseries_from_nc(config_dir, 'results', var_name, ctd_point, min_depth, max_depth)
    control_timeseries.append(timeseries)

    # add the lat/lon to the lists
    longitudes.append(longitude)
    latitudes.append(latitude)

    timeseries, _, _ = \
        read_and_integrate_timeseries_from_nc(config_dir, 'results_melange', var_name, ctd_point, min_depth, max_depth)
    melange_timeseries.append(timeseries)

    timeseries, _, _ = \
        read_and_integrate_timeseries_from_nc(config_dir, 'results_melange_plume', var_name, ctd_point, min_depth, max_depth)
    melange_plume_timeseries.append(timeseries)

    timeseries, _, _ = \
        read_and_integrate_timeseries_from_nc(config_dir, 'results_plume', var_name, ctd_point, min_depth, max_depth)
    plume_timeseries.append(timeseries)

# plt.plot(iterations,var_grid)
# plt.show()

plot_timeseries_comparison(project_dir, control_timeseries, melange_timeseries,
                           melange_plume_timeseries, plume_timeseries, longitudes, latitudes,
                           axctd_timeseries)


shutil.copyfile(os.path.join(project_dir, 'Figures', 'Ocean', 'Modeling', 'Temperature Timeseries Comparison.png'),
                 os.path.join(project_dir, 'Manuscript','Draft 1', 'Figure_3.png'))