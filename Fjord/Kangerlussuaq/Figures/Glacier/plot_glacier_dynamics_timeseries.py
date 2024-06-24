

import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.patches import Polygon
from pyproj import Transformer
from datetime import datetime, timedelta

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1, run_test = True):

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
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon

def get_velocity_timeseries(project_folder,glacier):

    file_path = os.path.join(project_folder,'Data','Glacier',
                             glacier+' Sample Point Velocity Timeseries.nc')

    ds = nc4.Dataset(file_path)

    sample_x = ds.sample_point_x
    sample_y = ds.sample_point_y

    point = reproject_polygon(np.array([[sample_x, sample_y]]), 3413, 4326)
    print('    sample point: '+str(point))

    # grp = ds.groups[source]
    time = ds.variables['dec_yr'][:]
    vel = ds.variables['velocity'][:]
    err = ds.variables['error'][:]
    src = ds.variables['source_id'][:]
    timeseries = np.column_stack([time, vel, err, src])
    # source_id_list = ds.variables['source_id'].source_keys
    # src_dict = {2:'MEaSUREs Optical', 4:'ITS_LIVE Annual', 3:'MEaSUREs Quarterly', 1:'MEaSUREs InSAR'}
    # print('    - Note: hard-coded the vel dict')
    source_id_list = ds.variables['source_id'].source_keys.split(', ')
    src_dict = {}
    for src in source_id_list:
        src_dict[int(src.split(':')[0])] = src.split(':')[1]

    ds.close()

    return(timeseries, src_dict, sample_x, sample_y)

def get_elevation_timeseries(project_folder,glacier):

    file_path = os.path.join(project_folder, 'Data', 'Glacier',
                             glacier + ' Sample Point Elevation Timeseries.nc')

    ds = nc4.Dataset(file_path)

    # grp = ds.groups[source]
    time = ds.variables['dec_yr'][:]
    elev = ds.variables['elevation'][:]
    src = ds.variables['source_id'][:]
    timeseries = np.column_stack([time, elev, src])
    source_id_list = ds.variables['source_id'].source_keys.split(', ')
    src_dict = {}
    for src in source_id_list:
        src_dict[int(src.split(':')[0])] = src.split(':')[1]

    ds.close()

    return(timeseries, src_dict)

def get_glacier_retreat_timeseries(project_folder,glacier):

    file_path = os.path.join(project_folder,'Data','Glacier',
                             glacier+' Retreat Timeseries.nc')

    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    retreat = ds.variables['retreat'][:]
    ds.close()

    timeseries = np.column_stack([time,retreat])

    return(timeseries)


def plot_comparison_figure(output_file, retreat_timeseries,
                           velocity_timeseries, velocity_source_dict,
                           elevation_timeseries, elevation_source_dict,
                           for_publication=False):

    min_year = 2015
    max_year = 2022

    dpi=300
    fig = plt.figure(figsize=(8, 8), dpi=dpi)

    ############################################################################################################
    # Velocity Timeseries

    plt.subplot(3,1,1)

    velocity_timeseries = velocity_timeseries[velocity_timeseries[:,0]>2014, :]
    velocity_timeseries = velocity_timeseries[velocity_timeseries[:,1]<6500, :]

    colors = ['steelblue','olivedrab','k','k',]
    symbols = ['.','s','.','^']
    for source_id in [1,2]:
        velocity_timeseries_subset = velocity_timeseries[velocity_timeseries[:,-1]==source_id, :]
        plt.plot(velocity_timeseries_subset[:, 0], velocity_timeseries_subset[:, 1],
                 color=colors[source_id-1], marker=symbols[source_id-1], label=velocity_source_dict[source_id],linewidth=0)
        for i in range(np.shape(velocity_timeseries_subset)[0]):
            # plt.plot([velocity_timeseries[i, 0], velocity_timeseries[i, 0]+1],
            #          [velocity_timeseries[i, 1], velocity_timeseries[i, 1]],'-',linewidth=0.5, color=colors[vt])
            plt.plot([velocity_timeseries_subset[i, 0], velocity_timeseries_subset[i, 0]],
                     [velocity_timeseries_subset[i, 1]-velocity_timeseries_subset[i, 2],
                      velocity_timeseries_subset[i, 1]+velocity_timeseries_subset[i, 2]], '-', linewidth=0.5, color=colors[source_id-1])

    ymin = np.min(velocity_timeseries[:, 1])
    ymax = np.max(velocity_timeseries[:, 1])
    y_range = ymax - ymin
    ymin -= y_range * 0.1
    ymax += y_range * 0.1
    plt.gca().set_ylim([ymin, ymax])

    plt.ylabel('Speed (m/yr)')

    plt.gca().set_xlim([min_year,max_year])
    plt.grid(linestyle='--',alpha=0.5)
    plt.legend(loc=4)
    plt.gca().set_xticklabels([])

    plt.text(min_year+0.1, ymax-y_range*0.1,
             'a) ', ha='left', va='top', fontsize=11)

    ############################################################################################################
    # Elevation Timeseries
    plt.subplot(3, 1, 2)
    color_dict = {1:'black', 4:'steelblue', 5:'black', 7:'olivedrab', 10:'olivedrab'}
    symbol_dict = {1:'^', 4:'s', 5:'s', 7:'^', 10:'.' }
    source_ids = [1,4,7,10]

    elevation_timeseries = elevation_timeseries[elevation_timeseries[:, 0] > 2014, :]

    for elevation_counter in range(len(source_ids)):
        elevation_timeseries_subset = elevation_timeseries[elevation_timeseries[:,-1]==source_ids[elevation_counter], :]
        plt.plot(elevation_timeseries_subset[:, 0], elevation_timeseries_subset[:, 1],
                 color=color_dict[source_ids[elevation_counter]], marker=symbol_dict[source_ids[elevation_counter]], label=elevation_source_dict[source_ids[elevation_counter]],linewidth=0)

    ymin = np.min(elevation_timeseries[:, 1])
    ymax = np.max(elevation_timeseries[:, 1])
    y_range = ymax - ymin
    ymin -= y_range * 0.1
    ymax += y_range * 0.1
    plt.gca().set_ylim([ymin, ymax])

    plt.gca().set_xlim([min_year, max_year])
    plt.grid(linestyle='--', alpha=0.5)
    plt.legend(loc=3, ncol=2)
    plt.ylabel('Ice Elevation (m)')
    plt.gca().set_xticklabels([])

    plt.text(min_year + 0.1, ymax - y_range * 0.1,
             'b) ', ha='left', va='top', fontsize=11)

    ############################################################################################################
    # Retreat Timeseries
    plt.subplot(3, 1, 3)

    set_int = interp1d(retreat_timeseries[:,0], retreat_timeseries[:,1])
    y_shift = set_int(2016)
    retreat_timeseries[:, 1] = (retreat_timeseries[:,1]-y_shift)/(-1000)


    plt.plot(retreat_timeseries[:, 0], retreat_timeseries[:, 1], 'k^', label='TermPicks',markersize=5)

    ymin = np.min(retreat_timeseries[:, 1])
    ymax = np.max(retreat_timeseries[:, 1])
    y_range = ymax - ymin
    ymin -= y_range * 0.1
    ymax += y_range * 0.1
    plt.gca().set_ylim([ymin, ymax])

    plt.gca().set_xlim([min_year, max_year])
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('Ice Front Position\n(km, relative to\n2016)')
    plt.legend(loc=3, ncol=2)

    plt.text(min_year + 0.1, ymax - y_range * 0.1,
             'c) ', ha='left', va='top', fontsize=11)

    if for_publication:
        plt.savefig(output_file[:-3]+'pdf',bbox_inches='tight', dpi=dpi)
    else:
        plt.savefig(output_file, bbox_inches='tight', dpi=dpi)
    plt.close(fig)


glacier = 'Kangerlussuaq'

project_folder = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

longitude = -28.65249
latitude = 71.90093

velocity_timeseries, velocity_source_dict, sample_x, sample_y = get_velocity_timeseries(project_folder,glacier)
print(velocity_source_dict)

elevation_timeseries, elevation_source_dict = get_elevation_timeseries(project_folder,glacier)
print(elevation_source_dict)

retreat_timeseries = get_glacier_retreat_timeseries(project_folder,glacier)

output_file = os.path.join(project_folder,'Figures','Glacier',glacier+' Dynamics Timeseries.png')
plot_comparison_figure(output_file, retreat_timeseries,
                           velocity_timeseries, velocity_source_dict,
                           elevation_timeseries, elevation_source_dict,
                           for_publication=False)















