
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import scipy.io

def get_hadley_objective_analysis_timeseries(project_folder):
    file_path = os.path.join(project_folder,'Data','In Situ','Hadley','Timeseries','Hadley_IPCC_Timeseries.nc')
    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    temp = ds.variables['Theta'][:]
    ds.close()
    timeseries = np.column_stack([time,temp])
    timeseries [timeseries[:,0]!=0,:]

    return(timeseries)

def get_regional_ECCO_temperature_timeseries(project_folder):
    file_path = os.path.join(project_folder,'Data','Modeling','ECCO','ECCOv5_THETA_IPCC_timeseries.nc')
    ds = nc4.Dataset(file_path)
    time = ds.variables['dec_yr'][:]
    temp = ds.variables['THETA_mean'][:]
    ds.close()
    timeseries = np.column_stack([time,temp])
    return(timeseries)

def get_ECCO_temperature_timeseries(project_folder,ECCO_locations):

    all_timeseries = []

    file_path = os.path.join(project_folder,'Data','Modeling','ECCO','ECCOv5_Scoresby_Sund_Timeseries.nc')

    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]

    for location_name in ECCO_locations:
        grp_name = '_'.join(location_name.lower().split())
        if grp_name in list(ds.groups.keys()):
            grp = ds.groups['_'.join(location_name.lower().split())]
            theta = grp.variables['THETA'][:]
            all_timeseries.append(np.column_stack([time,theta]))
        else:
            all_timeseries.append([])

    ds.close()

    return(all_timeseries)

def get_nested_model_temperature_timeseries(project_folder,model_name,locations):

    retrieval_locations = list(locations)
    all_timeseries = []

    # if we have the regional timeseries, then get it
    if model_name=='L1_CE_Greenland' and retrieval_locations[0]=='Regional_Ocean':
        file_path = os.path.join(project_folder, 'Data', 'Modeling', 'Downscaled',model_name, 'L1_CE_Greenland_IPCC_Timeseries.nc')
        ds = nc4.Dataset(file_path)
        time = ds.variables['time'][:]
        temp = ds.variables['Theta'][:]
        ds.close()
        timeseries = np.column_stack([time, temp])
        all_timeseries.append(timeseries)
        retrieval_locations.pop(0)

    file_path = os.path.join(project_folder, 'Data', 'Modeling', 'Downscaled', model_name,
                             model_name + '_dv_CTD_Timeseries.nc')
    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    for location_name in retrieval_locations:
        grp_name = '_'.join(location_name.lower().split())
        if grp_name in list(ds.groups.keys()):
            grp = ds.groups[grp_name]
            theta = grp.variables['THETA'][:]
            all_timeseries.append(np.column_stack([time, theta]))
        else:
            all_timeseries.append([])
    ds.close()

    return(all_timeseries)

def get_in_situ_temperature_data(project_folder,source,locations):

    all_timeseries = []

    for location_name in locations:

        file_path = os.path.join(project_folder,'Data','In Situ',source,'Data',source+'_'+location_name+'_Temperature_Timeseries.nc')

        if os.path.exists(file_path):

            ds = nc4.Dataset(file_path)
            time = ds.variables['time'][:]
            theta = ds.variables['temperature_mean'][:]
            theta_std = ds.variables['temperature_std'][:]
            ds.close()

            timeseries = np.column_stack([time,theta,theta_std])

        else:
            timeseries = []

        all_timeseries.append(timeseries)

    return(all_timeseries)

def create_error_envelope(temperature_timeseries):

    upper_error = np.column_stack(
        [temperature_timeseries[:, 0], temperature_timeseries[:, 1] + temperature_timeseries[:, 2]])
    lower_error = np.column_stack(
        [temperature_timeseries[:, 0], temperature_timeseries[:, 1] - temperature_timeseries[:, 2]])
    error_polygon = np.vstack([upper_error, np.flipud(lower_error)])
    return(error_polygon)

def get_slater_TF_estimate(project_folder):

    file_path = os.path.join(project_folder,'Data','Comparison','slater_data.mat')
    # print(file_path)
    #
    # f = h5py.File(file_path, 'r')
    # print(f.keys())
    # data = f.get('data/variable1')
    # data = np.array(data)  # For converting to a NumPy array

    mat = scipy.io.loadmat(file_path)

    #DJG = #27
    djg_data = mat['glaciers'][0,26]

    djg_EN4_timeseries = djg_data[10]
    # print(djg_data[10][0][0])

    djg_Racmo_Qsg_timeseries = np.column_stack([djg_data[9][0][0][0].ravel(),djg_data[9][0][0][1].ravel()])

    djg_EN4_timeseries = np.column_stack([djg_data[10][0][0][1].ravel(), djg_data[10][0][0][0].ravel()])

    djg_melt_timeseries = np.column_stack([djg_data[13][0][0][0].ravel(), djg_data[13][0][0][1].ravel()])

    # plt.subplot(3,1,1)
    # plt.plot(djg_Racmo_Qsg_timeseries[:,0],djg_Racmo_Qsg_timeseries[:,1])
    # plt.title('RACMO')
    #
    # plt.subplot(3, 1, 2)
    # plt.plot(djg_EN4_timeseries[:, 0], djg_EN4_timeseries[:, 1])
    # plt.title('EN4 TF')
    #
    # plt.subplot(3, 1, 3)
    # plt.plot(djg_melt_timeseries[:, 0], djg_melt_timeseries[:, 1])
    # plt.title('EN4 Melt')
    # plt.show()

    return(djg_EN4_timeseries)

def plot_comparison_figure(output_file, plot_titles, EN4_timeseries,
                           ECCO_temperature_timeseries, ECCO_locations,
                           L1_temperature_timeseries, L1_locations,
                           L2_temperature_timeseries, L2_locations,
                           hadley_temperature_timeseries, hadley_locations,
                           awi_temperature_timeseries, awi_locations,
                           omg_temperature_timeseries, omg_locations):

    fig = plt.figure(figsize=(8, 11))

    ############################################################################################################
    # Shelf Entrance Timeseries

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    sources = ['L0','L1','L2']

    ctd_colors = ['brown', 'black','green']
    ctd_sources = ['ARGO','AWI','OMG']

    plt.subplot(4, 1, 1)
    hadley_years = np.arange(1991,2023).tolist()
    for year in hadley_years:
        year_indices = np.logical_and(EN4_timeseries[:,0]>=year,EN4_timeseries[:,0]<year+1)
        year_mean = np.mean(EN4_timeseries[year_indices,1])
        if year==hadley_years[0]:
            plt.plot([year,year+1],[year_mean,year_mean],color='k',linewidth=1, label='EN4 Objective Analysis')
        else:
            plt.plot([year, year + 1], [year_mean, year_mean], color='k', linewidth=1)
        plt.plot(year + 0.5, year_mean,'.', color='k', markersize=5)
    # for i in range(len(EN4_TF_data)):
    #     plt.plot([EN4_TF_data[i, 0]-0.5,EN4_TF_data[i, 0]+0.5],
    #              [EN4_TF_data[i, 1] - 1.9, EN4_TF_data[i, 1] - 1.9], '-', color='k', linewidth=0.5)
    for s in range(len(sources)):
        source = sources[s]
        if source == 'L0':
            temperature_timeseries = ECCO_temperature_timeseries[ECCO_locations.index(plot_titles[0])]
        elif source == 'L1':
            temperature_timeseries = L1_temperature_timeseries[L1_locations.index(plot_titles[0])]
        else:
            temperature_timeseries = []

        if len(temperature_timeseries)>0:
            plt.plot(temperature_timeseries[:, 0], temperature_timeseries[:, 1], '-', color=colors[s], label=source,linewidth=1)
            # error_polygon = create_error_envelope(temperature_timeseries)
            # error_patch = Polygon(error_polygon, edgecolor='none', facecolor=colors[s], alpha=0.2)
            # plt.gca().add_patch(error_patch)

    plt.ylabel('a) '+plot_titles[0] + '\n Temperature ($^{\circ}$C)')
    plt.gca().set_xticklabels([])
    plt.title('Temperature Timeseries (200-400m)')

    plt.legend(loc=2, ncol=3)

    plt.gca().set_xlim([1990, 2023])
    plt.grid(linestyle='--', alpha=0.5)

    ############################################################################################################
    # Sound Entrance

    plt.subplot(4,1,2)
    for s in range(len(sources)):
        source = sources[s]
        # if source == 'ECCO':
        #     temperature_timeseries = ECCO_temperature_timeseries[ECCO_locations.index(plot_titles[1])]
        if source == 'L1':
            temperature_timeseries = L1_temperature_timeseries[L1_locations.index(plot_titles[1])]
        elif source == 'L2':
            temperature_timeseries = L2_temperature_timeseries[L2_locations.index(plot_titles[1])]
        else:
            temperature_timeseries = []

        if len(temperature_timeseries)>0:
            plt.plot(temperature_timeseries[:, 0], temperature_timeseries[:, 1], '-', color=colors[s], label=source,linewidth=1)
            # error_polygon = create_error_envelope(temperature_timeseries)
            # error_patch = Polygon(error_polygon, edgecolor='none', facecolor=colors[s], alpha=0.2)
            # plt.gca().add_patch(error_patch)

    for s in range(len(ctd_sources)):
        source = ctd_sources[s]
        if source == 'ARGO':
            temperature_timeseries = hadley_temperature_timeseries[hadley_locations.index(plot_titles[1])]
        elif source == 'AWI':
            temperature_timeseries = awi_temperature_timeseries[awi_locations.index(plot_titles[1])]
        elif source == 'OMG':
            temperature_timeseries = omg_temperature_timeseries[omg_locations.index(plot_titles[1])]
        else:
            temperature_timeseries = []


        if len(temperature_timeseries) > 0:
            plt.plot(temperature_timeseries[:, 0], temperature_timeseries[:, 1], '.', color='k', label=source)
            # for i in range(np.shape(temperature_timeseries)[0]):
                # plt.plot([temperature_timeseries[i, 0], temperature_timeseries[i, 0]],
                #          [temperature_timeseries[i, 1] - temperature_timeseries[i, 2],
                #           temperature_timeseries[i, 1] + temperature_timeseries[i, 2]],
                #          '-', color=ctd_colors[s], linewidth=0.8)

    plt.legend(loc=2,ncol=3)

    plt.ylabel('b) '+plot_titles[1]+'\n Temperature ($^{\circ}$C)')
    plt.gca().set_xticklabels([])
    # plt.title('Temperature Timeseries at Sample Sites Through Scoresby Sund')

    # plt.legend(loc=2,ncol=3)

    plt.gca().set_xlim([1990, 2023])
    plt.gca().set_ylim([-1,1.9])
    plt.grid(linestyle='--',alpha=0.5)

    ############################################################################################################
    # Mid-Sound
    plt.subplot(4, 1, 3)

    for s in range(len(sources)):
        source = sources[s]
        # if source == 'ECCO':
        #     temperature_timeseries = ECCO_temperature_timeseries[ECCO_locations.index(plot_titles[1])]
        if source == 'L1':
            temperature_timeseries = L1_temperature_timeseries[L1_locations.index(plot_titles[2])]
        elif source == 'L2':
            temperature_timeseries = L2_temperature_timeseries[L1_locations.index(plot_titles[2])]
        else:
            temperature_timeseries = []

        if len(temperature_timeseries) > 0:
            plt.plot(temperature_timeseries[:, 0], temperature_timeseries[:, 1], '-', color=colors[s], label=source,
                     linewidth=1)
            # error_polygon = create_error_envelope(temperature_timeseries)
            # error_patch = Polygon(error_polygon, edgecolor='none', facecolor=colors[s], alpha=0.2)
            # plt.gca().add_patch(error_patch)

    for s in range(len(ctd_sources)):
        source = ctd_sources[s]
        if source == 'AWI':
            temperature_timeseries = awi_temperature_timeseries[awi_locations.index(plot_titles[2])]
        elif source == 'OMG':
            temperature_timeseries = omg_temperature_timeseries[omg_locations.index(plot_titles[2])]
        else:
            temperature_timeseries = []

        if len(temperature_timeseries) > 0:
            plt.plot(temperature_timeseries[:, 0], temperature_timeseries[:, 1], '.', color=ctd_colors[s], label=source)
            # for i in range(np.shape(temperature_timeseries)[0]):
            #     plt.plot([temperature_timeseries[i, 0], temperature_timeseries[i, 0]],
            #              [temperature_timeseries[i, 1] - temperature_timeseries[i, 2],
            #               temperature_timeseries[i, 1] + temperature_timeseries[i, 2]],
            #              '-', color=ctd_colors[s], linewidth=0.8)

    plt.text(1991, 1, 'Not resolved by ECCO', color='blue')

    plt.gca().set_xlim([1990, 2023])
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('c) '+plot_titles[2] + '\n Temperature ($^{\circ}$C)')
    plt.gca().set_xticklabels([])

    ############################################################################################################
    # Near-Glacier Timeseries

    plt.subplot(4,1,4)

    for s in range(len(sources)):
        source = sources[s]
        if source == 'L1':
            temperature_timeseries = L1_temperature_timeseries[L1_locations.index(plot_titles[3])]
        elif source == 'L2':
            temperature_timeseries = L2_temperature_timeseries[L2_locations.index(plot_titles[3])]
        else:
            temperature_timeseries = []

        if len(temperature_timeseries) > 0:
            plt.plot(temperature_timeseries[:, 0], temperature_timeseries[:, 1], '-', color=colors[s], label=source,linewidth=0.5)
            # error_polygon = create_error_envelope(temperature_timeseries)
            # error_patch = Polygon(error_polygon, edgecolor='none', facecolor=colors[s], alpha=0.2)
            # plt.gca().add_patch(error_patch)

    for s in range(len(ctd_sources)):
        source = ctd_sources[s]
        if source == 'AWI':
            temperature_timeseries = awi_temperature_timeseries[awi_locations.index(plot_titles[3])]
        elif source == 'OMG':
            temperature_timeseries = omg_temperature_timeseries[omg_locations.index(plot_titles[3])]
        else:
            temperature_timeseries = []

        if len(temperature_timeseries) > 0:
            plt.plot(temperature_timeseries[:, 0], temperature_timeseries[:, 1], '.', color=ctd_colors[s], label=source)
            # for i in range(np.shape(temperature_timeseries)[0]):
            #     plt.plot([temperature_timeseries[i, 0], temperature_timeseries[i, 0]],
            #              [temperature_timeseries[i, 1] - temperature_timeseries[i, 2],
            #               temperature_timeseries[i, 1] + temperature_timeseries[i, 2]],
            #              '-', color=ctd_colors[s], linewidth=0.8)

    # plt.text(1991, 1, 'Not resolved by ECCO', color='blue')

    plt.gca().set_xlim([1990, 2023])
    plt.grid(linestyle='--', alpha=0.5)

    plt.ylabel('d) '+plot_titles[3] + '\n Temperature ($^{\circ}$C)')

    plt.savefig(output_file,bbox_inches='tight')
    plt.close(fig)

project_folder = '/Users/michwood/Documents/Research/Projects/Scoresby Sund'

plot_titles = ['Regional Ocean', 'Sound Entrance', 'Mid-Sound', 'Near-DJG']
output_names =['Regional_Ocean', 'Sound_Entrance', 'Mid_Sound', 'Near_Glacier']

min_depth = 200
max_depth = 400

print('    - Reading the regional Hadley timeseries')
EN4_timeseries = get_hadley_objective_analysis_timeseries(project_folder)

print('    - Reading in the ECCO timeseries')
ECCO_locations = plot_titles[:1]
ECCO_sample_locations = output_names[:1]
ECCO_temperature_timeseries = []
ECCO_temperature_timeseries.append(get_regional_ECCO_temperature_timeseries(project_folder))

print('    - Reading in the L1 timeseries')
L1_model_name = 'L1_CE_Greenland'
L1_locations = plot_titles
L1_temperature_timeseries = get_nested_model_temperature_timeseries(project_folder, L1_model_name, output_names)

print('    - Reading in the L2 timeseries')
L3_model_name = 'L3_Scoresby_Sund'
L2_locations = plot_titles
L2_temperature_timeseries = get_nested_model_temperature_timeseries(project_folder, L3_model_name, output_names)

print('    - Reading in Hadley timeseries')
hadley_locations = plot_titles
hadley_temperature_timeseries = get_in_situ_temperature_data(project_folder,'Hadley',output_names)

print('    - Reading in AWI timeseries')
awi_locations = plot_titles
awi_temperature_timeseries = get_in_situ_temperature_data(project_folder,'AWI',output_names)

print('    - Reading in OMG timeseries')
omg_locations = plot_titles
omg_temperature_timeseries = get_in_situ_temperature_data(project_folder,'OMG',output_names)


output_file = os.path.join(project_folder,'Manuscript','Figures','Figure 2','Scoresby_Sund_Temperature_Comparison.png')

print('    - Creating the plot')
plot_comparison_figure(output_file, plot_titles,EN4_timeseries,
                       ECCO_temperature_timeseries, ECCO_locations,
                       L1_temperature_timeseries, L1_locations,
                       L2_temperature_timeseries, L2_locations,
                       hadley_temperature_timeseries, hadley_locations,
                       awi_temperature_timeseries, awi_locations,
                       omg_temperature_timeseries, omg_locations)