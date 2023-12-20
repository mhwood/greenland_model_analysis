
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

def plot_comparison_figure(output_file, plot_titles, EN4_timeseries,
                           ECCO_temperature_timeseries, ECCO_locations,
                           L1_temperature_timeseries, L1_locations):

    fig = plt.figure(figsize=(10, 5))

    ############################################################################################################
    # Shelf Entrance Timeseries

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colors = ['purple','olivedrab']
    sources = ['L0','L1','L2']

    ctd_colors = ['brown', 'black','green']
    ctd_sources = ['ARGO','AWI','OMG']

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
            temperature_timeseries[np.array(temperature_timeseries[:,0]).astype(int)==2018,1]=np.nan
        elif source == 'L1':
            temperature_timeseries = L1_temperature_timeseries[L1_locations.index(plot_titles[0])]
        else:
            temperature_timeseries = []

        if len(temperature_timeseries)>0:
            plt.plot(temperature_timeseries[:, 0], temperature_timeseries[:, 1], '-', color=colors[s], label=source)
            # error_polygon = create_error_envelope(temperature_timeseries)
            # error_patch = Polygon(error_polygon, edgecolor='none', facecolor=colors[s], alpha=0.2)
            # plt.gca().add_patch(error_patch)

    plt.ylabel('a) '+plot_titles[0] + '\n Temperature ($^{\circ}$C)')
    # plt.gca().set_xticklabels([])
    plt.title('Temperature Timeseries (200-500m)')

    plt.legend(loc=2, ncol=3)

    plt.gca().set_xlim([1990, 2023])
    plt.grid(linestyle='--', alpha=0.5)

    plt.ylabel('Potential Temperature ($^{\circ}$C)')

    plt.savefig(output_file,bbox_inches='tight')
    plt.close(fig)

project_folder = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund'

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



output_file = os.path.join(project_folder,'Manuscript','Figures','SOM','L1_Regional_Comparison.png')

print('    - Creating the plot')
plot_comparison_figure(output_file, plot_titles,EN4_timeseries,
                       ECCO_temperature_timeseries, ECCO_locations,
                       L1_temperature_timeseries, L1_locations)