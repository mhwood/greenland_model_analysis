
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import scipy.io

def get_EN4_temperature_timeseries(project_folder, polygon_name):

    file_path = os.path.join(project_folder, 'Data', 'In Situ', 'Hadley','Timeseries','Hadley_'+polygon_name+'_Timeseries.nc')

    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    theta = ds.variables['Theta'][:]
    ds.close()
    timeseries = np.column_stack([time, theta])

    timeseries = timeseries[timeseries[:,0]!=0,:]

    return(timeseries)

def get_ECCO_temperature_timeseries(project_folder, polygon_name):

    file_path = os.path.join(project_folder,'Data','Modeling','ECCO','ECCOv5_THETA_'+polygon_name+'_timeseries.nc')

    ds = nc4.Dataset(file_path)
    time = ds.variables['dec_yr'][:]
    theta = ds.variables['THETA_mean'][:]
    ds.close()
    timeseries = np.column_stack([time, theta])

    return(timeseries)

def get_nested_model_temperature_timeseries(project_folder,model_name, polygon_name):

    file_path = os.path.join(project_folder, 'Data', 'Modeling', 'Downscaled',model_name,model_name+'_'+polygon_name+'_Timeseries.nc')

    ds = nc4.Dataset(file_path)
    time = ds.variables['time'][:]
    theta = ds.variables['Theta'][:]
    ds.close()
    timeseries = np.column_stack([time, theta])

    return(timeseries)

def plot_comparison_figure(output_file, EN4_temperature_timeseries,
                           ECCO_temperature_timeseries, L1_temperature_timeseries):

    fig = plt.figure(figsize=(9, 4))

    plt.plot(EN4_temperature_timeseries[:,0],EN4_temperature_timeseries[:,1],color='k',label='EN4 Objective Analysis')
    # plt.plot(ECCO_temperature_timeseries[:, 0], ECCO_temperature_timeseries[:, 1], color='g',label='ECCOv5')
    plt.plot(L1_temperature_timeseries[:, 0], L1_temperature_timeseries[:, 1],'--', color='b',label='L1 Model', alpha=0.5)
    L1_bias = np.mean(L1_temperature_timeseries[:, 1])-np.mean(EN4_temperature_timeseries[:,1])
    plt.plot(L1_temperature_timeseries[:, 0], L1_temperature_timeseries[:, 1]-L1_bias,color='b', label='L1 Model (Bias-Corrected)')
    plt.text(1998,-0.5,'L1 model bias relative to EN4 Objective Analysis: '+'{:.3f}'.format(L1_bias)+'$^{\circ}$C',color='b')

    plt.legend(loc=2,ncol=2)

    plt.ylabel('Temperature ($^{\circ}$C)')
    plt.title('Temperature Timeseries in IPCC Sample Area (200-500m)')

    plt.gca().set_xlim([1990.5, 2023])
    plt.grid(linestyle='--', alpha=0.5)

    plt.savefig(output_file,bbox_inches='tight')
    plt.close(fig)

project_folder = '/Users/michwood/Documents/Research/Projects/Scoresby Sund'

min_depth = 200
max_depth = 500

polygon_name = 'IPCC'

print('    - Reading in the EN4 timeseries')
EN4_temperature_timeseries = get_EN4_temperature_timeseries(project_folder, polygon_name)

print('    - Reading in the ECCO timeseries')
ECCO_temperature_timeseries = get_ECCO_temperature_timeseries(project_folder, polygon_name)

print('    - Reading in the L1 timeseries')
L1_model_name = 'L1_CE_Greenland'
L1_temperature_timeseries = get_nested_model_temperature_timeseries(project_folder, L1_model_name, polygon_name)

output_file = os.path.join(project_folder,'Figures','IPCC_Area_Temperature_Comparison.png')

print('    - Creating the plot')
plot_comparison_figure(output_file, EN4_temperature_timeseries,
                       ECCO_temperature_timeseries, L1_temperature_timeseries)