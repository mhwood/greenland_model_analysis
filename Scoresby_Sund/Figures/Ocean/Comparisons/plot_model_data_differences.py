


import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4

def read_model_data_temperatures(project_folder,year, min_depth, max_depth):

    file_path = os.path.join(project_folder,'Data','In Situ','Hadley','Model_Comparison','L1_CE_Greenland',
                               'L1_CE_Model_Data_Differences_'+str(year)+'_'+str(min_depth)+'m_'+str(max_depth)+'m.nc')

    ds = nc4.Dataset(file_path)
    model_theta = ds.variables['model_theta'][:]
    obs_theta = ds.variables['obs_theta'][:]
    ds.close()

    return(model_theta,obs_theta)

def pearson_correlation(x,y):
    meanX = np.mean(x)
    meanY = np.mean(y)
    numerator = np.sum((x-meanX)*(y-meanY))
    denominator1 = np.sum((x-meanX)*(x-meanX))
    denominator2 = np.sum((y - meanY) * (y - meanY))
    r = numerator/(denominator1**0.5*denominator2**0.5)
    return(r)

def plot_scatter_comparison(output_file,model_theta,obs_theta,min_theta,max_theta):

    r = pearson_correlation(model_theta,obs_theta)

    fig = plt.figure(figsize=(8,8))

    plt.plot(obs_theta,model_theta,'g.')
    plt.plot([min_theta,max_theta],[min_theta,max_theta],'k-')

    plt.gca().set_xlim([min_theta, max_theta])
    plt.gca().set_ylim([min_theta, max_theta])

    plt.text(min_theta+0.5,max_theta-0.5,'r = '+'{:0.4f}'.format(r),va = 'top', ha='left')
    plt.text(min_theta + 0.5, max_theta - 1, 'mean model-obs difference = ' + '{:0.4f}'.format(np.mean(model_theta-obs_theta)), va='top', ha='left')

    plt.xlabel('Observations (Hadley Centre EN4.2.2)')
    plt.ylabel('Model (L1_CE_Greenland)')

    plt.title('Model-Data Differences')

    plt.savefig(output_file)
    plt.close(fig)

def plot_annual_difference_comparison(project_folder, year, model_theta, obs_theta, plot_type = 'scatter'):

    min_theta = -1.9
    max_theta = 7
    theta_step = 0.5

    output_file = os.path.join(project_folder,'Figures','Ocean','L1_CE_Greenland','Model Data Comparison',
                               'L1_CE_Model_Data_Differences_'+str(year)+'_'+str(min_depth)+'m_'+str(max_depth)+'.png')

    if plot_type == 'scatter':
        plot_scatter_comparison(output_file, model_theta, obs_theta, min_theta, max_theta)

def plot_total_difference_comparison(project_folder, all_points, plot_type = 'scatter'):

    model_theta = all_points[:,0]
    obs_theta = all_points[:,1]

    min_theta = -1.9
    max_theta = 7
    theta_step = 0.5

    output_file = os.path.join(project_folder,'Figures','Ocean','L1_CE_Greenland','Model Data Comparison',
                               'L1_CE_Model_Data_Differences_All'+'_'+str(min_depth)+'m_'+str(max_depth)+'.png')

    if plot_type == 'scatter':
        plot_scatter_comparison(output_file, model_theta, obs_theta, min_theta, max_theta)

project_folder = '/Users/michwood/Documents/Research/Projects/Scoresby Sund'

min_depth = 200
max_depth = 500

all_points_started = False

for year in range(1992,2022):

    model_theta, obs_theta = read_model_data_temperatures(project_folder,year, min_depth, max_depth)

    plot_annual_difference_comparison(project_folder, year, model_theta, obs_theta)

    if not all_points_started:
        all_points_started = True
        all_points = np.column_stack([model_theta, obs_theta])
    else:
        all_points = np.vstack([all_points,np.column_stack([model_theta, obs_theta])])

plot_total_difference_comparison(project_folder, all_points, plot_type = 'scatter')

