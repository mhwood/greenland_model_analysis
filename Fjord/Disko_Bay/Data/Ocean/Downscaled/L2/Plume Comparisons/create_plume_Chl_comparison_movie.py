
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import netCDF4 as nc4



def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    drF = ds.variables['drF'][:]
    ds.close()
    Z_bottom = np.cumsum(drF)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2
    return(XC, YC, Z, Depth)


def read_Chl_from_file(config_dir, model_name, month, experiment_name):

    if experiment_name=='control':
        results_dir = os.path.join(config_dir,'L2',model_name,'results','BGC_daily_Chl')
    else:
        results_dir = os.path.join(config_dir,'L2',model_name,'results_iceplume','BGC_daily_Chl')

    file_name_to_use = ''
    for file_name in os.listdir(results_dir):
        if file_name[0]!='.' and file_name[-3:]=='.nc':
            y = int(file_name.split('.')[1][:4])
            m = int(file_name.split('.')[1][4:6])
            if month==m and y==1993:
                file_name_to_use = file_name

    ds = nc4.Dataset(os.path.join(results_dir, file_name_to_use))
    chl = ds.variables['Chl01'][:, :, :]
    chl += ds.variables['Chl02'][:, :, :]
    chl += ds.variables['Chl03'][:, :, :]
    chl += ds.variables['Chl04'][:, :, :]
    chl += ds.variables['Chl05'][:, :, :]
    iterations = ds.variables['iterations'][:]
    ds.close()

    return(iterations, chl)



def plot_Chl_panel(output_folder, month, day, XC, YC, Depth,
                   iterations_control, chl_control, iterations_iceplume, chl_iceplume):

    fig = plt.figure(figsize=(14, 5))
    plt.style.use('dark_background')

    vmin = 0
    vmax = 5

    plot_width = 8

    gs = GridSpec(1, 2*plot_width+3,  left=0.05, right=0.95, bottom=0.05, top=0.95)

    ax1 = fig.add_subplot(gs[:, :plot_width])
    plot_grid = np.ma.masked_where(Depth<=0, chl_control[day,:,:])
    ax1.pcolormesh(plot_grid, vmin=0, vmax=5, cmap='turbo')
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_title('No Subglacial Discharge')
    # print('control',np.min(chl_control[0,:,:]), np.max(chl_control[0,:,:]))


    ax2 = fig.add_subplot(gs[:, plot_width+1:2*plot_width+1])
    plot_grid = np.ma.masked_where(Depth <= 0, chl_iceplume[day, :, :])
    ax2.pcolormesh(plot_grid, vmin=0, vmax=5, cmap='turbo')
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title('Subglacial Discharge Plume')
    # print('iceplume',np.min(chl_iceplume[0, :, :]), np.max(chl_iceplume[0, :, :]))

    ax3 = fig.add_subplot(gs[:, -1])
    x = np.array([0,1])
    y = np.linspace(0,5,100)
    X, Y = np.meshgrid(x, y)
    ax3.pcolormesh(x,y,Y, vmin=vmin, vmax=vmax, cmap='turbo')
    ax3.set_xticks([])
    ax3.set_ylabel('Chl-a Concentration (mg/cm$^3$)')
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position("right")

    date_str = '1993'+'{:02d}'.format(month)+'{:02d}'.format(day+1)
    plt.savefig(os.path.join(output_folder,'panels','Chl_comparison_'+date_str+'.png'))

    plt.close(fig)



config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/configurations/' \
             'downscale_darwin/'
output_folder='/Users/mhwood/Documents/Research/Projects/Disko Bay/Data/Ocean/L2/Chlorophyll Iceplume Comparison'
model_name = 'L2_Disko_Bay'


XC, YC, Z, Depth = read_grid_geometry_from_nc(config_dir, model_name)

for month in range(8,9):
    iterations_control, chl_control = read_Chl_from_file(config_dir, model_name, month, experiment_name='control')

    iterations_iceplume, chl_iceplume = read_Chl_from_file(config_dir, model_name, month, experiment_name='iceplume')

    if month in [1,3,5,7,8,10,12]:
        n_days = 31
    elif month in [4,6,9,11]:
        n_days = 30
    else:
        n_days = 28
    for day in range(n_days):
        plot_Chl_panel(output_folder, month, day, XC, YC, Depth,
                       iterations_control, chl_control, iterations_iceplume, chl_iceplume)
