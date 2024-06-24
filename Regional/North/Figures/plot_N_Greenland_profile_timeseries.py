
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import cmocean.cm as cm


def read_nc_profile(project_folder, model_name, location, var_name):

    timeseries_file = os.path.join(project_folder, 'Data', 'Ocean', 'L1', model_name+'_'+location+'_'+var_name+'_timeseries.nc')

    ds = nc4.Dataset(timeseries_file)
    time = ds.variables['time'][:]
    depth = ds.variables['depth'][:]
    var_grid = ds.variables[var_name][:]
    ds.close()

    if var_name=='Total_Chl':
        depth = depth[depth<100]
        var_grid = var_grid[:len(depth),:]

    return(time, depth, var_grid)

def create_profile_figure(project_folder, var_name, location, time, depth, var_grid):



    if var_name=='Theta':
        vmin = 0
        vmax = 1
        cmap='turbo'
    elif var_name=='Salt':
        vmin = 32
        vmax = 35
        cmap=cm.haline
    elif var_name=='Total_Chl':
        vmin = 0
        vmax = 2
        cmap='turbo'
    else:
        vmin = np.min(var_grid)
        vmax = np.max(var_grid)
        cmap = 'turbo'


    fig = plt.figure(figsize=(10,4))

    plt.style.use('dark_background')

    # C = plt.pcolormesh(time,depth,var_grid, shading='nearest', cmap='turbo', vmin=0, vmax=4)
    C = plt.contourf(time, depth, var_grid, 100, cmap=cmap, vmin=vmin, vmax=vmax, extend='max')

    cbar = plt.colorbar(C)

    # m = plt.cm.ScalarMappable(cmap=cmap)
    # m.set_array(var_grid)
    # m.set_clim(vmin, vmax)
    # cbar = plt.colorbar(m, boundaries=np.linspace(vmin,vmax, 100))
    cbar.set_label(var_name)

    plt.gca().invert_yaxis()
    # plt.colorbar(C)


    output_file = os.path.join(project_folder,'Figures','L1_N_Greenland_'+'_'+location+'_'+var_name+'_Profile_Timeseries.png')
    plt.savefig(output_file)
    plt.close(fig)
    a=1



project_folder = '/Users/mhwood/Documents/Research/Projects/North Greenland'

model_name = 'L1_N_Greenland'
location = 'Nares'
var_name = 'Theta'

time, depth, var_grid = read_nc_profile(project_folder, model_name, location, var_name)

create_profile_figure(project_folder, location, var_name, time, depth, var_grid)
