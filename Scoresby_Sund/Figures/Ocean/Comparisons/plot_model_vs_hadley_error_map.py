
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import netCDF4 as nc4



def retrieve_error_grid(project_folder,var_name, min_depth, max_depth):

    output_file = os.path.join(project_folder, 'Data', 'Modeling', 'Downscaled', 'L1_CE_Greenland',
                               'L1_CE_Greenland_vs_Hadley_Analysis_' + str(min_depth) + 'm_' + str(max_depth) + 'm.nc')

    ds = nc4.Dataset(output_file)
    Lon = ds.variables['Lon'][:,:]
    Lat = ds.variables['Lat'][:,:]
    hadley_grid = ds.groups['Hadley'].variables[var_name][:,:,:]
    hadley_time = ds.groups['Hadley'].variables['time'][:]
    model_grid = ds.groups['Model'].variables[var_name][:, :, :]
    model_time = ds.groups['Model'].variables['time'][:]
    error_grid = ds.groups['Error'].variables[var_name][:,:,:]
    error_time = ds.groups['Error'].variables['time'][:]
    ds.close()

    return(Lon, Lat, hadley_time, hadley_grid, model_time, model_grid, error_time, error_grid)

def read_model_grid():
    grid_file = '/Volumes/mhwood/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
                'downscaled_greenland/nc_grids/L1_CE_Greenland_grid.nc'
    ds = nc4.Dataset(grid_file)
    xc = ds.variables['XC'][:, :]
    yc = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:,:]
    ds.close()

    return(xc,yc,Depth)


def plot_error_map(project_folder, var_name, XC, YC, Depth,
                   Lon, Lat, hadley_time, hadley_grid, model_time, model_grid, error_time, error_grid):

    fig = plt.figure(figsize=(10,12))

    gs1 = GridSpec(18, 13, left=0.09, right=0.95, wspace=0.05)

    ####################################################################################
    # Plot the ECCO Hadley Comparison

    ax1 = fig.add_subplot(gs1[:6, :6])

    # plot_error_map = np.sqrt(np.mean(error_grid ** 2, axis=0))
    # plot_error_map = np.ma.masked_where(plot_error_map == 0, plot_error_map)
    # C = plt.pcolormesh(Lon, Lat, plot_error_map, cmap='turbo', shading='nearest', vmin=0, vmax=3)
    # cbar = plt.colorbar(C)
    # cbar.set_label('RMSE ($^{\circ}$C)')

    plt.contour(XC, YC, Depth, linewidths=0.5, colors='k', levels=[1, 500])

    plt.plot(XC[:, 0], YC[:, 0], 'k-')
    plt.plot(XC[:, -1], YC[:, -1], 'k-')
    plt.plot(XC[0, :], YC[0, :], 'k-')
    plt.plot(XC[-1, :], YC[-1, :], 'k-')

    plt.ylabel('Longitude')
    plt.xlabel('Latitude')

    plt.gca().set_xlim([np.min(XC), np.max(XC)])
    plt.gca().set_ylim([np.min(YC), np.max(YC)])
    plt.title('ECCOv5 Error')

    ####################################################################################
    # Plot the L1 Hadley Comparison

    ax1 = fig.add_subplot(gs1[:6, 7:])

    plot_error_map = np.sqrt(np.mean(error_grid**2,axis=0))
    plot_error_map = np.ma.masked_where(plot_error_map==0,plot_error_map)
    C = plt.pcolormesh(Lon,Lat,plot_error_map, cmap='turbo', shading='nearest',vmin=0,vmax=3)
    cbar = plt.colorbar(C)
    cbar.set_label('RMSE ($^{\circ}$C)')

    plt.contour(XC,YC,Depth,linewidths=0.5,colors='k',levels=[1,500])

    plt.plot(XC[:, 0], YC[:, 0], 'k-')
    plt.plot(XC[:, -1], YC[:, -1], 'k-')
    plt.plot(XC[0, :], YC[0, :], 'k-')
    plt.plot(XC[-1, :], YC[-1, :], 'k-')

    plt.ylabel('Longitude')
    plt.xlabel('Latitude')

    plt.gca().set_xlim([np.min(XC),np.max(XC)])
    plt.gca().set_ylim([np.min(YC), np.max(YC)])
    plt.title('L1_CE_Greenland Error')

    ####################################################################################
    # Plot an example timeseries

    lon = -9
    lat = 76
    Dist = ((Lon-lon)**2 + (Lat-lat)**2)**0.5
    row, col = np.where(Dist==np.min(Dist))
    row = row[0]
    col = col[0]

    ax2 = fig.add_subplot(gs1[7:10, :])
    ax2.plot(model_time,model_grid[:,row,col],label='L1',color='purple')
    ax2.plot(hadley_time, hadley_grid[:, row, col], label='EN4 OA',color='k')
    plt.legend(loc=3)
    rmse = np.sqrt(np.mean((error_grid[:, row, col])**2))
    plt.title('Timeseries at '+str(lon)+'$^{\circ}$E, '+str(lat)+'$^{\circ}$N    RMSE = '+str(rmse)+'$^{\circ}$C')
    plt.ylabel('Temperature ($^{\circ}$C)')
    plt.grid(linewidth=0.5, alpha=0.5, linestyle='--')

    ####################################################################################
    # Plot another example timeseries

    lon = -18
    lat = 70
    Dist = ((Lon - lon) ** 2 + (Lat - lat) ** 2) ** 0.5
    row, col = np.where(Dist == np.min(Dist))
    row = row[0]
    col = col[0]

    ax2 = fig.add_subplot(gs1[11:14, :])
    ax2.plot(model_time, model_grid[:, row, col], label='L1', color='purple')
    ax2.plot(hadley_time, hadley_grid[:, row, col], label='EN4 OA', color='k')
    rmse = np.sqrt(np.mean((error_grid[:, row, col]) ** 2))
    plt.title('Timeseries at ' + str(lon) + '$^{\circ}$E, ' + str(lat) + '$^{\circ}$N    RMSE = '+str(rmse)+'$^{\circ}$C')
    plt.ylabel('Temperature ($^{\circ}$C)')
    plt.grid(linewidth=0.5, alpha=0.5, linestyle='--')

    ####################################################################################
    # Plot another example timeseries

    lon = -30
    lat = 65
    Dist = ((Lon - lon) ** 2 + (Lat - lat) ** 2) ** 0.5
    row, col = np.where(Dist == np.min(Dist))
    row = row[0]
    col = col[0]

    ax2 = fig.add_subplot(gs1[15:18, :])
    ax2.plot(model_time, model_grid[:, row, col], label='L1', color='purple')
    ax2.plot(hadley_time, hadley_grid[:, row, col], label='EN4 OA', color='k')
    rmse = np.sqrt(np.mean((error_grid[:, row, col]) ** 2))
    plt.title('Timeseries at ' + str(lon) + '$^{\circ}$E, ' + str(lat) + '$^{\circ}$N    RMSE = '+str(rmse)+'$^{\circ}$C')
    plt.ylabel('Temperature ($^{\circ}$C)')
    plt.grid(linewidth=0.5,alpha=0.5,linestyle='--')

    plt.suptitle('Model Comparison with EN4.2.2 Objective Analysis')

    plt.savefig(os.path.join(project_folder,'Figures','Ocean','L1_CE_Greenland','L1_CE_Greenland_'+str(var_name)+'Hadley_Error.png'))
    plt.close(fig)
    a=1

project_folder = '/Users/michwood/Documents/Research/Projects/Scoresby Sund'
var_name = 'Theta'
min_depth = 200
max_depth = 500

Lon, Lat, hadley_time, hadley_grid, model_time, model_grid, error_time, error_grid = \
    retrieve_error_grid(project_folder,var_name, min_depth, max_depth)

XC, YC, Depth = read_model_grid()

# plot it with coastline and everything
plot_error_map(project_folder, var_name, XC, YC, Depth,
                   Lon, Lat, hadley_time, hadley_grid, model_time, model_grid, error_time, error_grid)

