

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon


def read_transect_from_nc(config_dir, experiment, model_name, fjord_name, var_name):

    if experiment=='control':
        results_dir = 'results'
    else:
        results_dir = 'results_'+experiment

    print('  - Reading in the '+var_name+' transect from the '+experiment+' experiment')

    output_file = os.path.join(config_dir, 'L2', model_name, results_dir , 'dv',
                               model_name + '_' + fjord_name + '_Trough_' + var_name +'.nc')
    # print('Reading from '+output_file)

    ds = nc4.Dataset(output_file)
    distance = ds.variables['distance'][:]
    depth = ds.variables['depth'][:]
    # iterations = ds.variables['iterations'][:]
    # thickness = ds.variables['ice_thickness'][:]
    bathymetry = ds.variables['bathymetry'][:]
    var_grid = ds.variables[var_name.upper()][:, :, :]
    ds.close()

    distance = np.flip(distance)
    bathymetry = np.flip(bathymetry)
    var_grid = np.flip(var_grid, axis=-1)
    non_zero_indices = var_grid[:,0,0]!=0
    print('    - Limiting to '+str(np.sum(non_zero_indices))+' non-zero indices')
    var_grid = var_grid[non_zero_indices,:,:]
    var_grid = var_grid[:76, :, :]
    var_grid = np.nanmean(var_grid, axis=0)

    return(depth, distance, var_grid, bathymetry)


def plot_panel(project_dir, output_file_name, depth,
               distance_control, var_grid_control, bathymetry_control,
               distance_melange, var_grid_melange, bathymetry_melange,
               distance_plume, var_grid_plume, bathymetry_plume,
               distance_melange_plume, var_grid_melange_plume, bathymetry_melange_plume):

    fig = plt.figure(figsize=(10,12))

    plot_height = 8

    gs2 = GridSpec(plot_height*4+3, 23, left=0.1, right=0.95, top=0.95, bottom=0.05, hspace=0.05)

    vmin = -2
    vmax = 2.5

    dmin=-0.5
    dmax=0.5

    #############################################################################################################

    ax1 = fig.add_subplot(gs2[:plot_height, :-3])
    ax1.pcolormesh(distance_control, depth, var_grid_control, cmap='turbo', vmin=vmin, vmax=vmax)
    # ax1.plot(distance_control, thickness_control, 'k-')
    ax1.plot(distance_control, bathymetry_control, 'k-')

    bathy_outline_zach = np.vstack([np.column_stack([distance_control, bathymetry_control]),
                                    np.flipud(np.column_stack([distance_control, 2000*np.ones_like(bathymetry_control)]))])
    bathy_poly = Polygon(bathy_outline_zach, facecolor='silver')
    ax1.add_patch(bathy_poly)

    # ice_outline_zach = np.vstack([np.column_stack([distance_control, thickness_control]),
    #                                 np.flipud(np.column_stack([distance_control, -10 * np.ones_like(bathymetry_control)]))])
    # ice_poly = Polygon(ice_outline_zach, facecolor='white')
    # ax1.add_patch(ice_poly)
    ax1.text(np.max(distance_control)-2, 1000, 'Control', ha='right', va='bottom',
             color='k', fontsize=16)

    # ax1.invert_yaxis()
    # max_x_index = np.min(np.where(bathymetry_control == 0)[0] - 1)
    # ax1.set_xlim([distance_control[max_x_index], np.min(distance_control)+150])
    ax1.set_xticklabels([])
    ax1.set_ylim([1000, 0])
    ax1.set_ylabel('Depth (m)')
    ax1.set_title('Along-Fjord Transects')

    #############################################################################################################

    ax2 = fig.add_subplot(gs2[plot_height+1:1+2*plot_height, :-3])
    ax2.pcolormesh(distance_melange, depth, var_grid_melange, cmap='turbo', vmin=vmin, vmax=vmax)
    # ax2.plot(distance_melange, thickness_melange, 'k-')
    ax2.plot(distance_melange, bathymetry_melange, 'k-')

    bathy_outline_melange = np.vstack([np.column_stack([distance_melange, bathymetry_melange]),
                                    np.flipud(np.column_stack([distance_melange, 2000 * np.ones_like(bathymetry_melange)]))])
    bathy_poly = Polygon(bathy_outline_melange, facecolor='silver')
    ax2.add_patch(bathy_poly)

    # ice_outline_melange = np.vstack([np.column_stack([distance_melange, thickness_melange]),
    #                               np.flipud(np.column_stack([distance_melange, -10 * np.ones_like(bathymetry_melange)]))])
    # ice_poly = Polygon(ice_outline_melange, facecolor='white')
    # ax2.add_patch(ice_poly)
    # ax2.text(np.max(distance_melange),20,' 79N', ha='left', va='top', color='k')

    ax2.text(np.max(distance_control) - 2, 1000, 'Melange Only', ha='right', va='bottom',
             color='k', fontsize=16)

    # ax2.invert_yaxis()
    # ax2.set_xlim([distance_melange[max_x_index], np.min(distance_melange)])
    ax2.set_xticklabels([])
    ax2.set_ylim([1000, 0])
    ax2.set_ylabel('Depth (m)')

    #############################################################################################################

    ax3 = fig.add_subplot(gs2[2 + 2 * plot_height:2 + 3 * plot_height, :-3])
    ax3.pcolormesh(distance_plume, depth, var_grid_plume, cmap='turbo', vmin=vmin, vmax=vmax)
    # ax3.plot(distance_plume, thickness_plume, 'k-')
    ax3.plot(distance_plume, bathymetry_plume, 'k-')

    bathy_outline_plume = np.vstack([np.column_stack([distance_plume, bathymetry_plume]),
                                     np.flipud(np.column_stack(
                                         [distance_plume, 2000 * np.ones_like(bathymetry_plume)]))])
    bathy_poly = Polygon(bathy_outline_plume, facecolor='silver')
    ax3.add_patch(bathy_poly)

    # ice_outline_plume = np.vstack([np.column_stack([distance_plume, thickness_plume]),
    #                               np.flipud(np.column_stack([distance_plume, -10 * np.ones_like(bathymetry_plume)]))])
    # ice_poly = Polygon(ice_outline_plume, facecolor='white')
    # ax3.add_patch(ice_poly)
    # ax3.text(np.max(distance_plume), 20, ' 79N', ha='left', va='top', color='k')

    ax3.text(np.max(distance_control) - 2, 1000, 'Plume Only', ha='right', va='bottom',
             color='k', fontsize=16)

    # ax3.invert_yaxis()
    # ax3.set_xlim([distance_plume[max_x_index], np.min(distance_plume)])
    ax3.set_xticklabels([])
    ax3.set_ylim([1000, 0])
    ax3.set_ylabel('Depth (m)')

    #############################################################################################################

    ax4 = fig.add_subplot(gs2[3+3*plot_height:3+4*plot_height, :-3])
    ax4.pcolormesh(distance_melange_plume, depth, var_grid_melange_plume, cmap='turbo', vmin=vmin, vmax=vmax)
    # ax4.plot(distance_melange_plume, thickness_melange_plume, 'k-')
    ax4.plot(distance_melange_plume, bathymetry_melange_plume, 'k-')

    bathy_outline_melange_plume = np.vstack([np.column_stack([distance_melange_plume, bathymetry_melange_plume]),
                                       np.flipud(np.column_stack(
                                           [distance_melange_plume, 2000 * np.ones_like(bathymetry_melange_plume)]))])
    bathy_poly = Polygon(bathy_outline_melange_plume, facecolor='silver')
    ax4.add_patch(bathy_poly)

    # ice_outline_melange_plume = np.vstack([np.column_stack([distance_melange_plume, thickness_melange_plume]),
    #                               np.flipud(np.column_stack([distance_melange_plume, -10 * np.ones_like(bathymetry_melange_plume)]))])
    # ice_poly = Polygon(ice_outline_melange_plume, facecolor='white')
    # ax4.add_patch(ice_poly)
    # ax4.text(np.max(distance_melange_plume), 20, ' 79N', ha='left', va='top', color='k')

    ax4.text(np.max(distance_control) - 2, 1000, 'Melange and Plume', ha='right', va='bottom',
             color='k', fontsize=16)

    # ax4.invert_yaxis()
    # ax4.set_xlim([distance_melange_plume[max_x_index], np.min(distance_melange_plume)])
    ax4.set_ylim([1000, 0])
    ax4.set_ylabel('Depth (m)')
    ax4.set_xlabel('Distance From Kangerlussuaq Glacier (km)')

    #############################################################################################################

    ax3 = fig.add_subplot(gs2[1+plot_height:2+3*plot_height, -2])
    x = np.array([0, 1])
    y = np.linspace(vmin, vmax, 100)
    X, Y = np.meshgrid(x, y)
    ax3.pcolormesh(X, Y, Y, cmap='turbo')
    ax3.yaxis.tick_right()
    ax3.set_ylabel('Temperature ($^{\circ}$C)')
    ax3.yaxis.set_label_position("right")
    ax3.set_xticks([])

    #############################################################################################################

    output_file = os.path.join(project_dir, 'Figures','Ocean','Modeling',output_file_name)

    plt.savefig(output_file)
    plt.close(fig)


project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'


config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/' \
             'MITgcm/configurations/downscaled_greenland'

field_name = 'Theta'
L2_model_name = 'L2_Kanger'
output_file_name = 'Kangerlussuaq Fjord Mean '+str(field_name)+' Comparison.png'

depth, distance_control, var_grid_control, bathymetry_control = \
    read_transect_from_nc(config_dir, 'control', L2_model_name, 'Kanger', field_name.upper())
distance_control *= 1e-3
distance_control = np.max(distance_control)-distance_control

_, distance_melange, var_grid_melange, bathymetry_melange = \
    read_transect_from_nc(config_dir, 'melange', L2_model_name, 'Kanger', field_name.upper())
distance_melange *= 1e-3
distance_melange = np.max(distance_melange)-distance_melange

_, distance_plume, var_grid_plume, bathymetry_plume = \
    read_transect_from_nc(config_dir, 'plume', L2_model_name, 'Kanger', field_name.upper())
distance_plume *= 1e-3
distance_plume = np.max(distance_plume)-distance_plume

_, distance_melange_plume, var_grid_melange_plume, bathymetry_melange_plume = \
    read_transect_from_nc(config_dir, 'melange_plume', L2_model_name, 'Kanger', field_name.upper())
distance_melange_plume *= 1e-3
distance_melange_plume = np.max(distance_melange_plume)-distance_melange_plume

plot_panel(project_dir, output_file_name, depth,
               distance_control, var_grid_control, bathymetry_control,
               distance_melange, var_grid_melange, bathymetry_melange,
           distance_plume, var_grid_plume, bathymetry_plume,
           distance_melange_plume, var_grid_melange_plume, bathymetry_melange_plume)





