
import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
# from PIL import Image, ImageDraw
# import cmocean.cm as cm
# import argparse


def read_geometry_from_grid_nc(config_dir,level_name,model_name):

    grid_path = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')

    ds = nc4.Dataset(grid_path)
    Depth = ds.variables['Depth'][:, :]
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    ds.close()

    bottom = np.column_stack([XC[0, :], YC[0, :]])
    top = np.column_stack([XC[-1, :], YC[-1, :]])
    left = np.column_stack([XC[:, 0], YC[:, 0]])
    right = np.column_stack([XC[:, -1], YC[:, -1]])
    polygon = np.vstack([bottom, right, np.flipud(top), np.flipud(left)])

    return(XC, YC, Depth, polygon)

def generate_global_plot(output_path,
                  center_lon, rotation_lon, center_lat, rotation_lat,
                  polygon):

    plot_mode = 'dark'

    with_land = True

    fig = plt.figure(figsize=(12,12))
    if plot_mode=='dark':
        plt.style.use("dark_background")

    m = Basemap(projection='ortho', resolution=None,
                lat_0=center_lat + rotation_lat, lon_0=center_lon + rotation_lon)
    if with_land:
        m.bluemarble(scale=0.5)

    polygon_lon, polygon_lat = m(polygon[:, 0], polygon[:, 1])
    m.plot(polygon_lon, polygon_lat, 'r-', linewidth=3)

    axicon = fig.add_axes([0.22, 0.9, 0.5, 0.1])
    # plt.text(0, 0, 'Global Source Model (1/3$^{\circ}$)', fontsize=20)
    axicon.axis('off')
    axicon.set_xticks([])
    axicon.set_yticks([])

    plt.savefig(output_path)
    plt.close(fig)

def create_globe_domain_plot(output_dir, config_dir, level_name, model_name):

    XC, YC, Depth, polygon = read_geometry_from_grid_nc(config_dir, level_name, model_name)

    ##################################
    # some metadata

    center_lon = np.mean(XC)
    center_lat = np.mean(YC)

    rotation_lon = 45
    rotation_lat = -20

    # ################################################################################################
    # # read in the global grid
    #
    # file_path = os.path.join(config_dir,'L0_540','input','bathy_llc540')
    # L0_var_grid_faces = read_L0_var(file_path,'EtaN')
    #
    ##################################
    # make the plot

    output_path = os.path.join(output_dir,'L0_global_map.png')

    generate_global_plot(output_path,
                  center_lon, rotation_lon, center_lat, rotation_lat,
                  polygon)

config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/configurations/' \
             'downscale_darwin'

output_dir = '/Users/mhwood/Documents/Research/Projects/Disko Bay/Figures'

level_name = 'L2'
model_name = 'L2_Disko_Bay'

create_globe_domain_plot(output_dir, config_dir, level_name, model_name)