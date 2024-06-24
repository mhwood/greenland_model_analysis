
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

def create_L2_basemap_domain_plot(output_dir, config_dir, level_name, model_name):

    XC, YC, Depth, polygon = read_geometry_from_grid_nc(config_dir, level_name, model_name)


    output_path = os.path.join(output_dir,'L2_domain_basemap.png')

    plot_mode = 'dark'

    with_land = True

    fig = plt.figure(figsize=(12, 12))
    if plot_mode == 'dark':
        plt.style.use("dark_background")

    # m = Basemap(projection='ortho', resolution=None,
    #             lat_0=center_lat + rotation_lat, lon_0=center_lon + rotation_lon)
    m = Basemap(width=12000000,height=8000000,
            resolution='l',projection='stere',\
            lat_ts=70,lat_0=90,lon_0=-45.)
    if with_land:
        m.bluemarble(scale=20)

    polygon_lon, polygon_lat = m(polygon[:, 0], polygon[:, 1])
    m.plot(polygon_lon, polygon_lat, 'r-', linewidth=3)

    print(np.min(polygon[:,0]), np.max(polygon[:,0]))
    print(np.min(polygon[:, 1]), np.max(polygon[:, 1]))

    min_x1, min_y1 = m(np.min(polygon[:,0]), np.min(polygon[:,1]))
    max_x1, max_y1 = m(np.max(polygon[:, 0]), np.max(polygon[:, 1]))

    min_x2, min_y2 = m(np.min(polygon[:, 0]), np.max(polygon[:, 1]))
    max_x2, max_y2 = m(np.max(polygon[:, 0]), np.min(polygon[:, 1]))

    min_x = np.min([min_x1, min_x2])
    min_y = np.min([min_y1, min_y2])
    max_x = np.max([max_x1, max_x2])
    max_y = np.max([max_y1, max_y2])

    plt.gca().set_xlim([np.min(polygon_lon), np.max(polygon_lon)])
    plt.gca().set_ylim([np.min(polygon_lat), np.max(polygon_lat)])

    # axicon = fig.add_axes([0.22, 0.9, 0.5, 0.1])
    # # plt.text(0, 0, 'Global Source Model (1/3$^{\circ}$)', fontsize=20)
    # axicon.axis('off')
    # axicon.set_xticks([])
    # axicon.set_yticks([])

    plt.savefig(output_path)
    plt.close(fig)

config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/configurations/' \
             'downscale_darwin'

output_dir = '/Users/mhwood/Documents/Research/Projects/Disko Bay/Figures'

level_name = 'L2'
model_name = 'L2_Disko_Bay'

create_L2_basemap_domain_plot(output_dir, config_dir, level_name, model_name)