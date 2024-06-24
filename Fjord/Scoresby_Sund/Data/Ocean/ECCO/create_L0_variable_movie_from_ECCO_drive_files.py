
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from pyproj import Transformer
import cmocean.cm as cm
import shutil
import argparse
import sys

def read_ecco_geometry(ecco_dir, sNx, sNy,
                       ordered_ecco_tiles, ordered_tile_min_rows, ordered_tile_min_cols, ordered_ecco_tile_rotations):

    XC_grid = np.zeros((len(ordered_ecco_tiles) * sNx, len(ordered_ecco_tiles[1]) * sNy))
    YC_grid = np.zeros((len(ordered_ecco_tiles) * sNx, len(ordered_ecco_tiles[1]) * sNy))

    print('   - Reading geometry')
    for r in range(len(ordered_ecco_tiles)):
        for c in range(len(ordered_ecco_tiles[0])):
            # print('     - Reading grid '+str(r)+', '+str(c))
            tile = ordered_ecco_tiles[r][c]

            file_path = os.path.join(ecco_dir, 'GRID', 'GRID.' + '{:04d}'.format(tile+1) + '.nc')
            # print('Reading from grid GRID.' + '{:04d}'.format(tile+1) + '.nc')
            ds = nc4.Dataset(file_path)
            XC = ds.variables['XC'][:, :]
            YC = ds.variables['YC'][:, :]
            ds.close()

            min_row = ordered_tile_min_rows[r][c]
            min_col = ordered_tile_min_cols[r][c]
            rotations = ordered_ecco_tile_rotations[r][c]

            XC = XC[min_row:min_row + sNy, min_col:min_col + sNx]
            YC = YC[min_row:min_row + sNy, min_col:min_col + sNx]

            for rotation in range(rotations):
                XC = np.rot90(XC)
                YC = np.rot90(YC)

            XC_grid[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = XC
            YC_grid[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = YC

    # plt.subplot(1,2,1)
    # C = plt.imshow(XC_grid,origin='lower')
    # plt.colorbar(C)
    #
    # plt.subplot(1, 2, 2)
    # C = plt.imshow(YC_grid,origin='lower')
    # plt.colorbar(C)
    #
    # plt.show()

    return(XC_grid, YC_grid)


def read_annual_ecco_field(ecco_dir, field_name, year, sNx, sNy, depth_level,
                    ordered_ecco_tiles, ordered_tile_min_rows, ordered_tile_min_cols, ordered_ecco_tile_rotations):

    total_grid = np.zeros((12, len(ordered_ecco_tiles)*sNx, len(ordered_ecco_tiles[1])*sNy))
    total_time = np.zeros((12,))

    print('   - Reading file in year '+str(year))
    file_path = os.path.join(ecco_dir,field_name,field_name+'_'+str(year)+'.nc')
    ds = nc4.Dataset(file_path)
    grid = ds.variables[field_name][:,:,:,:,:]
    grid = grid[:,depth_level,:,:,:]
    for r in range(len(ordered_ecco_tiles)):
        for c in range(len(ordered_ecco_tiles[0])):
            # print('     - Reading grid '+str(r)+', '+str(c))
            tile = ordered_ecco_tiles[r][c]
            min_row = ordered_tile_min_rows[r][c]
            min_col = ordered_tile_min_cols[r][c]
            rotations = ordered_ecco_tile_rotations[r][c]
            grid_subset = grid[:, tile, min_row:min_row+sNy, min_col:min_col+sNx]

            for rotation in range(rotations):
                grid_subset = np.rot90(grid_subset,axes=(1,2))

            total_grid[:,r*sNy:(r+1)*sNy,c*sNx:(c+1)*sNx] = grid_subset


    time = year+np.arange(0.5, 1.5, 1/12)
    total_time[:12] = time

    ds.close()

    # plt.imshow(total_grid[0,:,:],origin='lower')
    # plt.show()

    return (total_time, total_grid)

def reproject_points(points,inputCRS,outputCRS,x_column=0,y_column=1):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(points[:, x_column], points[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
        run_test = False
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon = np.copy(points)
    output_polygon[:, x_column] = x2
    output_polygon[:, y_column] = y2
    return output_polygon


def compile_panels_to_movie(config_dir,config_name,field_name):
    pwd = os.getcwd()

    panels_dir = os.path.join(config_dir,'L0','plots',config_name,'ECCO',field_name)

    file_name = 'L0_'+field_name+'.mp4'

    os.chdir(panels_dir)
    os.system("ffmpeg -r 5 -i CE_Greenland_"+field_name+"_%04d.png -vcodec mpeg4 -b 1M -y " + file_name)
    os.rename(file_name, os.path.join('..', file_name))

    os.chdir(pwd)


def plot_mnc_fields(config_dir,ecco_dir,field_name,remove_old,skip):

    config_name = 'CE_Greenland'

    sNx = 90
    sNy = 90
    depth_level = 18

    ordered_ecco_tiles = [[10, 2, 2], [10, 6, 6]]
    ordered_tile_min_rows = [[2 * sNy, 2 * sNy, 2 * sNy], [2 * sNy, 2 * sNy, sNy]]
    ordered_tile_min_cols = [[0, 0, sNx], [0, 0, 0]]
    ordered_ecco_tile_rotations = [[1, 0, 0], [2, 3, 3]]  # rotations are counter-clockwise

    XC, YC = read_ecco_geometry(ecco_dir, sNx, sNy,
                       ordered_ecco_tiles, ordered_tile_min_rows, ordered_tile_min_cols, ordered_ecco_tile_rotations)

    if config_name == 'CE_Greenland':
        XC = XC[:,60:]
        YC = YC[:,60:]

    points = np.column_stack([XC.ravel(), YC.ravel()])
    points = reproject_points(points, inputCRS=4326, outputCRS=3413)
    XC = np.reshape(points[:, 0], np.shape(XC))
    YC = np.reshape(points[:, 1], np.shape(YC))

    years = np.arange(1992,2018)

    for year in years:

        counter = (year-1992)*12

        time, field_grid = read_annual_ecco_field(ecco_dir, field_name, year, sNx, sNy, depth_level,
                                     ordered_ecco_tiles, ordered_tile_min_rows, ordered_tile_min_cols, ordered_ecco_tile_rotations)

        if config_name == 'CE_Greenland':
            field_grid = field_grid[:,:,60:]

        min_x = 407815
        max_x = 1419944
        min_y = -2510694
        max_y = -923218

        output_dir = os.path.join(config_dir,'L0','plots',config_name,'ECCO')
        if field_name not in os.listdir(output_dir):
            os.mkdir(os.path.join(output_dir,field_name))
        output_dir = os.path.join(output_dir,field_name)

        if remove_old:
            os.system('rm -rf '+output_dir+'/*')

        if field_name == 'THETA':
            vmin = -1.9
            vmax = 8
            cmap = cm.thermal
        # if field_name == 'SALT':
        #     vmin = 32.5
        #     vmax = 35.2
        #     cmap = cm.haline
        # if field_name == 'UVEL' or field_name == 'VVEL':
        #     vmin = -1
        #     vmax = 1
        #     cmap = cm.balance
        # if field_name == 'ETAN':
        #     vmin = -3
        #     vmax = 0
        #     cmap = 'viridis'
        # if field_name == 'SPEED':
        #     vmin = 0
        #     vmax = 1
        #     cmap = cm.tempo_r
        # if field_name == 'VORTICITY':
        #     vmin = -0.25
        #     vmax = 0.25
        #     cmap = cm.curl
        # if field_name == 'SIarea':
        #     vmin = 0
        #     vmax = 1
        #     cmap = cm.ice

        panel_numbers = np.arange(0,np.shape(field_grid)[0],skip)

        for i in panel_numbers:

            fig = plt.figure(figsize=(6,8))
            plt.style.use('dark_background')

            gs = GridSpec(12, 10, left=0.05, right=0.99, hspace=0.15, top = 0.95, bottom = 0.05)
            ax1 = fig.add_subplot(gs[:-1, :])
            C = ax1.pcolormesh(XC,YC,field_grid[i, :, :],shading='nearest',cmap = cm.thermal,vmin = -1, vmax = 4)
            plt.colorbar(C)

            ax1.set_xlim([min_x, max_x])
            ax1.set_ylim([min_y, max_y])
            ax1.set_xticks([])
            ax1.set_yticks([])

            plt.title(field_name + ', depth = 257 m')

            ax2 = fig.add_subplot(gs[-1, :-2])

            rect = Rectangle((1992,0),time[i]-1992,1,fc='silver')
            ax2.add_patch(rect)

            ax2.set_ylim([0,1])
            ax2.set_xlim([1992,2018])
            ax2.set_xticks(np.arange(1992,2018,5))
            ax2.set_yticks([])
            plt.xlabel('Year')

            output_path = os.path.join(output_dir,config_name+'_'+field_name+'_'+'{:04d}'.format(counter)+'.png')
            plt.savefig(output_path)
            plt.close(fig)
            counter += 1

    compile_panels_to_movie(config_dir, config_name, field_name)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L0, L0, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_dir", action="store",
                        help="The directory where the ECCO data are stored.", dest="ecco_dir",
                        type=str, required=True)

    parser.add_argument("-f", "--field_name", action="store",
                        help="The name of the field to plot.", dest="field_name",
                        type=str, required=True)

    parser.add_argument("-a", "--plot_anomaly", action="store",
                        help="Choose whether to plot the anomaly or the actual field.", dest="plot_anomaly",
                        type=int, required=False)

    parser.add_argument("-r", "--remove_old", action="store",
                        help="Choose whether to remove old files (1 is true, 0 is false).", dest="remove_old",
                        type=int, required=False, default = 0)

    parser.add_argument("-s", "--skip", action="store",
                        help="Choose how many panels to skip at a time.", dest="skip",
                        type=int, required=False, default=1)

    args = parser.parse_args()
    config_dir = args.config_dir
    ecco_dir = args.ecco_dir
    field_name = args.field_name
    remove_old = args.remove_old
    skip = args.skip

    if remove_old==0:
        remove_old = False
    else:
        remove_old = True

    plot_mnc_fields(config_dir,ecco_dir,field_name,remove_old,skip)
   

