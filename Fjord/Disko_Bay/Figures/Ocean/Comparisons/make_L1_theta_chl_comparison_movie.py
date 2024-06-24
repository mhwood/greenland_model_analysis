
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
import shutil
import cmocean.cm as cm
import argparse
from matplotlib.patches import Rectangle
import datetime as dt
from datetime import datetime, timedelta
from osgeo import gdal
import sys
import moviepy.video.io.ImageSequenceClip

def read_geometry_from_grid_nc(config_dir,level_name,model_name):

    grid_path = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')

    ds = nc4.Dataset(grid_path)
    Depth = ds.variables['Depth'][:, :]
    ds.close()

    return(Depth)

def read_background_imagery(file_path):
    ds = gdal.Open(file_path)
    R = np.array(ds.GetRasterBand(1).ReadAsArray())
    G = np.array(ds.GetRasterBand(2).ReadAsArray())
    B = np.array(ds.GetRasterBand(3).ReadAsArray())
    rows = np.shape(R)[0]
    cols = np.shape(R)[1]
    R = R.reshape((rows,cols, 1))
    G = G.reshape((rows, cols, 1))
    B = B.reshape((rows, cols, 1))
    image = np.concatenate([B,G,R],axis=2)
    brightness_factor = 0.1 # 0 - 1
    image = (np.max(image)-np.max(image)*(brightness_factor))/(np.max(image)-np.min(image))*(image-np.min(image))+np.max(image)*(brightness_factor)
    image[image<0]=0
    image[image>1]=1
    transform = ds.GetGeoTransform()
    extents = [transform[0],transform[0]+transform[1]*np.shape(image)[1],transform[3]+ transform[5] * np.shape(image)[0], transform[3]]
    # x_resolution = transform[1]
    # y_resolution = transform[5]
    return(image,extents)

def iter_number_to_date(iter_number):
    seconds_per_iter = 300
    total_seconds = iter_number*seconds_per_iter
    date = datetime(1992,1,1) + timedelta(seconds=total_seconds)
    return(date)

def date_to_iter_number(date):
    seconds_per_iter = 300
    total_seconds = (date-datetime(1992,1,1)).total_seconds()
    iter_number = total_seconds/seconds_per_iter
    return(iter_number)

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = dt.datetime(year,month,day,hour,minute,second)
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_monthly_file_name(config_dir, year, month):

    level_name = 'L1'
    model_name = 'L1_W_Greenland'

    theta_file_name = ''
    for file_name in os.listdir(os.path.join(config_dir,level_name,model_name,
                                             'results','TS_AW_daily_snap')):
        if file_name[0]!='.':
            if int(file_name.split('.')[1][:4]) == year and int(file_name.split('.')[1][4:6]) == month:
                theta_file_name = file_name

    chl_file_name = ''
    for file_name in os.listdir(os.path.join(config_dir, level_name, model_name,
                                             'results', 'BGC_daily_Chl')):
        if file_name[0] != '.':
            if int(file_name.split('.')[1][:4]) == year and int(file_name.split('.')[1][4:6]) == month:
                chl_file_name = file_name

    ds = nc4.Dataset(os.path.join(config_dir,level_name,model_name,
                                             'results','TS_AW_daily_snap', theta_file_name))
    iterations = ds.variables['iterations'][:]
    theta = ds.variables['Theta'][:]
    ds.close()

    ds = nc4.Dataset(os.path.join(config_dir, level_name, model_name,
                                  'results', 'BGC_daily_Chl', chl_file_name))
    chl = ds.variables['Chl01'][:]
    chl += ds.variables['Chl02'][:]
    chl += ds.variables['Chl03'][:]
    chl += ds.variables['Chl04'][:]
    chl += ds.variables['Chl05'][:]
    ds.close()

    return(iterations, theta, chl)

def create_panel_plot(output_dir, file_name, theta, chl,
                      iter_number, depth, add_background_imagery = False):

    date = iter_number_to_date(iter_number)
    dec_yr = YMD_to_DecYr(date.year, date.month, date.day)
    min_iter = date_to_iter_number(datetime(date.year, 1, 1))
    max_iter = date_to_iter_number(datetime(date.year + 1, 1, 1))

    if add_background_imagery:
        file_path = "/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/" \
                    "configurations/downscale_darwin/L1/L1_W_Greenland/plots/basemap/L1_W_Greenland_MODIS_20220720_row_col.tif"
        background_image, extents = read_background_imagery(file_path)

    ############################################################################
    # make the plot

    fig = plt.figure(figsize=(14, 10))
    plt.style.use('dark_background')

    gs2 = GridSpec(17, 29, left=0.05, right=0.95, top = 0.95, bottom = 0.05, hspace=0.05)

    ax1 = fig.add_subplot(gs2[:-4, :14])
    ax2 = fig.add_subplot(gs2[:-4, 15:])

    if add_background_imagery:
        rect = Rectangle((0,0),np.shape(depth)[1], np.shape(depth)[0],
                         facecolor='white', edgecolor='white', zorder=1)
        ax1.add_patch(rect)
        rect = Rectangle((0, 0), np.shape(depth)[1], np.shape(depth)[0],
                         facecolor='white', edgecolor='white', zorder=1)
        ax2.add_patch(rect)
        ax1.imshow(background_image, extent=extents, alpha=0.7, zorder=1)
        ax2.imshow(background_image, extent=extents, alpha=0.7, zorder=1)

    plot_grid_theta = np.copy(theta)
    # print('        - range: ',np.min(plot_grid),np.max(plot_grid))
    plot_grid_theta = np.ma.masked_where(depth <= 0, plot_grid_theta)
    # if 'AW' in field_name:
    #     plot_grid = np.ma.masked_where(plot_grid==0,plot_grid)
    #     # land_mask = np.ma.masked_where(land_mask,plot_grid==0)
    # else:
    #     plot_grid = np.ma.masked_where(depth == 0, plot_grid)
    C = ax1.imshow(plot_grid_theta, origin='lower', vmin=0, vmax=6, cmap='turbo', zorder=1)

    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_title('Potential Temperature (250 m)')

    plot_grid_chl = np.copy(chl)
    # print('        - range: ',np.min(plot_grid),np.max(plot_grid))
    plot_grid_chl = np.ma.masked_where(depth <= 0, plot_grid_chl)
    # if 'AW' in field_name:
    #     plot_grid = np.ma.masked_where(plot_grid==0,plot_grid)
    #     # land_mask = np.ma.masked_where(land_mask,plot_grid==0)
    # else:
    #     plot_grid = np.ma.masked_where(depth == 0, plot_grid)
    C = ax2.imshow(plot_grid_chl, origin='lower', vmin=0, vmax=4, cmap='turbo', zorder=1)

    # if 'AW' in field_name and add_background_imagery:
    #     ax1.imshow(aw_depth_mask, origin='lower', vmin=0, vmax=1, cmap='Greys')
    #     ax1.contour(depth,levels=[0.1],colors='w',linewidths=0.4)
    #
    # if not add_background_imagery:
    #     ax1.imshow(land_mask, origin='lower', vmin=0, vmax=2, cmap='Greys')

    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title('Chlorophyll-a (Surface)')

    # ax1.text(5, 5, 'Timestep: ' + str(int(iter_number)), ha='left', va='bottom', color='k')

    min_year = 1992
    max_year = 2022
    ax2 = fig.add_subplot(gs2[-1, 1:-1])
    width = (dec_yr - min_year)
    rect = Rectangle((min_year, 0), width, 1, fc='silver', ec='white')
    ax2.add_patch(rect)
    ax2.set_xlim([min_year, max_year])
    ax2.set_ylim([0, 1])
    for i in range(min_year,max_year,5):
        plt.plot([i,i],[0,1],'w-',linewidth=0.5)
    ax2.set_xticks(np.arange(min_year,max_year,5))
    ax2.set_yticks([])

    ax3 = fig.add_subplot(gs2[-3, 1:-1])
    width = (iter_number - min_iter) / (max_iter - min_iter)
    rect = Rectangle((date.year, 0), width, 1, fc='silver', ec='white')
    ax3.add_patch(rect)
    ax3.set_xlim([date.year, date.year + 1])
    ax3.set_ylim([0, 1])
    for i in range(1,12):
        plt.plot([date.year+i/12,date.year+i/12],[0,1],'w-',linewidth=0.5)
    ax3.set_xticks(np.arange(date.year+1/24, date.year + 1, 1 / 12))
    ax3.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
    ax3.set_yticks([])

    output_path = os.path.join(output_dir, 'panels', file_name)
    plt.savefig(output_path)
    plt.close(fig)


def create_movie(output_dir):
    # make a list of files for each movie panel
    file_list = []
    for file_name in os.listdir(output_dir + '/panels'):
        if file_name[-3:] == 'png':
            file_list.append(output_dir + '/panels/' + file_name)

    # sort the panels
    file_list.sort()

    # set the frames per second
    fps = 5

    # use the ImageSequenceClip module to set up the clip
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(file_list, fps=fps)

    # write the video to a file
    clip.write_videofile(output_dir+'/L1_Theta_Chl.mp4')


config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/' \
             'darwin3/configurations/downscale_darwin'

output_dir = '/Users/mhwood/Documents/Research/Projects/Disko Bay/Figures/L1/Theta Chlorophyll Comparison'

depth = read_geometry_from_grid_nc(config_dir,'L1','L1_W_Greenland')

year = 1997

# for month in range(1,13):
#     iterations, theta_grid, chl_grid = read_monthly_file_name(config_dir, year, month)
#
#     for i in range(len(iterations)):
#
#         iter_number = iterations[i]
#         theta = theta_grid[i, :, :]
#         chl = chl_grid[i, :, :]
#
#         file_name = 'L1_theta_chl_'+str(year)+'{:02d}'.format(month)+'{:02d}'.format(i+1)+'.png'
#         print(file_name)
#
#         create_panel_plot(output_dir, file_name, theta, chl,
#                           iter_number, depth, add_background_imagery=True)

create_movie(output_dir)





