

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime
from matplotlib.gridspec import GridSpec
import shutil

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)


def iter_number_to_dec_yr(iter_number,seconds_per_iter=60):

    total_seconds = iter_number*seconds_per_iter
    date = datetime.datetime(1992,1,1) + datetime.timedelta(seconds=total_seconds)

    dec_yr = YMD_to_DecYr(date.year,date.month,date.day)
    # print(date)
    return(dec_yr)


def read_melt_timeseries_from_nc(config_dir, results_dir, glacier):

    plume_index = 3

    file_name = os.path.join(config_dir,'L2','L2_Kanger',results_dir,'dv','L2_Kanger_iceplume_'+glacier+'.nc')
    ds = nc4.Dataset(file_name)
    depth = ds.variables['depth'][:]
    melt = ds.variables['ICEFRNTM'][:, :, :]
    iterations = ds.variables['iterations'][:]
    if 'melange' not in results_dir:
        max_depth = ds.max_depth
    else:
        max_depths = ds.variables['max_depths'][:]
        melange_drafts = ds.variables['drafts'][:]
    ds.close()


    melt = melt[:, :, plume_index].T
    if 'melange' in results_dir:
        max_depth = max_depths[plume_index]
        melange_draft = -1*melange_drafts[plume_index]
    else:
        melange_draft = 0

    dec_yrs = np.zeros((len(iterations),))
    for i in range(len(iterations)):
        dec_yrs[i] = iter_number_to_dec_yr(iterations[i])

    return(depth, dec_yrs, melt, max_depth, melange_draft)


def plot_melt_timeseries(project_dir, glacier, depth, dec_yrs_plume, melt_plume,
                     dec_yrs_melange_plume, melt_melange_plume, max_depth, melange_draft):

    min_year = 2015
    max_year = 2022

    vmin = 0
    vmax = 10

    dmin = -1
    dmax = 1

    # max_depth = 850
    # melange_draft = 240

    month_labels = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    # x_ticks = []
    # x_grid_locations = []
    # x_tick_labels = []
    # for year in range(min_year,max_year):
    #     for month in range(1,13):
    #         x_ticks.append(YMD_to_DecYr(year,month,15))
    #         x_tick_labels.append(month_labels[month-1])
    #         if month in [1,3,5,7,8,10,12]:
    #             x_grid_locations.append(YMD_to_DecYr(year,month,31))
    #         elif month in [4,6,9,11]:
    #             x_grid_locations.append(YMD_to_DecYr(year, month, 30))
    #         else:
    #             if year%4==0:
    #                 x_grid_locations.append(YMD_to_DecYr(year, month, 29))
    #             else:
    #                 x_grid_locations.append(YMD_to_DecYr(year, month, 28))



    fig = plt.figure(figsize=(8,10))

    gs = GridSpec(4, 23, top=0.96, bottom=0.05, left = 0.08, right=0.91)

    ax1 = fig.add_subplot(gs[0,:-2])
    melt_plume_plot = melt_plume#np.ma.masked_where(melt_plume==0, melt_plume)
    ax1.pcolormesh(dec_yrs_plume, depth, melt_plume_plot, vmin=vmin, vmax=vmax, cmap='turbo')
    # ax1.plot([min_year, max_year], [melange_draft, melange_draft], 'w-')
    # ax1.plot([min_year, max_year], [melange_draft, melange_draft],'k--',linewidth=1)
    ax1.set_xlim([min_year, max_year])
    ax1.set_ylim([max_depth,0])
    ax1.set_ylabel('Depth (m)')
    ax1.grid(linestyle='--',linewidth=0.5,alpha=0.5)
    ax1.set_title('Glacier Submarine Melt Rate (Without M\u00e9lange)')
    ax1.set_xticklabels([])
    ax1.text(2015.1, 30, 'a)', ha='left', va='top', fontsize=12)
    # ax1.set_xticks(x_ticks)
    # ax1.set_xticklabels(x_tick_labels)
    # for loc in range(len(x_grid_locations)):
    #     ax1.plot([x_grid_locations[loc], x_grid_locations[loc]], [max_depth,0], '-', linewidth=0.5, color='silver')

    ax1c = fig.add_subplot(gs[:2,-1])
    # ax1c.tick_params(axis='both', which='major', labelsize=12)
    cx = np.array([0, 1])
    cy = np.arange(vmin, vmax,0.01)
    CX, CY = np.meshgrid(cx, cy)
    ax1c.pcolormesh(CX, CY, CY, cmap='turbo')
    ax1c.set_xticks([])
    ax1c.set_ylabel('Melt Rate (m/day)')
    ax1c.yaxis.tick_right()
    ax1c.yaxis.set_label_position("right")

    ax2 = fig.add_subplot(gs[1, :-2])
    melt_melange_plume_plot = melt_melange_plume#np.ma.masked_where(melt_melange_plume == 0, melt_melange_plume)
    ax2.pcolormesh(dec_yrs_melange_plume, depth, melt_melange_plume_plot, vmin=vmin, vmax=vmax, cmap='turbo')
    ax2.plot([min_year, max_year], [melange_draft, melange_draft], 'w-')
    ax2.plot([min_year, max_year], [melange_draft, melange_draft], 'k--', linewidth=1)
    ax2.set_xlim([min_year, max_year])
    ax2.set_ylim([max_depth, 0])
    ax2.set_ylabel('Depth (m)')
    ax2.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax2.set_title('Glacier Submarine Melt Rate (With M\u00e9lange)')
    ax2.set_xticklabels([])
    ax2.text(2021.9, melange_draft-20, 'm\u00e9lange draft', ha='right', va='bottom', fontsize=12, color='white')
    ax2.text(2015.1, 30, 'b)', ha='left', va='top', fontsize=12,
                   bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    # ax2.set_xticks(x_ticks)
    # ax2.set_xticklabels(x_tick_labels)
    # for loc in range(len(x_grid_locations)):
    #     ax2.plot([x_grid_locations[loc], x_grid_locations[loc]], [max_depth, 0], '-', linewidth=0.5, color='silver')

    ax3 = fig.add_subplot(gs[2, :-2])
    difference = melt_melange_plume - melt_plume
    difference[melt_plume==0] = 0
    difference[melt_melange_plume==0] = 0
    difference_stats = difference[depth>300,:]
    mean_difference = np.mean(difference_stats[difference_stats!=0])
    std_difference = np.std(difference_stats[difference_stats != 0])
    print('Difference: '+str(mean_difference)+' +/- '+str(std_difference))
    ax3.pcolormesh(dec_yrs_melange_plume, depth, difference, vmin=dmin, vmax=dmax, cmap='seismic')
    ax3.text(2021.9, melange_draft - 20, 'm\u00e9lange draft', ha='right', va='bottom', fontsize=12)
    ax3.plot([min_year, max_year], [melange_draft, melange_draft], 'w-')
    ax3.plot([min_year, max_year], [melange_draft, melange_draft], 'k--', linewidth=1)
    ax3.set_xlim([min_year, max_year])
    ax3.set_ylim([max_depth, 0])
    ax3.set_ylabel('Depth (m)')
    ax3.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax3.set_title('Difference (With Melange - Without Melange)')
    ax3.set_xticklabels([])
    ax3.text(2015.1, 30, 'c)', ha='left', va='top', fontsize=12)
    # ax3.set_xticks(x_ticks)
    # ax3.set_xticklabels(x_tick_labels)
    # for loc in range(len(x_grid_locations)):
    #     ax3.plot([x_grid_locations[loc], x_grid_locations[loc]], [max_depth, 0], '-', linewidth=0.5, color='silver')

    ax2c = fig.add_subplot(gs[2, -1])
    # ax1c.tick_params(axis='both', which='major', labelsize=12)
    cx = np.array([0, 1])
    cy = np.arange(dmin, dmax, 0.01)
    CX, CY = np.meshgrid(cx, cy)
    ax2c.pcolormesh(CX, CY, CY, cmap='seismic')
    ax2c.set_xticks([])
    ax2c.set_ylabel('Melt Rate (m/day)')
    ax2c.yaxis.tick_right()
    ax2c.yaxis.set_label_position("right")

    max_melt_plume = np.max(melt_plume,axis=0)
    max_melt_melange_plume = np.max(melt_melange_plume, axis=0)

    ax4 = fig.add_subplot(gs[3,:-2])
    ax4.plot(dec_yrs_plume[max_melt_plume != 0], max_melt_plume[max_melt_plume != 0], label='plume only',color='purple')
    ax4.plot(dec_yrs_melange_plume[max_melt_melange_plume!=0],
             max_melt_melange_plume[max_melt_melange_plume!=0], label='m\u00e9lange + plume', color='darkorange')
    ax4.set_xlim([min_year, max_year])
    ax4.set_ylabel('Melt Rate (m/day)')
    ax4.set_title('Maximal Melt Rate')
    ax4.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax4.text(2015.1, 10.5, 'd)', ha='left', va='top', fontsize=12)
    # ax4.set_xticks(x_ticks)
    # ax4.set_xticklabels(x_tick_labels)
    # for loc in range(len(x_grid_locations)):
    #     plt.plot([x_grid_locations[loc], x_grid_locations[loc]], [6.5,10], '-', linewidth=0.5, color='silver')
    plt.legend()

    output_file = os.path.join(project_dir, 'Figures','Ocean', 'Modeling', 'Kangerlussuaq Modeled Ice Front Melt.png')
    plt.savefig(output_file)
    plt.close(fig)

    shutil.copyfile(os.path.join(project_dir, 'Figures', 'Ocean', 'Modeling', 'Kangerlussuaq Modeled Ice Front Melt.png'),
                    os.path.join(project_dir, 'Manuscript', 'Draft 1', 'Figure_6.png'))


def plot_melt_rate_comparison():
    config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
                 'configurations/downscaled_greenland'

    project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

    glacier = 'Kangerlussuaq'


    depth, dec_yrs_melange_plume, melt_melange_plume, max_depth, melange_draft = read_melt_timeseries_from_nc(config_dir,
                                                                                  'results_melange_plume', glacier)
    print(max_depth,melange_draft)

    _, dec_yrs_plume, melt_plume, _, _ = read_melt_timeseries_from_nc(config_dir,
                                                                                  'results_plume', glacier)

    plot_melt_timeseries(project_dir, glacier, depth, dec_yrs_plume, melt_plume,
                         dec_yrs_melange_plume, melt_melange_plume, max_depth, melange_draft)


if __name__ == '__main__':
    plot_melt_rate_comparison()
