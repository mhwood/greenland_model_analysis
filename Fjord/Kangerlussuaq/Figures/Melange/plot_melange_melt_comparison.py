
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import shutil
from matplotlib.gridspec import GridSpec

def read_melt_rates_from_nc(project_dir):

    output_file = os.path.join(project_dir, 'Data', 'Modeling', 'L2_Kanger_melange_melt_rates.nc')

    ds = nc4.Dataset(output_file)

    grp = ds.groups['no_plume']
    time = grp.variables['time'][:]
    mean_melt = grp.variables['melt_mean'][:]
    melt_timeseries_no_plume = np.column_stack([time,mean_melt])

    grp = ds.groups['plume']
    time = grp.variables['time'][:]
    mean_melt = grp.variables['melt_mean'][:]
    melt_timeseries_plume = np.column_stack([time, mean_melt])

    ds.close()

    return(melt_timeseries_no_plume, melt_timeseries_plume)


def plot_melange_melt_timeseries(project_dir, melt_timeseries_no_plume, melt_timeseries_plume):

    fontsize = 14

    fig = plt.figure(figsize=(10,5))

    gs = GridSpec(1,1, left=0.09, right=0.95, bottom=0.08, top=0.95)

    ax = fig.add_subplot(gs[0, 0])

    melt_timeseries_no_plume[melt_timeseries_no_plume[:,1]==0,1] = np.nan
    melt_timeseries_plume[melt_timeseries_plume[:, 1] == 0, 1] = np.nan

    ax.plot(melt_timeseries_no_plume[:,0], melt_timeseries_no_plume[:,1]*86400/917, '-',
             color='green', label='Without Plume')
    ax.plot(melt_timeseries_plume[:, 0], melt_timeseries_plume[:, 1]*86400/917, '-',
             color='darkorange', label='With Plume')
    ax.set_ylabel('Vertical Melt Rate (m/day)', fontsize=fontsize)
    ax.set_title('Melange Submarine Melt Rate', fontsize=fontsize)
    ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)

    ax.text(2018.6,0.18, 'Model Output Error, Missing Data, Re-running Now', fontsize=16, color='darkorange',
            ha='center',va='center')

    ax.legend(loc='best', fontsize=fontsize)

    ax.set_xlim([2015,2022])

    output_file = os.path.join(project_dir,'Figures','Modeling','Kangerlussuaq Modeled Melange Melt.png')
    plt.savefig(output_file)
    plt.close(fig)


def plot_melange_comparison():
    project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

    melt_timeseries_no_plume, melt_timeseries_plume = read_melt_rates_from_nc(project_dir)


    plot_melange_melt_timeseries(project_dir, melt_timeseries_no_plume, melt_timeseries_plume)

    shutil.copyfile(os.path.join(project_dir,'Figures','Modeling','Kangerlussuaq Modeled Melange Melt.png'),
                    os.path.join(project_dir, 'Manuscript', 'Draft 1', 'Figure_7.png'))

if __name__ == '__main__':
    plot_melange_comparison()