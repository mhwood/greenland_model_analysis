

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime

def read_model_grid(config_dir):
    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids','L2_Kanger_grid.nc'))
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    ds.close()
    return(XC, YC)

def read_iceplume_mask(config_dir,XC, YC):
    mask_file = os.path.join(config_dir, 'L2', 'L2_Kanger','input', 'dv', 'melange_mask.bin')
    mask = np.fromfile(mask_file,'>f4')
    mask = np.reshape(mask, np.shape(XC))
    return(mask)

def read_shelfice_melt_from_nc(config_dir, results_dir):
    nc_file = os.path.join(config_dir,'L2','L2_Kanger',results_dir,'dv','L2_Kanger_shelfice_SHFFWFLX.nc')
    ds = nc4.Dataset(nc_file)
    melange_melt = ds['SHFFWFLX'][:, :]
    melange_X = ds['X'][:, :]
    melange_Y = ds['Y'][:, :]
    iterations = ds['iterations'][:]
    ds.close()
    return(iterations, melange_X, melange_Y, melange_melt)

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

def create_melt_timeseries(iterations, melange_melt):

    mask = melange_melt[-50, :, :]!=0
    print(np.sum(mask))
    print(np.sum(melange_melt))

    time = np.zeros((len(iterations),))
    melt_mean = np.zeros((len(iterations),))
    melt_median = np.zeros((len(iterations),))
    melt_std = np.zeros((len(iterations),))

    for i in range(len(iterations)):
        time[i] = iter_number_to_dec_yr(iterations[i])
        melt_points = melange_melt[i, :, :]
        melt_points = -1*melt_points[mask]
        melt_mean[i] = np.mean(melt_points)
        melt_median[i] = np.median(melt_points)
        melt_std[i] = np.std(melt_points)

    return(time, melt_mean, melt_median, melt_std)

def write_timeseries_to_nc(project_dir,
                           time_no_plume, melt_mean_no_plume, melt_median_no_plume, melt_std_no_plume,
                           time_plume, melt_mean_plume, melt_median_plume, melt_std_plume):

    output_file = os.path.join(project_dir,'Data','Modeling','L2_Kanger_melange_melt_rates.nc')

    ds = nc4.Dataset(output_file, 'w')

    grp = ds.createGroup('no_plume')
    grp.createDimension('time',len(time_no_plume))
    t = grp.createVariable('time','f4',('time',))
    t[:] = time_no_plume
    t = grp.createVariable('melt_mean', 'f4', ('time',))
    t[:] = melt_mean_no_plume
    t = grp.createVariable('melt_median', 'f4', ('time',))
    t[:] = melt_median_no_plume
    t = grp.createVariable('melt_std', 'f4', ('time',))
    t[:] = melt_std_no_plume

    grp = ds.createGroup('plume')
    grp.createDimension('time', len(time_plume))
    t = grp.createVariable('time', 'f4', ('time',))
    t[:] = time_plume
    t = grp.createVariable('melt_mean', 'f4', ('time',))
    t[:] = melt_mean_plume
    t = grp.createVariable('melt_median', 'f4', ('time',))
    t[:] = melt_median_plume
    t = grp.createVariable('melt_std', 'f4', ('time',))
    t[:] = melt_std_plume

    ds.close()


def collect_melange_melt():
    config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
                 'configurations/downscaled_greenland'

    project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

    XC, YC = read_model_grid(config_dir)

    melange_mask = read_iceplume_mask(config_dir,XC, YC)

    results_dir = 'results_melange'
    iterations_no_plume, melange_X, melange_Y, melange_melt_no_plume = read_shelfice_melt_from_nc(config_dir, results_dir)
    time_no_plume, melt_mean_no_plume, melt_median_no_plume, melt_std_no_plume =\
        create_melt_timeseries(iterations_no_plume, melange_melt_no_plume)

    # plt.plot(time_no_plume,melt_mean_no_plume)
    # plt.show()

    results_dir = 'results_melange_plume'
    iterations_plume, _, _, melange_melt_plume = read_shelfice_melt_from_nc(config_dir, results_dir)
    time_plume, melt_mean_plume, melt_median_plume, melt_std_plume =\
        create_melt_timeseries(iterations_plume, melange_melt_plume)

    # plt.plot(time_plume,melt_mean_plume)
    # plt.show()

    write_timeseries_to_nc(project_dir,
                               time_no_plume, melt_mean_no_plume, melt_median_no_plume, melt_std_no_plume,
                               time_plume, melt_mean_plume, melt_median_plume, melt_std_plume)


if __name__ == '__main__':
    collect_melange_melt()





