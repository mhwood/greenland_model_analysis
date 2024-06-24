
import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from scipy.interpolate import interp1d
import datetime

def iter_number_to_date(iter_number,model_level='L1'):

    if model_level=='L0':
        seconds_per_iter = 1200
    if model_level=='L1':
        seconds_per_iter = 300
    if model_level=='L1_150':
        seconds_per_iter = 150
    if model_level=='L1_100':
        seconds_per_iter = 100
    if model_level=='N1':
        seconds_per_iter = 120
    if model_level=='L3':
        seconds_per_iter = 60

    total_seconds = iter_number*seconds_per_iter
    date = datetime.datetime(1992,1,1) + datetime.timedelta(seconds=total_seconds)
    return(date)

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def find_profile_files_in_box(results_dir, min_lon, min_lat, max_lon, max_lat, var_name):

    first_file = True
    years = []
    months = []
    days = []
    dec_yrs = []

    if var_name=='Nitrate':
        read_name = 'NO3'
    if var_name=='Phosphate':
        read_name = 'PO4'

    for file_name in sorted(os.listdir(results_dir)):
        if file_name.endswith('.nc') and file_name[0]!='.':
            print('    - Reading from '+file_name)

            ds = nc4.Dataset(os.path.join(results_dir,file_name))
            depth = ds.variables['depths'][:]
            Lon = ds.variables['longitude'][:, :]
            Lat = ds.variables['latitude'][:, :]
            var_grid = ds.variables[read_name][:, :, :, :]
            iterations = ds.variables['iterations'][:]
            ds.close()

            bbox = np.array([[min_lon, min_lat],
                             [max_lon, min_lat],
                             [max_lon, max_lat],
                             [min_lon, max_lat],
                             [min_lon, min_lat]])

            p = mplPath.Path(bbox)
            full_grid = np.column_stack([Lon.ravel(), Lat.ravel()])
            inside = p.contains_points(full_grid)

            profile = np.zeros((len(depth),))
            for d in range(len(depth)):
                layer_grid = var_grid[0,d,:,:].ravel()
                inside_points = layer_grid[inside]
                if np.any(inside_points!=0):
                    profile[d] = np.mean(inside_points[inside_points!=0])

            # plt.plot(profile,depth)
            # plt.gca().invert_yaxis()
            # plt.show()

            date = iter_number_to_date(int(iterations[0]))
            dec_yr = YMD_to_DecYr(date.year,date.month,date.day)
            years.append(date.year)
            months.append(date.month)
            days.append(date.day)
            dec_yrs.append(dec_yr)

            if first_file:
                first_file = False
                output_grid = np.reshape(profile, (len(profile),1))
            else:
                output_grid = np.hstack([output_grid, np.reshape(profile, (len(profile),1))])

    return(depth, years, months, days, dec_yrs, output_grid)

def outputting_transects_to_nc(folder, var_name, years, months, days, dec_yrs,
                               depth, min_lon, min_lat, max_lon, max_lat, var_grid):

    ds = nc4.Dataset(os.path.join(folder,'Data','Ocean','L1','Disko_Bay', 'L1_'+var_name+'_Profiles_Disko_Bay.nc'),'w')

    ds.createDimension('time',len(years))
    ds.createDimension('depth',len(depth))

    var = ds.createVariable('year','i4',('time',))
    var[:] = years

    var = ds.createVariable('month', 'i4', ('time',))
    var[:] = months

    var = ds.createVariable('days', 'i4', ('time',))
    var[:] = days

    var = ds.createVariable('dec_yrs', 'i4', ('time',))
    var[:] = dec_yrs

    var = ds.createVariable('depth', 'f4', ('depth',))
    var[:] = depth

    var = ds.createVariable(var_name, 'f4', ('depth','time'))
    var[:,:] = var_grid

    ds.min_lon = min_lon
    ds.min_lat = min_lat
    ds.max_lon = max_lon
    ds.max_lat = max_lat

    ds.close()

folder = '/Users/mhwood/Documents/Research/Projects/Disko Bay'
config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/' \
             'configurations/downscale_darwin'
model_name = 'L1_W_Greenland'
results_dir = os.path.join(config_dir,'L1',model_name,'results','BGC_consts_mon_mean')

folder = '/Users/mhwood/Documents/Research/Projects/Disko Bay'

# just the trough in Disko
min_lat = 68.8750
min_lon = -53.3744
max_lat = 68.9840
max_lon = -53.1117

# all of Disko
min_lat = 68.5252
min_lon = -53.8670
max_lat = 69.8216
max_lon = -50.5690

var_name = 'Phosphate'

depth, years, months, days, dec_yrs, var_grid = find_profile_files_in_box(results_dir, min_lon, min_lat, max_lon, max_lat, var_name)

outputting_transects_to_nc(folder, var_name, years, months, days, dec_yrs,
                               depth, min_lon, min_lat, max_lon, max_lat, var_grid)


