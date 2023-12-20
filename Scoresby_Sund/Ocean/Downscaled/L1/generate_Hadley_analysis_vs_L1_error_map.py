
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import netCDF4 as nc4
from scipy.interpolate import interp1d
import datetime

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_model_grid():
    grid_file = '/Volumes/mhwood/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
                'downscaled_greenland/nc_grids/L1_CE_Greenland_grid.nc'
    ds = nc4.Dataset(grid_file)
    xc = ds.variables['XC'][:, :]
    yc = ds.variables['YC'][:, :]
    drF = ds.variables['drF'][:]
    Depth = ds.variables['Depth'][:,:]
    ds.close()

    Z_bottom = np.cumsum(drF)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2
    return(xc,yc,Z,Depth)

def read_indices_of_hadley_analysis_covering_domain(XC, YC, hadley_analysis_dir):

    ds = nc4.Dataset(os.path.join(hadley_analysis_dir,'EN.4.2.2.analyses.g10.1993',
                                  'EN.4.2.2.f.analysis.g10.199301.nc'))
    lon = ds.variables['lon'][:]
    lat = ds.variables['lat'][:]
    ds.close()

    Lon, Lat = np.meshgrid(lon,lat)

    Lon = np.hstack([Lon[:,-180:],Lon[:,:180]])
    Lon[Lon>180]-=360
    Lat = np.hstack([Lat[:, -180:], Lat[:, :180]])

    # C = plt.imshow(Lon,origin='lower')
    # plt.colorbar(C)
    # plt.show()

    ll_dist = ((Lon - XC[0,0]) ** 2 + (Lat - YC[0,0]) ** 2)
    ll_row, ll_col = np.where(ll_dist == np.min(ll_dist))

    lr_dist = ((Lon - XC[0,-1]) ** 2 + (Lat - YC[0,-1]) ** 2)
    lr_row, lr_col = np.where(lr_dist == np.min(lr_dist))

    ur_dist = ((Lon - XC[-1,-1]) ** 2 + (Lat - YC[-1,-1]) ** 2)
    ur_row, ur_col = np.where(ur_dist == np.min(ur_dist))

    ul_dist = ((Lon - XC[-1,0]) ** 2 + (Lat - YC[-1,0]) ** 2)
    ul_row, ul_col = np.where(ul_dist == np.min(ul_dist))

    min_row = np.min([ll_row[0], lr_row[0], ur_row[0], ul_row[0]])-1
    max_row = np.max([ll_row[0], lr_row[0], ur_row[0], ul_row[0]])+1
    min_col = np.min([ll_col[0], lr_col[0], ur_col[0], ul_col[0]])-1
    max_col = np.max([ll_col[0], lr_col[0], ur_col[0], ul_col[0]])+1

    Lon_subset = Lon[min_row:max_row+1,min_col:max_col+1]
    Lat_subset = Lat[min_row:max_row+1,min_col:max_col+1]

    # plt.subplot(1,2,1)
    # C = plt.imshow(Lon_subset, origin='lower')
    # plt.colorbar(C)
    # plt.subplot(1, 2, 2)
    # C = plt.imshow(Lat_subset, origin='lower')
    # plt.colorbar(C)
    # plt.show()

    return(min_row,max_row,min_col,max_col, Lon_subset, Lat_subset)

def get_year_hadley_timeseries(hadley_analysis_dir, year,
                               min_row, max_row, min_col, max_col, min_depth, max_depth):

    intp_depths = np.arange(min_depth, max_depth, +1).astype(float)

    year_T_timseries = np.zeros((12, max_row - min_row + 1, max_col - min_col + 1))
    year_S_timseries = np.zeros((12, max_row - min_row + 1, max_col - min_col + 1))

    for month in range(1,13):
        file_path = os.path.join(hadley_analysis_dir,'EN.4.2.2.analyses.g10.'+str(year),
                                 'EN.4.2.2.f.analysis.g10.'+str(year)+'{:02d}'.format(month)+'.nc')
        ds = nc4.Dataset(file_path)
        temp = ds.variables['temperature'][:,:,:,:]
        salt = ds.variables['salinity'][:, :, :, :]
        depth = ds.variables['depth'][:]
        ds.close()

        temp -= 273.15

        max_depth_index = np.argmin(np.abs(depth-max_depth))+1

        temp = np.concatenate([temp[0, :, : ,-180:],temp[0,:,:,:180]], axis = 2)
        salt = np.concatenate([salt[0, :, :, -180:], salt[0,:,:, :180]], axis = 2)

        for row in range(min_row,max_row+1):
            for col in range(min_col,max_col+1):

                T_profile = temp[:max_depth_index,row,col].astype(float)
                S_profile = salt[:max_depth_index, row, col].astype(float)

                if np.sum(np.isnan(T_profile))==0 and np.sum(T_profile<-2)==0:
                    set_int = interp1d(depth[:max_depth_index].astype(float),T_profile)
                    year_T_timseries[month-1,row-min_row,col-min_col] = np.mean(set_int(intp_depths))
                    set_int = interp1d(depth[:max_depth_index].astype(float), S_profile)
                    year_S_timseries[month - 1, row-min_row, col-min_col] = np.mean(set_int(intp_depths))

    # plt.subplot(1,2,1)
    # C = plt.imshow(year_T_timseries[2,:,:], origin='lower')
    # plt.colorbar(C)
    # plt.subplot(1, 2, 2)
    # C = plt.imshow(year_S_timseries[2,:,:], origin='lower')
    # plt.colorbar(C)
    # plt.show()

    time = np.zeros((12,))
    for month in range(1,13):
        time[month-1] = YMD_to_DecYr(year,month,15)

    return(time, year_T_timseries, year_S_timseries)

def get_year_model_timeseries(state_mon_dir,year, XC, YC, Lon, Lat, min_depth, max_depth):

    intp_depths = np.arange(min_depth, max_depth, +1).astype(float)

    # for yr in range(year-1,year+2):
    file_prefix = 'L1_state_3D_mon_mean.' + str(year)
    for file_name in os.listdir(state_mon_dir):
        if file_name[:len(file_prefix)]==file_prefix:
            ds = nc4.Dataset(os.path.join(state_mon_dir,file_name))
            time = ds.variables['time'][:]
            depth = ds.variables['depth'][:]
            theta = ds.variables['Theta'][:,:,:,:]
            salt = ds.variables['Salt'][:, :, :, :]
            ds.close()

    time *= 1 / (24 * 60 * 60 * 365.25)
    time += 1992

    top = np.column_stack([XC[-1, :], YC[-1,:]])
    bottom = np.column_stack([XC[0, :], YC[0, :]])
    left = np.column_stack([XC[:,0], YC[:,0]])
    right = np.column_stack([XC[:, -1], YC[:, -1]])
    box = np.vstack([bottom,right,np.flipud(top),np.flipud(left)])
    p = mplPath.Path(box)

    year_T_timseries = np.zeros((len(time), np.shape(Lon)[0], np.shape(Lon)[1]))
    year_S_timseries = np.zeros((len(time), np.shape(Lon)[0], np.shape(Lon)[1]))

    for row in range(np.shape(Lon)[0]):
        if row%5==0:
            print('             - Working on row '+str(row)+' to '+str(np.min([row+5,np.shape(Lon)[0]]))+' of '+str(np.shape(Lon)[0]))
        for col in range(np.shape(Lat)[1]):
            if p.contains_point((Lon[row,col],Lat[row,col])):
                lon_indices = np.logical_and(XC>Lon[row,col]-0.5, XC<Lon[row,col]+0.5)
                lat_indices = np.logical_and(YC > Lat[row, col] - 0.5, YC < Lat[row, col] + 0.5)
                depth_indices = Depth>max_depth
                all_indices = np.logical_and(np.logical_and(lon_indices,lat_indices),depth_indices)

                if np.sum(all_indices)>0:
                    # fill in the profiles with spatial means
                    T_profile = np.zeros((np.shape(theta)[0],np.shape(theta)[1]))
                    S_profile = np.zeros((np.shape(theta)[0],np.shape(theta)[1]))
                    for t in range(np.shape(theta)[0]):
                        for d in range(np.shape(theta)[1]):
                            theta_slice = theta[t,d,:,:]
                            T_profile[t,d] = np.mean(theta_slice[all_indices])
                            salt_slice = salt[t,d,:,:]
                            S_profile[t, d] = np.mean(salt_slice[all_indices])

                    # average over depths
                    for t in range(np.shape(theta)[1]):
                        set_int = interp1d(depth, T_profile[:,t])
                        year_T_timseries[t, row, col] = np.mean(set_int(intp_depths))
                        set_int = interp1d(depth, S_profile[:,t])
                        year_S_timseries[t, row, col] = np.mean(set_int(intp_depths))

    # plt.subplot(1,2,1)
    # C = plt.imshow(year_T_timseries[2,:,:], origin='lower')
    # plt.colorbar(C)
    # plt.subplot(1, 2, 2)
    # C = plt.imshow(year_S_timseries[2,:,:], origin='lower')
    # plt.colorbar(C)
    # plt.show()

    return(time,year_T_timseries,year_S_timseries)


def save_year_timeseries_as_nc(project_folder, year, min_depth, max_depth,
                               Lon, Lat,
                               hadley_time, hadley_year_T, hadley_year_S,
                               model_time, model_year_T, model_year_S):

    output_file = os.path.join(project_folder,'Data', 'Modeling', 'Downscaled', 'L1_CE_Greenland', 'Hadley Analysis Comparison',
                               'L1_CE_Greenland_vs_Hadley_Analysis_' + str(year) + '_' + str(min_depth) + 'm_' + str(max_depth) + 'm.nc')

    ds = nc4.Dataset(output_file, 'w')

    ds.createDimension('lon',np.shape(Lon)[1])
    ds.createDimension('lat', np.shape(Lon)[0])

    lvar = ds.createVariable('Lon','f4',('lat','lon'))
    lvar[:,:] = Lon

    lvar = ds.createVariable('Lat', 'f4', ('lat', 'lon'))
    lvar[:, :] = Lat

    grp = ds.createGroup('Hadley')

    grp.createDimension('time', len(hadley_time))

    tvar = grp.createVariable('time', 'f4', ('time',))
    tvar[:] = hadley_time

    Tvar = grp.createVariable('Theta', 'f4', ('time','lat','lon'))
    Tvar[:, :, :] = hadley_year_T

    Svar = grp.createVariable('Salt', 'f4', ('time', 'lat', 'lon'))
    Svar[:, :, :] = hadley_year_S

    grp = ds.createGroup('Model')

    grp.createDimension('time', len(model_time))

    tvar = grp.createVariable('time', 'f4', ('time',))
    tvar[:] = model_time

    Tvar = grp.createVariable('Theta', 'f4', ('time', 'lat', 'lon'))
    Tvar[:, :, :] = model_year_T

    Svar = grp.createVariable('Salt', 'f4', ('time', 'lat', 'lon'))
    Svar[:, :, :] = model_year_S

    ds.close()

def calculate_errors_from_year_files(project_folder,years,min_depth, max_depth):

    stack_started = False
    year_dir = os.path.join(project_folder,'Data', 'Modeling', 'Downscaled',
                                                'L1_CE_Greenland', 'Hadley Analysis Comparison')

    for year in years:
        year_file = 'L1_CE_Greenland_vs_Hadley_Analysis_' + str(year) + '_' + str(min_depth) + 'm_' + str(max_depth) + 'm.nc'
        ds = nc4.Dataset(os.path.join(year_dir,year_file))
        hadley_year_time = ds.groups['Hadley'].variables['time'][:]
        hadley_year_T = ds.groups['Hadley'].variables['Theta'][:,:,:]
        hadley_year_S = ds.groups['Hadley'].variables['Salt'][:, :, :]
        model_year_time = ds.groups['Model'].variables['time'][:]
        model_year_T = ds.groups['Model'].variables['Theta'][:, :, :]
        model_year_S = ds.groups['Model'].variables['Salt'][:, :, :]
        ds.close()

        if not stack_started:
            stack_started = True
            hadley_time = np.reshape(hadley_year_time,(len(hadley_year_time),1))
            hadley_T = hadley_year_T
            hadley_S = hadley_year_S
            model_time = np.reshape(model_year_time, (len(model_year_time), 1))
            model_T = model_year_T
            model_S = model_year_S
        else:
            hadley_time = np.vstack([hadley_time,
                                     np.reshape(hadley_year_time, (len(hadley_year_time), 1))])
            hadley_T = np.concatenate([hadley_T,hadley_year_T],axis=0)
            hadley_S = np.concatenate([hadley_S,hadley_year_S],axis=0)
            model_time = np.vstack([model_time,
                                    np.reshape(model_year_time, (len(model_year_time), 1))])
            model_T = np.concatenate([model_T,model_year_T],axis=0)
            model_S = np.concatenate([model_S,model_year_S],axis=0)

    error_T = np.zeros_like(hadley_T)
    error_S = np.zeros_like(hadley_S)

    for t in range(np.shape(hadley_T)[0]):
        model_time_index = np.argmin(np.abs(model_time-hadley_time[t]))
        error_T[t,:,:] = model_T[model_time_index,:,:] - hadley_T[t,:,:]
        error_S[t, :, :] = model_S[model_time_index, :, :] - hadley_S[t, :, :]

    error_time = hadley_time

    return(hadley_time, hadley_T, hadley_S,
           model_time, model_T, model_S,
           error_time, error_T, error_S)

def save_error_timeseries_as_nc(project_folder, min_depth, max_depth,
                                Lon, Lat,
                                hadley_time, hadley_T, hadley_S,
                                model_time, model_T, model_S,
                                error_time, error_T, error_S):

    output_file = os.path.join(project_folder,'Data', 'Modeling', 'Downscaled', 'L1_CE_Greenland',
                               'L1_CE_Greenland_vs_Hadley_Analysis_' + str(min_depth) + 'm_' + str(max_depth) + 'm.nc')

    ds = nc4.Dataset(output_file, 'w')

    ds.createDimension('lon',np.shape(Lon)[1])
    ds.createDimension('lat', np.shape(Lon)[0])

    lvar = ds.createVariable('Lon','f4',('lat','lon'))
    lvar[:,:] = Lon

    lvar = ds.createVariable('Lat', 'f4', ('lat', 'lon'))
    lvar[:, :] = Lat

    grp = ds.createGroup('Hadley')

    grp.createDimension('time', len(hadley_time))

    tvar = grp.createVariable('time', 'f4', ('time',))
    tvar[:] = hadley_time

    Tvar = grp.createVariable('Theta', 'f4', ('time','lat','lon'))
    Tvar[:, :, :] = hadley_T

    Svar = grp.createVariable('Salt', 'f4', ('time', 'lat', 'lon'))
    Svar[:, :, :] = hadley_S

    grp = ds.createGroup('Model')

    grp.createDimension('time', len(model_time))

    tvar = grp.createVariable('time', 'f4', ('time',))
    tvar[:] = model_time

    Tvar = grp.createVariable('Theta', 'f4', ('time', 'lat', 'lon'))
    Tvar[:, :, :] = model_T

    Svar = grp.createVariable('Salt', 'f4', ('time', 'lat', 'lon'))
    Svar[:, :, :] = model_S

    grp = ds.createGroup('Error')

    grp.createDimension('time', len(error_time))

    tvar = grp.createVariable('time', 'f4', ('time',))
    tvar[:] = error_time

    Tvar = grp.createVariable('Theta', 'f4', ('time', 'lat', 'lon'))
    Tvar[:, :, :] = error_T

    Svar = grp.createVariable('Salt', 'f4', ('time', 'lat', 'lon'))
    Svar[:, :, :] = error_S

    ds.close()




project_folder = '/Users/michwood/Documents/Research/Projects/Scoresby Sund'

state_mon_dir = '/Volumes/mhwood/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
                'downscaled_greenland/L1/grid/L1_CE_Greenland/results/monthly_nc/state_3D_mon_mean'

hadley_analysis_dir = '/Volumes/mhwood/Data_Repository/Hadley Analysis'

years = np.arange(1992,2022)
min_depth = 200
max_depth = 500

# step 0: get the model grid
XC, YC, Z, Depth = read_model_grid()

# find the indices of the Hadley analysis in the domain
min_row, max_row, min_col, max_col, Lon, Lat = \
    read_indices_of_hadley_analysis_covering_domain(XC, YC, hadley_analysis_dir)

for year in years:

    year_file = 'L1_CE_Greenland_vs_Hadley_Analysis_' + str(year) + '_' + str(min_depth) + 'm_' + str(max_depth) + 'm.nc'

    if year_file not in os.listdir(os.path.join(project_folder,'Data', 'Modeling', 'Downscaled',
                                                'L1_CE_Greenland', 'Hadley Analysis Comparison')):

        print('    - Working on the comparison in year '+str(year))

        print('        - Gathering the Hadley analysis data')
        # for each cell, fill in the mean temperature (200-500 m) from the objective analysis
        hadley_time, hadley_year_T, hadley_year_S = get_year_hadley_timeseries(hadley_analysis_dir, year,
                                                                              min_row, max_row, min_col, max_col,
                                                                              min_depth, max_depth)

        print('        - Gathering the model data')
        # for each cell, fill in the mean temperature (200-500 m) from the downscaled model
        model_time, model_year_T, model_year_S = get_year_model_timeseries(state_mon_dir, year, XC, YC, Lon, Lat, min_depth, max_depth)

        print('        - Saving data for this year as an nc')
        # save both T timeseries and S timeseries
        save_year_timeseries_as_nc(project_folder, year, min_depth, max_depth,
                                   Lon, Lat,
                                   hadley_time, hadley_year_T, hadley_year_S,
                                   model_time, model_year_T, model_year_S)

    else:
        print('    - '+year_file+' already created')

hadley_time, hadley_T, hadley_S, model_time, model_T, model_S, error_time, error_T, error_S = \
    calculate_errors_from_year_files(project_folder,years,min_depth, max_depth)

save_error_timeseries_as_nc(project_folder, min_depth, max_depth,
                                Lon, Lat,
                                hadley_time, hadley_T, hadley_S,
                                model_time, model_T, model_S,
                                error_time, error_T, error_S)
