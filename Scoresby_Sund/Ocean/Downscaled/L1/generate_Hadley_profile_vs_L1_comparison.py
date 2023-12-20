
import os
import numpy as np
import matplotlib.pyplot as plt
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
    ds.close()

    Z_bottom = np.cumsum(drF)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2
    return(xc,yc,Z)

def get_year_model_timeseries(state_mon_dir,year):

    # for yr in range(year-1,year+2):
    file_prefix = 'L1_state_3D_mon_mean.' + str(year)
    for file_name in os.listdir(state_mon_dir):
        if file_name[:len(file_prefix)]==file_prefix:
            ds = nc4.Dataset(os.path.join(state_mon_dir,file_name))
            time = ds.variables['time'][:]
            theta = ds.variables['Theta'][:,:,:,:]
            ds.close()

                # for t in range(len(time)):
                #     dec_yr = 1992+time[t]/(24*60*60*365.25)
                #     min_dec_yr = dec_yr - 2/(365.25)
                #     max_dec_yr = dec_yr + 28 / (365.25)
                #     print(min_dec_yr, dec_yr, max_dec_yr)
                #
                # # mean_etan = np.zeros((len(time),1))
                # # time = np.reshape(time,(len(time),1))
                # # for i in range(len(time)):
                # #     subset = etan[i,:,:]
                # #     mean_etan[i] = np.mean(subset[subset!=0])

    time *= 1/(24*60*60*365.25)
    time += 1992

    return(time,theta)

def get_observation_timeseries(project_folder, year, min_depth, max_depth):

    obs_dir = os.path.join(project_folder,'Data','In Situ','Hadley','Data',str(year))

    xytT_obs = np.zeros((len(os.listdir(obs_dir)),4))

    interp_depths = np.arange(min_depth,max_depth)

    counter = 0
    for file_name in os.listdir(obs_dir):
        if file_name[0]!='.' and file_name[-3:]=='.nc':
            ds= nc4.Dataset(os.path.join(obs_dir,file_name))
            depth = ds.variables['depth'][:]
            if len(depth)>0:
                if np.min(depth)<min_depth and np.max(depth)>max_depth:
                    xytT_obs[counter,0] = ds.longitude
                    xytT_obs[counter,1] = ds.latitude
                    try:
                        dec_yr = YMD_to_DecYr(ds.year, ds.month, ds.day, ds.hour, ds.minute)
                        xytT_obs[counter,2] = dec_yr
                        set_int = interp1d(depth,ds.variables['potential_temperature'][:])
                        theta = np.mean(set_int(interp_depths))
                        xytT_obs[counter,3] = theta
                        counter += 1
                    except:
                        print('Error found in '+file_name+' ('+str(ds.year)+'/'+str(ds.month)+'/'+str(ds.day)+')')
            ds.close()

    print('        - Found '+str(counter)+' profiles in year '+str(year))

    xytT_obs = xytT_obs[xytT_obs[:,0]!=0,:]

    return(xytT_obs)

def great_circle_distance(lon_ref, lat_ref, Lon, Lat):
    earth_radius = 6371000
    lon_ref_radians = np.radians(lon_ref)
    lat_ref_radians = np.radians(lat_ref)
    lons_radians = np.radians(Lon)
    lats_radians = np.radians(Lat)
    lat_diff = lats_radians - lat_ref_radians
    lon_diff = lons_radians - lon_ref_radians
    d = np.sin(lat_diff * 0.5) ** 2 + np.cos(lat_ref_radians) * np.cos(lats_radians) * np.sin(lon_diff * 0.5) ** 2
    h = 2 * earth_radius * np.arcsin(np.sqrt(d))
    return(h)

def calculate_model_data_differences(XC, YC, Z, model_decyr, model_theta, xytT_obs, min_depth, max_depth):

    indices = model_theta[0,0,:,:]!=0

    wet_X = XC[indices]
    wet_Y = YC[indices]

    interp_depths = np.arange(min_depth, max_depth)

    # model x, model y, model time, model theta,
    # obs x, obs y, obs time, obs, theta,
    # separation_distance, time_difference (model - obs), theta_difference (model-obs)
    model_data_differences = np.zeros((np.shape(xytT_obs)[0],11))

    for i in range(np.shape(xytT_obs)[0]):
        # for now, just get the closest point at the closest time and interpolate the model depths
        x = xytT_obs[i,0]
        y = xytT_obs[i,1]

        # calculate the distance stats
        dist = great_circle_distance(x,y,wet_X,wet_Y)
        dist_index = np.where(dist==np.min(dist))[0][0]
        row, col = np.where(np.logical_and(XC==wet_X[dist_index], YC==wet_Y[dist_index]))

        # store the distance stats
        model_data_differences[i, 0] = XC[row[0],col[0]]
        model_data_differences[i, 1] = YC[row[0], col[0]]
        model_data_differences[i, 4] = x
        model_data_differences[i, 5] = y
        model_data_differences[i, 8] = dist[dist_index]

        # calculate the time stats
        time_index = np.argmin(np.abs(model_decyr - xytT_obs[i, 2]))

        # store the time stats
        model_data_differences[i, 2] = model_decyr[time_index]
        model_data_differences[i, 6] = xytT_obs[i, 2]
        model_data_differences[i, 9] = model_decyr[time_index] - xytT_obs[i, 2]

        # calculate the theta stats
        model_profile = model_theta[:,time_index,row[0],col[0]]
        set_int = interp1d(Z,model_profile)
        theta = np.mean(set_int(interp_depths))

        model_data_differences[i, 3] = theta
        model_data_differences[i, 7] = xytT_obs[i, 3]
        model_data_differences[i, 10] = theta - xytT_obs[i, 3]

    # plt.plot(model_data_differences[:,2],model_data_differences[:,3],'k.')
    #     # plt.plot([-1,7],[-1,7])
    #     # plt.xlabel('Model')
    #     # plt.ylabel('Data')
    #     # plt.show()

    return(model_data_differences)

def save_model_data_differences_as_nc(project_folder, model_data_differences, year, min_depth, max_depth):

    output_file = os.path.join(project_folder,'Data','In Situ','Hadley','Model_Comparison','L1_CE_Greenland',
                               'L1_CE_Model_Data_Differences_'+str(year)+'_'+str(min_depth)+'m_'+str(max_depth)+'m.nc')

    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('n_samples',np.shape(model_data_differences)[0])

    var = ds.createVariable('model_x','f4',('n_samples',))
    var[:] = model_data_differences[:, 0]

    var = ds.createVariable('model_y', 'f4', ('n_samples',))
    var[:] = model_data_differences[:, 1]

    var = ds.createVariable('model_decyr', 'f4', ('n_samples',))
    var[:] = model_data_differences[:, 2]

    var = ds.createVariable('model_theta', 'f4', ('n_samples',))
    var[:] = model_data_differences[:, 3]

    var = ds.createVariable('obs_x', 'f4', ('n_samples',))
    var[:] = model_data_differences[:, 4]

    var = ds.createVariable('obs_y', 'f4', ('n_samples',))
    var[:] = model_data_differences[:, 5]

    var = ds.createVariable('obs_decyr', 'f4', ('n_samples',))
    var[:] = model_data_differences[:, 6]

    var = ds.createVariable('obs_theta', 'f4', ('n_samples',))
    var[:] = model_data_differences[:, 7]

    var = ds.createVariable('separation_dist', 'f4', ('n_samples',))
    var[:] = model_data_differences[:, 8]

    var = ds.createVariable('time_difference', 'f4', ('n_samples',))
    var[:] = model_data_differences[:, 9]
    var.note = 'model - obs'

    var = ds.createVariable('theta_difference', 'f4', ('n_samples',))
    var[:] = model_data_differences[:, 10]
    var.note = 'model - obs'

    ds.close()

project_folder = '/Users/michwood/Documents/Research/Projects/Scoresby Sund'

state_mon_dir = '/Volumes/mhwood/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
            'downscaled_greenland/L1/grid/L1_CE_Greenland/results/monthly_nc/state_3D_mon_mean'

min_depth = 200
max_depth = 500

# step 0: get the model grid
XC, YC, Z = read_model_grid()

for year in range(2013,2022):

    # step 1: get the theta timeseries
    model_decyr, model_theta = get_year_model_timeseries(state_mon_dir,year)

    # step 2: get the list of observations for this year (which cover the depth spans)
    xytT_obs = get_observation_timeseries(project_folder, year, min_depth, max_depth)

    # step 3: calculate model-data differences over the depth range
    model_data_differences = calculate_model_data_differences(XC, YC, Z, model_decyr, model_theta, xytT_obs, min_depth, max_depth)

    # step 4: output the differences as an annual nc
    save_model_data_differences_as_nc(project_folder, model_data_differences, year, min_depth, max_depth)
