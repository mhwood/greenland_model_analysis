
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime as dt

def read_model_grid(config_folder, model_name):

    ds = nc4.Dataset(os.path.join(config_folder,'nc_grids',model_name+'_grid.nc'))
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    drF = ds.variables['drF'][:]
    ds.close()

    return(XC, YC, drF)

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = dt.datetime(year,month,day,hour,minute,second)
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_timeseries_from_monthly_nc(model_folder,var_name,drF, XC, YC, x_index, y_index):

    if var_name=='Theta' or var_name=='Salt':
        subset_name = 'state_3D_mon_mean'
    if var_name=='Total_Chl':
        subset_name = 'BGC_Chl_mon_mean'
    if var_name=='NO3' or var_name=='PO4':
        subset_name = 'BGC_consts_mon_mean'

    data_folder = os.path.join(model_folder,'results',subset_name)

    file_names = []
    for file_name in os.listdir(data_folder):
        if file_name[0]!='.' and file_name[-3:]=='.nc':
            file_names.append(file_name)
    file_names = sorted(file_names)

    counter = 0
    dec_yrs = []
    profiles = np.zeros((len(drF),len(file_names)))
    for file_name in file_names:

        print('        - Reading from '+file_name)

        year = int(file_name.split('.')[1][:4])
        month = int(file_name.split('.')[1][4:6])
        dec_yr = YMD_to_DecYr(year,month,15)
        dec_yrs.append(dec_yr)

        ds = nc4.Dataset(os.path.join(data_folder,file_name))
        depths = ds.variables['depths'][:]
        # lon = ds.variables['longitude'][:, :]
        # lat = ds.variables['latitude'][:, :]
        if var_name=='Total_Chl':
            var_grid = ds.variables['Chl01'][:, :, :, :]
            var_grid += ds.variables['Chl02'][:, :, :, :]
            var_grid += ds.variables['Chl03'][:, :, :, :]
            var_grid += ds.variables['Chl04'][:, :, :, :]
            var_grid += ds.variables['Chl05'][:, :, :, :]
        else:
            var_grid = ds.variables[var_name][:, :, :, :]
        ds.close()

        profiles[:, counter] = var_grid[0,:,y_index, x_index]

        counter +=1

    depths = np.array(depths)
    dec_yrs = np.array(dec_yrs)
    profiles = profiles[:,:len(dec_yrs)]

    depth_indices = profiles[:,0]!=0
    profiles = profiles[depth_indices,:]
    depths = depths[depth_indices]

    # C = plt.imshow(profiles)
    # plt.colorbar(C)
    # plt.show()

    return(dec_yrs, depths, profiles)

def write_profile_timeseries_to_nc(project_folder, model_name, dec_yrs, depths, profiles, var_name):

    ds = nc4.Dataset(os.path.join(project_folder, 'Data', 'Ocean', 'L1', model_name+'_'+var_name+'_timeseries.nc'),'w')

    ds.createDimension('time',len(dec_yrs))
    ds.createDimension('depth',len(depths))

    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = dec_yrs

    dvar = ds.createVariable('depth','f4',('depth',))
    dvar[:] = depths

    vvar = ds.createVariable(var_name,'f4',('depth','time'))
    vvar[:, :] = profiles

    ds.close()
    a=1

def create_model_profile_timeseres(project_folder, config_folder, var_name):

    model_name = 'L1_W_Greenland'

    XC, YC, drF = read_model_grid(config_folder, model_name)

    latitude = 69.209
    longitude = -52.147
    dist = (XC-longitude)**2 + (YC-latitude)**2
    y_index, x_index = np.where(dist==np.min(dist))
    x_index = x_index[0]
    y_index = y_index[0]
    # print(x_index, y_index, XC[y_index, x_index], YC[y_index, x_index])

    model_folder = os.path.join(config_folder,'L1',model_name)
    dec_yrs, depths, profiles = read_timeseries_from_monthly_nc(model_folder, var_name, drF, XC, YC, x_index, y_index)

    write_profile_timeseries_to_nc(project_folder, model_name, dec_yrs, depths, profiles, var_name)


project_folder = '/Users/mhwood/Documents/Research/Projects/Disko Bay'


config_folder = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/' \
               'darwin3/configurations/downscale_darwin'

var_name = 'Total_Chl'

create_model_profile_timeseres(project_folder, config_folder, var_name)




