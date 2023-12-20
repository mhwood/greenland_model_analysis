
import os
import numpy as np
import netCDF4 as nc4


def create_mean_fields(ecco_dir,field_name):

    for year in range(1992,2018):
        print('Reading in data for year '+str(year))
        file_name = os.path.join(ecco_dir,field_name,field_name+'_'+str(year)+'.nc')
        ds = nc4.Dataset(file_name)
        if year==1992:
            i = ds.variables['i'][:]
            j = ds.variables['j'][:]
            k = ds.variables['k'][:]
            var = ds.variables[field_name][:, :, :, :, :]
            counter = np.shape(var)[0]
            sum_var = np.sum(var,axis=0)
            tile = ds.variables['tile'][:]
        else:
            var = ds.variables[field_name][:, :, :, :, :]
            counter += np.shape(var)[0]
            sum_var += np.sum(var, axis=0)
        ds.close()

    mean_var = sum_var/counter

    file_name = os.path.join(ecco_dir, 'Mean_Fields', field_name + '_Mean.nc')
    ds = nc4.Dataset(file_name,'w')

    ds.createDimension('tile', 13)
    ds.createDimension('i', 270)
    ds.createDimension('j', 270)
    ds.createDimension('k', 50)

    ivar = ds.createVariable('i','f4',('i',))
    ivar[:] = i
    jvar = ds.createVariable('j', 'f4', ('j',))
    jvar[:] = j
    kvar = ds.createVariable('k', 'f4', ('k',))
    kvar[:] = k
    tvar = ds.createVariable('tile', 'f4', ('tile',))
    tvar[:] = tile
    fvar = ds.createVariable(field_name,'f4',('k','tile','j','i'))
    fvar[:, :, :, :] = mean_var

    ds.close()



ecco_dir = '/Volumes/ifenty/Wood/Data Repository/ECCO'
field_name ='THETA'

create_mean_fields(ecco_dir,field_name)

