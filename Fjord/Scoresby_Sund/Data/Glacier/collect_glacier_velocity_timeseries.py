
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4

def get_glacier_velocity_timeseries(data_folder,glacier,source='ITS_LIVE Annual Mosaic Velocity Grids'):

    file_path = os.path.join(data_folder,'Data','Glaciers',glacier,'Velocity','Timeseries','Sample Point','Sources',
                             glacier+' '+source+' Sample Point Velocity Timeseries.csv')
    timeseries = np.genfromtxt(file_path,delimiter=',')
    print(np.shape(timeseries))

    return(timeseries)

def save_timeseries_as_nc(output_file,its_live_velocity_timeseries,measures_velocity_timeseries):

    ds = nc4.Dataset(output_file,'w')

    grp = ds.createGroup('ITS_LIVE')
    grp.createDimension('time',np.shape(its_live_velocity_timeseries)[0])
    tvar = grp.createVariable('time','f4',('time',))
    tvar[:] = its_live_velocity_timeseries[:,0]
    vvar = grp.createVariable('velocity', 'f4', ('time',))
    vvar[:] = its_live_velocity_timeseries[:, 1]
    tevar = grp.createVariable('time_span', 'f4', ('time',))
    tevar[:] = its_live_velocity_timeseries[:, 2]
    vevar = grp.createVariable('velocity_error', 'f4', ('time',))
    vevar[:] = its_live_velocity_timeseries[:, 3]

    grp = ds.createGroup('MEaSUREs')
    grp.createDimension('time', np.shape(measures_velocity_timeseries)[0])
    tvar = grp.createVariable('time', 'f4', ('time',))
    tvar[:] = measures_velocity_timeseries[:, 0]
    vvar = grp.createVariable('velocity', 'f4', ('time',))
    vvar[:] = measures_velocity_timeseries[:, 1]
    tevar = grp.createVariable('time_span', 'f4', ('time',))
    tevar[:] = measures_velocity_timeseries[:, 2]
    vevar = grp.createVariable('velocity_error', 'f4', ('time',))
    vevar[:] = measures_velocity_timeseries[:, 3]

    ds.close()



data_folder = '/Volumes/upernavik/Research/Projects/Glacier Variability'

glacier = 'Daugaard Jensen'

its_live_velocity_timeseries = get_glacier_velocity_timeseries(data_folder,glacier,source='ITS_LIVE Annual Mosaic Velocity Grids')
measures_velocity_timeseries = get_glacier_velocity_timeseries(data_folder,glacier,source='MEaSUREs Annual Mosaic Velocity Grids')


output_file = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/Data/Remote Sensing/Velocity/' \
              'Daugaard Jensen Velocity Timeseries.nc'
save_timeseries_as_nc(output_file,its_live_velocity_timeseries,measures_velocity_timeseries)





