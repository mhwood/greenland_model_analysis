
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4

folder = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/Data/Modeling/' \
         'ECCO/'

sample_area = 'IPCC'
sample_area = 'Scoresby_Sund'

if sample_area=='Scoresby_Sund':
    ds =nc4.Dataset(folder+'/ECCOv5_THETA_'+sample_area+'_v5_timeseries.nc')
    time = ds.variables['time'][:]
    theta = ds.groups['sound_entrance'].variables['THETA'][:]
    ds.close()
    timeseries_1 = np.column_stack([time,theta])
    timeseries_1=timeseries_1[:12*(2018-1992),:]
else:
    ds = nc4.Dataset(folder + '/ECCOv5_THETA_' + sample_area + '_v5_timeseries.nc')
    time = ds.variables['dec_yr'][:]
    theta = ds.variables['THETA_mean'][:]
    ds.close()
    timeseries_1 = np.column_stack([time, theta])

ds =nc4.Dataset(folder+'/ECCOv5_THETA_'+sample_area+'_extension_timeseries.nc')
time = ds.variables['dec_yr'][:]
theta = ds.variables['THETA_mean'][:]
ds.close()
timeseries_2 = np.column_stack([time,theta])

joint_timeseries = np.zeros((np.shape(timeseries_1)[0]+np.shape(timeseries_2)[0]-6,2))
joint_timeseries[:np.shape(timeseries_1)[0]-6,:]=timeseries_1[:-6,:]
joint_timeseries[np.shape(timeseries_1)[0]-6:,:]=timeseries_2

output_timeseries = np.copy(joint_timeseries)

smoothing_radius = 2
for i in range(np.shape(timeseries_1)[0]-6-smoothing_radius,np.shape(timeseries_1)[0]-6+smoothing_radius+1):
    output_timeseries[i,1] = np.mean(joint_timeseries[i-smoothing_radius:i+smoothing_radius,1])


# fill in some averaged values
# joint_timeseries[np.shape(timeseries_1)[0]-5,1] = (timeseries_1[-5,1]+timeseries_2[0,1])/2
# joint_timeseries[np.shape(timeseries_1)[0]-4,1] = (timeseries_1[-4,1]+timeseries_2[1,1])/2

plt.plot(timeseries_1[:,0], timeseries_1[:,1],linewidth=4, alpha=0.5)
plt.plot(timeseries_2[:,0], timeseries_2[:,1],linewidth=4, alpha=0.5)
plt.plot(output_timeseries[:,0], output_timeseries[:,1], 'k-')
plt.show()

output_file = folder+'/ECCOv5_THETA_'+sample_area+'_timeseries.nc'
ds = nc4.Dataset(output_file,'w')

ds.createDimension('time',np.shape(output_timeseries)[0])

tvar = ds.createVariable('dec_yr','f4',('time',))
tvar[:] = output_timeseries[:,0]
mvar = ds.createVariable('THETA_mean', 'f4', ('time',))
mvar[:] = output_timeseries[:, 1]

ds.close()

