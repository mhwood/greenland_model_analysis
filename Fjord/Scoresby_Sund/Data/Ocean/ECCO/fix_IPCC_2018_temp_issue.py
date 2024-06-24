

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import CubicSpline



folder = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/Data/Modeling/' \
         'ECCO'

ds = nc4.Dataset(folder+'/ECCOv5_THETA_IPCC_timeseries_2018_issue.nc')
time = ds.variables['dec_yr'][:]
theta = ds.variables['THETA_mean'][:]
ds.close()

indices_2018 = np.floor(time).astype(int)==2018
time_subset = time[~indices_2018]
theta_subset = theta[~indices_2018]

# plt.plot(time_subset, theta_subset, 'k.')
# plt.show()

cs = CubicSpline(time_subset,theta_subset)

new_theta = np.copy(theta)
new_theta[indices_2018] = cs(time[indices_2018])

# plt.plot(time,theta)
# plt.plot(time,new_theta,'g-')
# plt.show()

output_file = folder+'/ECCOv5_THETA_IPCC_timeseries.nc'
ds = nc4.Dataset(output_file,'w')

ds.createDimension('time',np.shape(time)[0])

tvar = ds.createVariable('dec_yr','f4',('time',))
tvar[:] = time
mvar = ds.createVariable('THETA_mean', 'f4', ('time',))
mvar[:] = new_theta

ds.close()

