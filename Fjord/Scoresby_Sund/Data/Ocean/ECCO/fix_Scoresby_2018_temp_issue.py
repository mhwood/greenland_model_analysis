

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import CubicSpline


folder = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/Data/Modeling/' \
         'ECCO'

ds = nc4.Dataset(folder+'/ECCOv5_Scoresby_Sund_Timeseries_2018_issue.nc')
time = ds.variables['time'][:]
theta = ds.groups['sound_entrance'].variables['THETA'][:]
ds.close()

indices_2018 = np.floor(time).astype(int)==2018
time_subset = time[~indices_2018]
theta_subset = theta[~indices_2018]

plt.plot(time_subset, theta_subset, 'k.')
plt.show()

cs = CubicSpline(time_subset,theta_subset)

new_theta = np.copy(theta)
new_theta[indices_2018] = cs(time[indices_2018])

plt.plot(time,theta)
plt.plot(time,new_theta,'g-')
plt.show()

output_file = folder+'/ECCOv5_Scoresby_Sund_Timeseries.nc'
ds = nc4.Dataset(output_file,'w')

ds.createDimension('time',np.shape(time)[0])

tvar = ds.createVariable('time','f4',('time',))
tvar[:] = time

grp = ds.createGroup('sound_entrance')
mvar = grp.createVariable('THETA', 'f4', ('time',))
mvar[:] = new_theta

ds.close()

