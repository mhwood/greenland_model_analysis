
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
# import cm.ocean as cm
import datetime


def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

config_folder = '/Volumes/helheim/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
                'configurations/downscaled_greenland/L3/L3_Scoresby_Sund'



file_path = os.path.join(config_folder,'results_iceplume_ks22','Scoresby_Sund_transects','Scoresby_Sund_transects.nc')

ds = nc4.Dataset(file_path)
years = ds.variables['years'][:]
months = ds.variables['months'][:]
depth = ds.variables['depth'][:]
dist = ds.variables['transect_distance'][:]
temp = ds.variables['transect_temp'][:,:]
ds.close()

time = np.arange(len(years))
time = []
for t in range(len(years)):
    time.append(YMD_to_DecYr(int(years[t]),int(months[t]),15))
time = np.array(time)

plot_anomaly = False

if plot_anomaly:
    output_file = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/Figures/Ocean/Transects/' \
                  'Scoresby_Sund_THETA_anomaly_profile.png'
else:
    output_file = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/Figures/Ocean/Transects/' \
                  'Scoresby_Sund_THETA_profile.png'

fig = plt.figure(figsize=(14,6))

if plot_anomaly:
    mean_profile = np.mean(temp[:,:,-10], axis=0)
    temp_anomaly = np.copy(temp[:,:,-10])
    for t in range(np.shape(temp)[0]):
        temp_anomaly[t,:] -= mean_profile
    C = plt.pcolormesh(time,depth,temp_anomaly.T, cmap = 'seismic',vmin=-1.5, vmax=1.5)
    plt.colorbar(C)
    plt.title('Fjrod Temperature Anomaly near Daugaard Jensen Glacier')
else:
    C = plt.pcolormesh(time, depth, temp[:, :, -10].T, cmap='turbo')
    plt.colorbar(C)
    plt.title('Fjrod Temperature near Daugaard Jensen Glacier')
plt.gca().set_ylim([600,0])
plt.ylabel('Depth (m)')
plt.gca().set_xticks(np.arange(2001,2021,2))

plt.savefig(output_file, bbox_inches='tight')
plt.close(fig)