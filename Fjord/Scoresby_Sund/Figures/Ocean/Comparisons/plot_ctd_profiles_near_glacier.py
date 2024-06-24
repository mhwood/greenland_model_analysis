
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4

folder = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund'
model_output_folder = '/Volumes/helheim/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
                      'downscaled_greenland/L3/L3_Scoresby_Sund/results_iceplume_ks22/state_3D_mon_mean'

def find_closest_file_path(data_folder):

    djg_lon = -28.4851
    djg_lat = 71.9259

    dist=1e22
    closest_file_path = ''
    closest_lon = 0
    closest_lat = 0
    for file_name in os.listdir(data_folder):
        if file_name[-3:]=='.nc':
            ds=nc4.Dataset(data_folder+'/'+file_name)
            lon = ds.longitude
            lat = ds.latitude
            ds.close()

            dist_check = ((lon-djg_lon)**2 + (lat-djg_lat)**2)**0.5
            if dist_check<dist:
                dist=dist_check
                closest_file_path = data_folder+'/'+file_name
                closest_lon = lon
                closest_lat = lat

    # print('lon comparison: ',djg_lon,closest_lon)
    # print('lat comparison: ',djg_lat, closest_lat)

    return(closest_file_path)

def retrieve_CTD_profile(file_path):
    ds = nc4.Dataset(file_path)
    depth = ds.variables['depth'][:]
    T = ds.variables['potential_temperature'][:]
    S = ds.variables['practical_salinity'][:]
    ds.close()

    T_profile = np.column_stack([depth,T])
    S_profile = np.column_stack([depth,S])

    return(T_profile, S_profile)

def retieve_model_data(model_output_folder,year,month):

    for file_name in os.listdir(model_output_folder):
        if file_name[-3:]=='.nc' and file_name[0]!='.':
            yr = int(file_name.split('.')[1][:4])
            mo = int(file_name.split('.')[1][4:6])
            if yr==year and month==mo:
                file_path = model_output_folder+'/'+file_name

    djg_lon = -28.4851
    djg_lat = 71.9259

    ds=nc4.Dataset(file_path)
    depth = ds.variables['depths'][:]
    longitude = ds.variables['longitude'][:,:]
    latitude = ds.variables['latitude'][:,:]
    T = ds.variables['Theta'][:, :, :, :]
    S = ds.variables['Salt'][:, :, :, :]
    ds.close()

    dist_err = ((longitude-djg_lon)**2 + (latitude-djg_lat)**2)**0.5
    row, col = np.where(dist_err==np.min(dist_err))
    row = row[0]
    col = col[0]

    T = T[0, :, row,col]
    S = S[0, :, row, col]

    T_profile = np.column_stack([depth,T])
    S_profile = np.column_stack([depth,S])

    # print('lon comparison: ',djg_lon, longitude[row,col])
    # print('lat comparison: ',djg_lat, latitude[row,col])
    return(T_profile, S_profile)


def generate_CTD_comparison_plot(folder,
                                 hadley_T_profile, hadley_S_profile,
                                 omg_T_profiles, omg_S_profiles,
                                 first_model_T_profile, first_model_S_profile,
                                 model_T_profiles, model_S_profiles):

    output_file_path = folder+'/Figures/Ocean/Near-Glacier CTD Comparison/Near-Glacier CTD Comparison.png'

    fig = plt.figure(figsize=(10,10))

    plt.subplot(2,2,1)
    plt.plot(hadley_T_profile[:,1], hadley_T_profile[:,0], 'k-',label='1990')
    for p in range(len(omg_T_profiles)):
        omg_T_profile = omg_T_profiles[p]
        if p==0:
            plt.plot(omg_T_profile[:,1], omg_T_profile[:,0], 'r-',label='2020')
        else:
            plt.plot(omg_T_profile[:, 1], omg_T_profile[:, 0], 'r-')
    plt.gca().set_ylim([700, 0])
    plt.ylabel('Depth')
    plt.xlabel('Potential Temperature (($^{\circ}$C)')
    plt.legend(loc=3)
    plt.gca().set_xlim([-1.5,1.5])
    plt.text(1.45,10,'a)',ha='right',va='top',fontsize=14)
    plt.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    plt.title('CTD Observations')

    plt.subplot(2, 2, 2)
    plt.plot(first_model_T_profile[:, 1], first_model_T_profile[:, 0], 'k-', label='2000')
    for p in range(len(model_T_profiles)):
        model_T_profile = model_T_profiles[p]
        if p == 0:
            plt.plot(model_T_profile[:, 1], model_T_profile[:, 0], 'r-', label='2020')
        else:
            plt.plot(model_T_profile[:, 1], model_T_profile[:, 0], 'r-')
    plt.gca().set_ylim([700, 0])
    plt.ylabel('Depth')
    plt.xlabel('Potential Temperature (($^{\circ}$C)')
    plt.text(1.45, 10, 'b)', ha='right', va='top', fontsize=14)
    plt.legend(loc=3)
    plt.gca().set_xlim([-1.5, 1.5])
    plt.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    plt.title('L3 Model Output')

    plt.subplot(2, 2, 3)
    plt.plot(hadley_S_profile[:, 1], hadley_S_profile[:, 0], 'k-',label='1990')
    for p in range(len(omg_S_profiles)):
        omg_S_profile = omg_S_profiles[p]
        if p==0:
            plt.plot(omg_S_profile[:,1], omg_S_profile[:,0], 'b-',label='2020')
        else:
            plt.plot(omg_S_profile[:, 1], omg_S_profile[:, 0], 'b-')
    plt.gca().set_ylim([700, 0])
    plt.ylabel('Depth')
    plt.xlabel('Salinity (psu)')
    plt.text(34.9, 10, 'c)', ha='right', va='top', fontsize=14)
    plt.legend(loc=3)
    plt.gca().set_xlim([30,35])
    plt.grid(linestyle='--', linewidth=0.5, alpha=0.5)

    plt.subplot(2, 2, 4)
    plt.plot(first_model_S_profile[:, 1], first_model_S_profile[:, 0], 'k-', label='2000')
    for p in range(len(model_S_profiles)):
        model_S_profile = model_S_profiles[p]
        if p == 0:
            plt.plot(model_S_profile[:, 1], model_S_profile[:, 0], 'b-', label='2020')
        else:
            plt.plot(model_S_profile[:, 1], model_S_profile[:, 0], 'b-')
    plt.gca().set_ylim([700, 0])
    plt.ylabel('Depth')
    plt.xlabel('Salinity (psu)')
    plt.legend(loc=3)
    plt.gca().set_xlim([30, 35])
    plt.text(34.9, 10, 'd)', ha='right', va='top', fontsize=14)
    plt.grid(linestyle='--',linewidth=0.5, alpha=0.5)

    plt.suptitle('Temperature and Salinity Near Daugaard-Jensen Glacier\n(28.49$^{\circ}$W, 71.93$^{\circ}$N)')

    plt.savefig(output_file_path, bbox_inches='tight')
    plt.close(fig)



hadley_path = find_closest_file_path(folder+'/Data/In Situ/AWI/Data/1990')
hadley_T_profile, hadley_S_profile = retrieve_CTD_profile(hadley_path)

omg_T_profiles = []
omg_S_profiles = []
for year in ['2020']:#'2017','2018','2019'
    omg_path = find_closest_file_path(folder + '/Data/In Situ/OMG/Data/'+year)
    omg_T_profile, omg_S_profile = retrieve_CTD_profile(omg_path)
    omg_T_profiles.append(omg_T_profile)
    omg_S_profiles.append(omg_S_profile)

first_model_T_profile, first_model_S_profile = retieve_model_data(model_output_folder,2000,8)

model_T_profiles = []
model_S_profiles = []
for year in range(2020,2021):
    model_T_profile, model_S_profile = retieve_model_data(model_output_folder,year,8)
    model_T_profiles.append(model_T_profile)
    model_S_profiles.append(model_S_profile)

generate_CTD_comparison_plot(folder,
                             hadley_T_profile, hadley_S_profile,
                             omg_T_profiles, omg_S_profiles,
                             first_model_T_profile, first_model_S_profile,
                             model_T_profiles, model_S_profiles)


