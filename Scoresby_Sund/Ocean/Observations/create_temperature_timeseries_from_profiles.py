
import os
import netCDF4 as nc4
import numpy as np
import shutil
import datetime
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import gsw

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def get_hadley_profile_list(project_dir, source, location, argo_only):

    metadata_file = os.path.join(project_dir,'Data','In Situ',source,'Metadata',source+'_CTD_Locations_'+location+'.csv')
    if os.path.exists(metadata_file):

        f = open(metadata_file)
        lines = f.read()
        f.close()
        lines = lines.split('\n')
        lines.pop(0)

        profile_list = []
        for line in lines:
            line = line.split(',')
            if len(line)>2:
                if argo_only:
                    if 'ARGO' in line[0]:
                        profile_list.append([line[2], 'CTD_'+line[0]])
                else:
                    profile_list.append([line[2], 'CTD_' + line[0]])
    else:
        profile_list = []

    return(profile_list)

def get_awi_profile_list(location):

    if location=='Sound_Entrance':
        profile_list =[['1990','CTD_19900925_0000_PS1944-1']]
    if location=='Mid_Sound':
        profile_list =[['1990','CTD_19900922_0000_PS1943-1']]
    if location=='Mid_Fjord':
        profile_list =[['1990','CTD_19900921_0000_PS1940-1']]
    if location=='Near_Glacier':
        profile_list =[['1990','CTD_19900921_0000_PS1938-1']]

    return(profile_list)

def get_omg_profile_list(project_dir,location):

    omg_ctd_dir = '/Volumes/zachariae/Research/Data Repository/Ocean Properties/OMG_CTDs/Processed'

    if location=='sound_entrance':
        profile_list =[['2016','CTD_20161010_145001'],
                       ['2017','CTD_20171015_123129'],
                       ['2018','CTD_20180825_140502'],
                       ['2019','CTD_20190819_134458'],
                       ['2020','CTD_20200904_155014'],
                       ['2021','CTD_20210825_125131']]
    if location=='mid_sound':
        profile_list = [['2016', 'CTD_20161010_140634'],
                        ['2017', 'CTD_20171016_130337'],
                        ['2018', 'CTD_20180827_152500'],
                        ['2019', 'CTD_20190819_143308'],
                        ['2020', 'CTD_20200904_175044'],
                        ['2021', 'CTD_20210820_134800']]
    if location=='mid_fjord':
        profile_list = [['2018', 'CTD_20180827_145527'],
                        ['2019', 'CTD_20190819_145857'],
                        ['2020', 'CTD_20200904_181248'],
                        ['2021', 'CTD_20210820_132835']]
    if location=='near_glacier':
        profile_list = [['2017', 'CTD_20171016_121232'],
                        ['2018', 'CTD_20180827_143507'],
                        ['2019', 'CTD_20190819_151838'],
                        ['2020', 'CTD_20200904_183011'],
                        ['2021', 'CTD_20210820_131357']]

    for file_set in profile_list:
        if file_set[1]+'.nc' not in os.listdir(os.path.join(project_dir,'Data','In Situ','OMG','Data',file_set[0])):
            source_file = os.path.join(omg_ctd_dir,file_set[0],file_set[1]+'.nc')
            dest_file = os.path.join(project_dir,'Data','In Situ','OMG','Data', file_set[0], file_set[1] + '.nc')
            shutil.copyfile(source_file,dest_file)

    return(profile_list)

def read_temperature_timeseries_from_profiles_by_depth(project_dir, source, profile_list, min_depth, max_depth):

    profile_dir = os.path.join(project_dir, 'Data', 'In Situ', source, 'Data')

    dec_yrs = []
    temperature_means = []
    temperature_stds = []

    interp_dpeths = np.arange(min_depth,max_depth)

    for profile_set in profile_list:
        year = profile_set[0]
        if int(year)<=2019:
            file_name = profile_set[1]+'.nc'
            ds = nc4.Dataset(os.path.join(profile_dir,year,file_name))
            depth = ds.variables['depth'][:]
            temp = ds.variables['potential_temperature'][:]
            year = ds.year
            month = ds.month
            day = ds.day
            hour = ds.hour
            minute = int(str(ds.minute)[:2])
            print(year,month,day,hour,minute)
            if day==32:
                day = 31
            dec_yr = YMD_to_DecYr(year, month, day, hour, minute)
            ds.close()
            if np.min(depth)<=min_depth and np.max(depth)>=max_depth:
                dec_yrs.append(dec_yr)
                set_int = interp1d(depth,temp)
                interp_temp = set_int(interp_dpeths)
                temperature_means.append(np.mean(interp_temp))
                temperature_stds.append(np.std(interp_temp))

    timeseries = np.column_stack([np.array(dec_yrs),
                                  np.array(temperature_means),
                                  np.array(temperature_stds)])

    return(timeseries)

def read_temperature_timeseries_from_profiles_by_density(project_dir, source, profile_list, min_density, max_density):

    profile_dir = os.path.join(project_dir, 'Data', 'In Situ', source, 'Data')

    dec_yrs = []
    temperature_means = []
    temperature_stds = []

    interp_densities = np.linspace(min_density,max_density,100)

    for profile_set in profile_list:
        year = profile_set[0]
        if int(year)<=2021:
            file_name = profile_set[1]+'.nc'
            ds = nc4.Dataset(os.path.join(profile_dir,year,file_name))
            depth = ds.variables['depth'][:]
            temp = ds.variables['potential_temperature'][:]
            salt = ds.variables['practical_salinity'][:]
            CT = gsw.CT_from_pt(salt, temp)
            pressure = gsw.p_from_z(-1*depth, ds.latitude)
            density = gsw.rho(salt, CT, pressure)
            year = ds.year
            month = ds.month
            day = ds.day
            hour = ds.hour
            minute = int(str(ds.minute)[:2])
            if day==32:
                day = 31
            dec_yr = YMD_to_DecYr(year, month, day, hour, minute)
            ds.close()

            dec_yrs.append(dec_yr)
            set_int = interp1d(density,temp)
            interp_temp = set_int(interp_densities)
            temperature_means.append(np.mean(interp_temp))
            temperature_stds.append(np.std(interp_temp))

    timeseries = np.column_stack([np.array(dec_yrs),
                                  np.array(temperature_means),
                                  np.array(temperature_stds)])

    return(timeseries)


def output_timeseries_to_nc(output_file,timeseries,min_depth,max_depth):

    ds = nc4.Dataset(output_file,'w')
    ds.createDimension('time', np.shape(timeseries)[0])

    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = timeseries[:,0]
    tvar = ds.createVariable('temperature_mean', 'f4', ('time',))
    tvar[:] = timeseries[:, 1]
    tvar = ds.createVariable('temperature_std', 'f4', ('time',))
    tvar[:] = timeseries[:, 2]

    ds.min_depth = min_depth
    ds.max_depth = max_depth

    ds.close()

project_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund'

# source = 'Hadley'
# source = 'AWI'
source = 'OMG'

locations = ['near_glacier','sound_entrance','mid_sound']
densities = [[1029,1030.1],[1029,1030.1],[1029,1029.7]]

min_depth = 200
max_depth = 400

for ll in range(len(locations)):
    location = locations[ll]

    if source=='Hadley':
        argo_only = True
        profile_list = get_hadley_profile_list(project_dir, source, location, argo_only)
    if source=='AWI':
        profile_list = get_awi_profile_list(location)
    if source=='OMG':
        profile_list = get_omg_profile_list(project_dir, location)

    min_density = densities[ll][0]
    max_density = densities[ll][1]

    #timeseries = read_temperature_timeseries_from_profiles_by_depth(project_dir, source, profile_list, min_depth, max_depth)
    timeseries = read_temperature_timeseries_from_profiles_by_density(project_dir, source, profile_list, min_density, max_density)

    output_file = os.path.join(project_dir,'Data','In Situ',source,'Data',source+'_'+location+'_Temperature_Timeseries.nc')
    output_timeseries_to_nc(output_file,timeseries,min_depth,max_depth)

    plt.plot(timeseries[:,0],timeseries[:,1],'.',color='purple')
    for i in range(np.shape(timeseries)[0]):
        plt.plot([timeseries[i,0],timeseries[i,0]],
                 [timeseries[i,1]-timeseries[i,2],timeseries[i,1]+timeseries[i,2]],
                 '-',color='purple')
    plt.show()









