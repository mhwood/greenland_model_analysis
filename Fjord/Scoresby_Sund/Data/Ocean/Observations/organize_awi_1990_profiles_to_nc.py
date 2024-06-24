

import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt



def split_lines_to_individual_profiles(lines):
    a=1
    file_names = []
    dates = []
    longitudes = []
    latitudes =[]
    profiles = []

    current_id = lines[0].split()[0]
    year = 0
    month = 0
    day = 0
    hour = 0
    minute = 0
    longitude = 0
    latitude = 0
    profile = []
    for line in lines:
        line = line.split()
        if len(line)>2:
            id = line[0]
            if id!=current_id:
                if len(profile)>0:
                    file_name = 'CTD_'+str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)+'_'+\
                                '{:02d}'.format(hour)+'{:02d}'.format(minute)+'_'+id
                    file_names.append(file_name)
                    dates.append([year,month,day,hour,minute])
                    longitudes.append(longitude)
                    latitudes.append(latitude)
                    profiles.append(np.array(profile))
                    profile = []
                    current_id = id
            else:
                print(current_id)
                year = int(line[1].split('-')[0])
                month = int(line[1].split('-')[1])
                if month==6:
                    month=9
                day = int(line[1].split('-')[2][:2])
                hour = int(line[1].split('-')[2][3:5])
                minute = int(line[1].split('-')[2][6:8])
                longitude = float(line[3])
                latitude = float(line[2])
                depth = float(line[5])
                temp = float(line[6])
                salt = float(line[7])
                profile.append([depth,temp,salt])

    return(file_names, dates, longitudes, latitudes, profiles)

def output_file_as_nc(output_file, lon, lat,
                      year,month,day,hour,minute,
                      binned_depth,binned_potential_temperature,binned_practical_salinity,source,project_name):

    ds = nc4.Dataset(output_file,'w')
    ds.createDimension('records', len(binned_depth))
    dep = ds.createVariable('depth', 'f4', 'records')
    dep[:] = binned_depth
    salt = ds.createVariable('practical_salinity', 'f4', 'records')
    salt[:] = binned_practical_salinity
    ptemp = ds.createVariable('potential_temperature', 'f4', 'records')
    ptemp[:] = binned_potential_temperature

    ds.longitude = lon
    ds.latitude = lat
    ds.year = year
    ds.month = month
    ds.day = day
    ds.hour = hour
    ds.minute = minute
    ds.source = source
    ds.project_name = project_name

    ds.close()


ctd_dir = '/Users/michwood/Documents/Research/Projects/Scoresby Sund/Data/In Situ/AWI'
file_name = ctd_dir+'/Data/Raw/Stacked_CTDs_MW.txt'

f = open(file_name)
lines = f.read()
f.close()
lines = lines.split('\n')


file_names, dates, longitudes, latitudes, profiles = split_lines_to_individual_profiles(lines)

metadata_output = 'File_ID,Year,Month,Day,Hour,Minute,Longitude,Latitude'

for f in range(len(file_names)):
    print('    - Outputting '+file_names[f]+' shape: '+str(np.shape(profiles[f])))
    print('        - Location: '+str(longitudes[f])+', '+str(latitudes[f]))
    print('        - Date: '+str(dates[f]))
    output_file = ctd_dir+'/Data/1990/'+file_names[f]+'.nc'
    output_file_as_nc(output_file, longitudes[f], latitudes[f],
                      dates[f][0],dates[f][1],dates[f][2],dates[f][3],dates[f][4],
                      profiles[f][:,0],profiles[f][:,1],profiles[f][:,2],'AWI','AWI')
    metadata_output += '\n'+file_names[f]+','+str(dates[f][0])+','+str(dates[f][1])+','+str(dates[f][2])+','+\
                       str(dates[f][3])+','+str(dates[f][4])+','+str(longitudes[f])+','+str(latitudes[f])


f = open(ctd_dir+'/Metadata/AWI_1990_Scoresby_CTD_Locations.csv','w')
f.write(metadata_output)
f.close()