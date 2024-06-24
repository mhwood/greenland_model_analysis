

import os
import numpy as np
import netCDF4 as nc4
import shapefile

folder = '/Users/mike/Documents/Research/Projects/Disko Bay'

lons = []
lats = []
dates = []

csv_output = 'file_name,longitude,latitude,date'

for file_name in os.listdir(folder+'/Data/Nitrate/OSD'):
    if file_name[-3:]=='.nc':

        ds=nc4.Dataset(os.path.join(folder,'Data/Nitrate/OSD',file_name))
        if 'date' in list(ds.variables.keys()):
            lon = ds.geospatial_lon_min
            lat = ds.geospatial_lat_min
            date = ds.variables['date'][0]
            lons.append(lon)
            lats.append(lat)
            dates.append(date)
            csv_output+='\n'+file_name+','+str(lon)+','+str(lat)+','+str(date)
        ds.close()


output_file = os.path.join(folder,'Map','Shapefiles','Nitrate_OSD_locations')
sf = shapefile.Writer(output_file)

sf.field('date','N')

for r in range(len(lons)):
    sf.point(lons[r], lats[r])
    sf.record(dates[r])

sf.close()

f = open(output_file+'.csv','w')
f.write(csv_output)
f.close()

prj_file = output_file+'.prj'
f = open(prj_file,'w')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],'
        'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()



