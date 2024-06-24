
import requests
import shutil
from datetime import datetime

def date_to_timestep(date):
    total_days = (((date-datetime(2002,7,4)).total_seconds())/(24*60*60))+915
    print(total_days)
    timestep = (total_days-915)/8
    return(timestep)

date = datetime(2022,7,20)
# date = datetime(2002,7,12)
timestep = date_to_timestep(date)
print(timestep)

data_set = 'MYD09A1.006'
tile_id = 'h17v01'

def download_file(url,output_file):
    print(url)
    # local_filename = url.split('/')[-1]
    with requests.get(url, stream=True) as r:
        with open(output_file, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


url = 'https://opendap.cr.usgs.gov/opendap/hyrax/MYD09A1.006/h17v01.ncml.nc4?' \
      'Latitude%5B0:1:2399%5D%5B0:1:2399%5D,' \
      'Longitude%5B0:1:2399%5D%5B0:1:2399%5D,' \
      'sur_refl_b01%5B0:1:1%5D%5B0:1:2399%5D%5B0:1:2399%5D,' \
      'sur_refl_b03%5B0:1:1%5D%5B0:1:2399%5D%5B0:1:2399%5D,' \
      'sur_refl_b04%5B0:1:1%5D%5B0:1:2399%5D%5B0:1:2399%5D,' \
      'time%5B0:1:1%5D'

output_file = '/Users/michwood/Documents/Research/Projects/Scoresby Sund/Data/Remote Sensing/Raw/MODIS/test.nc4'

download_file(url,output_file)

file_url  = 'https://opendap.cr.usgs.gov/opendap/hyrax/'
file_url += data_set+'/'
file_url += tile_id
file_url += '.ncml.nc4?'
file_url += 'Latitude[0:1:2399][0:1:2399],'
file_url += 'Longitude[0:1:2399][0:1:2399],'
file_url += 'sur_refl_b01[0:1:3][0:1:2399][0:1:2399],'
file_url += 'sur_refl_b03[0:1:3][0:1:2399][0:1:2399],'
file_url += 'sur_refl_b04[0:1:3][0:1:2399][0:1:2399],'
file_url += 'time[0:1:3]'

# links
#https://modis.gsfc.nasa.gov/data/dataprod/mod09.php
#https://opendap.cr.usgs.gov/opendap/hyrax/MYD09A1.006/contents.html
#https://opendap.cr.usgs.gov/opendap/hyrax/MYD09A1.006/h17v01.ncml.html
#https://modis-land.gsfc.nasa.gov/MODLAND_grid.html # this is the grid ref

