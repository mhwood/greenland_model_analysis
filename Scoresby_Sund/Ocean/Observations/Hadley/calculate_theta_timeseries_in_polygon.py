
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import shapefile
import netCDF4 as nc4
import datetime
from scipy.interpolate import interp1d
from pyproj import Proj, Transformer

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1):
    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))
    # inProj = Proj(init='epsg:'+str(inputCRS))
    # outProj = Proj(init='epsg:'+str(outputCRS))
    x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon

def read_sample_area_from_shapefile(project_dir, shapefile_name):
    shp_path = os.path.join(project_dir,'Map','Shapefiles','Ocean_Sample_Areas',shapefile_name)

    r = shapefile.Reader(shp_path)
    shapes = r.shapes()
    records = r.records()

    sample_area = np.array(shapes[0].points)

    # sample_area = reproject_polygon(sample_area,3413,4326,x_column=1, y_column=0)

    return(sample_area)

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def model_time_to_dec_yr(time):
    file_decyrs = np.zeros((len(time),1))
    for t in range(len(time)):
        date = datetime.datetime(1992,1,1) + datetime.timedelta(seconds = int(time[t]))
        dec_yr = YMD_to_DecYr(date.year,date.month,date.day)
        file_decyrs[t] = dec_yr
    return(file_decyrs)

def generate_mean_hadley_timeseries(hadley_analysis_dir, min_depth, max_depth, sample_area):

    p = mplPath.Path(sample_area)
    integration_depths = np.arange(min_depth,max_depth)

    years = np.arange(1991,2023).tolist()
    time = np.zeros((12*len(years)))
    timeseries = np.zeros((12 * len(years)))

    counter = 0
    for year in years:
        if year==2022:
            max_month = 9
        else:
            max_month = 13
        for month in range(1,max_month):
            time[counter] = YMD_to_DecYr(year,month,15)

            month_file = 'EN.4.2.2.f.analysis.g10.'+str(year)+'{:02d}'.format(month)+'.nc'
            # month_file = 'EN.4.2.1.f.analysis.l09.' + str(year) + '{:02d}'.format(month) + '.nc'
            print('    - Adding data from '+month_file)
            ds = nc4.Dataset(os.path.join(hadley_analysis_dir,'g10', 'EN.4.2.2.analyses.g10.'+str(year), month_file))
            # ds = nc4.Dataset(os.path.join(hadley_analysis_dir,'l09', 'EN.4.2.1.analyses.l09.' + str(year), month_file))
            lon = ds.variables['lon'][:]
            lat = ds.variables['lat'][:]
            Theta = ds.variables['temperature'][:,:,:,:]
            depth = ds.variables['depth'][:]
            ds.close()

            Lon, Lat = np.meshgrid(lon, lat)
            Lon = np.hstack([Lon[:, -180:], Lon[:, :180]])
            Lon[Lon > 180] -= 360
            Lat = np.hstack([Lat[:, -180:], Lat[:, :180]])
            Theta = np.concatenate([Theta[0, :, :, -180:], Theta[0, :, :, :180]], axis=2)

            depth = np.reshape(depth,(len(depth),1))

            inside = p.contains_points(np.column_stack([Lon.ravel(),Lat.ravel()]))
            inside = np.reshape(inside,np.shape(Lon))

            mean_profile = np.zeros((len(depth), 1))
            std_profile = np.zeros((len(depth), 1))

            for d in range(len(depth)):
                # if depth_indices[d]:
                depth_slice = Theta[d,:,:]
                inside_points = depth_slice[inside]
                if np.any(inside_points>-2):
                    inside_points = inside_points[inside_points>-2]
                    mean_profile[d] = np.mean(inside_points)
                    std_profile[d] = np.std(inside_points)

            mean_profile -= 273.15

            set_int = interp1d(depth.ravel(),mean_profile.ravel())
            mean_temp = np.mean(set_int(integration_depths))
            timeseries[counter] = mean_temp

            # plt.subplot(1, 2, 1)
            # plt.plot(mean_profile[depth < 2000], depth[depth < 2000])
            # plt.gca().invert_yaxis()
            # plt.title('Temp mean ('+str(mean_temp)+')')
            # plt.subplot(1, 2, 2)
            # plt.plot(std_profile[depth < 2000], depth[depth < 2000])
            # plt.gca().invert_yaxis()
            # plt.title('Temp std')
            # plt.show()

            counter +=1

    return(time, timeseries)

def store_timeseries_as_nc(output_dir,output_file,time, timeseries):

    ds = nc4.Dataset(os.path.join(output_dir,output_file),'w')
    ds.createDimension('time',np.shape(timeseries)[0])

    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = time

    svar = ds.createVariable('Theta','f4',('time',))
    svar[:] = timeseries

    ds.close()




project_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund'

hadley_analysis_dir = '/Volumes/helheim/Data_Repository/Hadley Analysis'
polygon_name = 'IPCC'

# get sample area from shapefile
shapefile_name = 'IPCC Sample Area Off Shelf'
sample_area = read_sample_area_from_shapefile(project_dir, shapefile_name)

# loop through the model results annually to make a timeseries of the mean temperature 200-500 m
min_depth = 200
max_depth = 500
time, timeseries = generate_mean_hadley_timeseries(hadley_analysis_dir, min_depth, max_depth, sample_area)

# store timeseries as nc
output_dir = os.path.join(project_dir,'Data','In Situ','Hadley','Timeseries')
output_file = 'Hadley_'+polygon_name+'_Timeseries.nc'
store_timeseries_as_nc(output_dir,output_file,time, timeseries)

