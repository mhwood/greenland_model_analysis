
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

# def read_sample_area_from_shapefile(project_dir, polygon_name):
#     shp_path = os.path.join(project_dir,'Map','Shapefiles','Ocean_Sample_Areas','Ocean_Sample_Areas')
#
#     r = shapefile.Reader(shp_path)
#     shapes = r.shapes()
#     records = r.records()
#
#     for r in range(len(records)):
#         if records[r][0]==polygon_name:
#             sample_area = np.array(shapes[r].points)
#
#     sample_area = reproject_polygon(sample_area,3413,4326,x_column=1, y_column=0)
#
#     return(sample_area)

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

def generate_mean_model_timeseries(model_dir, min_depth, max_depth, sample_area):

    state_dir = os.path.join(model_dir,'results','monthly_nc','state_3D_mon_mean')
    annual_files = []
    for file_name in os.listdir(state_dir):
        if file_name[0]!='.' and file_name[-3:]=='.nc':
            annual_files.append(file_name)
    annual_files = sorted(annual_files)

    p = mplPath.Path(sample_area)
    integration_depths = np.arange(min_depth,max_depth)

    timeseries_started = False

    for ff in range(len(annual_files)):
        print('    - Adding data from '+annual_files[ff])
        ds = nc4.Dataset(os.path.join(state_dir, annual_files[ff]))
        lon = ds.variables['longitude'][:,:]
        lat = ds.variables['latitude'][:,:]
        Theta = ds.variables['Theta'][:,:,:,:]
        time = ds.variables['time'][:]
        depth = ds.variables['depth'][:]
        ds.close()

        depth = np.reshape(depth,(len(depth),1))

        file_timeseries = np.zeros((len(time),1))
        file_decyrs = model_time_to_dec_yr(time)

        inside = p.contains_points(np.column_stack([lon.ravel(),lat.ravel()]))
        inside = np.reshape(inside,np.shape(lon))

        for timestep in range(len(time)):
            time_slice = Theta[:,timestep,:,:]
            mean_profile = np.zeros((len(depth), 1))
            std_profile = np.zeros((len(depth), 1))

            for d in range(len(depth)):
                # if depth_indices[d]:
                depth_slice = time_slice[d,:,:]
                inside_points = depth_slice[inside]
                if np.any(inside_points!=0):
                    inside_points = inside_points[inside_points!=0]
                    mean_profile[d] = np.mean(inside_points)
                    std_profile[d] = np.std(inside_points)

            # plt.subplot(1,2,1)
            # plt.plot(mean_profile,depth)
            # plt.subplot(1, 2, 2)
            # plt.plot(std_profile, depth)
            # plt.show()

            set_int = interp1d(depth.ravel(),mean_profile.ravel())
            file_timeseries[timestep] = np.mean(set_int(integration_depths))

        # plt.plot(file_decyrs, file_timeseries)
        # plt.show()

        if not timeseries_started:
            full_timeseries = np.column_stack([file_decyrs,file_timeseries])
            timeseries_started = True
        else:
            full_timeseries = np.vstack([full_timeseries,
                                         np.column_stack([file_decyrs,file_timeseries])])

    # plt.plot(full_timeseries[:,0], full_timeseries[:,1])
    # plt.show()

    return(full_timeseries)

def store_timeseries_as_nc(output_dir,output_file,timeseries):

    ds = nc4.Dataset(os.path.join(output_dir,output_file),'w')
    ds.createDimension('time',np.shape(timeseries)[0])

    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = timeseries[:,0]

    svar = ds.createVariable('Theta','f4',('time',))
    svar[:] = timeseries[:,1]

    ds.close()




project_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund'

model_dir = '/Volumes/helheim/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
            'downscaled_greenland/L1/grid/L1_CE_Greenland'

shapefile_name = 'IPCC Sample Area Off Shelf'
polygon_name = 'IPCC'

# get sample area from shapefile
sample_area = read_sample_area_from_shapefile(project_dir, shapefile_name)

# loop through the model results annually to make a timeseries of the mean temperature 200-500 m
min_depth = 200
max_depth = 500
timeseries = generate_mean_model_timeseries(model_dir, min_depth, max_depth, sample_area)

# store timeseries as nc
output_dir = os.path.join(project_dir,'Data','Modeling','Downscaled','L1_CE_Greenland')
output_file = 'L1_CE_Greenland_'+polygon_name+'_Timeseries.nc'
store_timeseries_as_nc(output_dir,output_file,timeseries)

