
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import shapefile
from scipy.interpolate import griddata
from pyproj import Transformer

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1,run_test = True):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(polygon_array[:, x_column], polygon_array[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
        run_test = False
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon

def read_model_geometry_onto_line(config_dir, model_name, line):
    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids',model_name+'_grid.nc'))
    drF = ds.variables['drF'][:]
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    hFacC = ds.variables['HFacC'][:, :, :]
    ds.close()

    points = np.column_stack([XC.ravel(), YC.ravel()])

    line_depth = griddata(points, Depth.ravel(), (line[:, 0], line[:, 1]), fill_value=0, method='nearest')

    line_hFacC = np.zeros((len(drF),np.shape(line)[0]))
    for d in range(len(drF)):
        interp_values = griddata(points, hFacC[d, :, :].ravel(), (line[:, 0], line[:, 1]), fill_value=0, method='nearest')
        line_hFacC[d, :] = interp_values

    return(drF, line_depth, line_hFacC)

def series_to_N_points(series,N):
    #find the total length of the series
    totalDistance=0
    for s in range(len(series[:,0])-1):
        totalDistance+=((series[s,0]-series[s+1,0])**2+(series[s,1]-series[s+1,1])**2)**0.5
    intervalDistance=totalDistance/(N-1)

    #make the list of points
    newSeries=series[0,:]
    currentS = 0
    currentPoint1=series[currentS,:]
    currentPoint2=series[currentS+1,:]
    for p in range(N-2):
        distanceAccrued = 0
        while distanceAccrued<intervalDistance:
            currentLineDistance=((currentPoint1[0]-currentPoint2[0])**2+(currentPoint1[1]-currentPoint2[1])**2)**0.5
            if currentLineDistance<intervalDistance-distanceAccrued:
                distanceAccrued+=currentLineDistance
                currentS+=1
                currentPoint1 = series[currentS, :]
                currentPoint2 = series[currentS + 1, :]
            else:
                distance=intervalDistance-distanceAccrued
                newX=currentPoint1[0]+(distance/currentLineDistance)*(currentPoint2[0]-currentPoint1[0])
                newY = currentPoint1[1] + (distance / currentLineDistance) * (currentPoint2[1] - currentPoint1[1])
                distanceAccrued=intervalDistance+1
                newSeries=np.vstack([newSeries,np.array([newX,newY])])
                currentPoint1=np.array([newX,newY])
    newSeries = np.vstack([newSeries, series[-1,:]])
    return(newSeries)

def get_depth():
    depth = np.array([5.0000000e+00, 1.5000000e+01, 2.5000000e+01, 3.5000000e+01, 4.5000000e+01,
                      5.5000000e+01, 6.5000000e+01, 7.5005005e+01, 8.5025002e+01, 9.5095001e+01,
                      1.0531000e+02, 1.1587000e+02, 1.2715000e+02, 1.3973999e+02, 1.5447000e+02,
                      1.7239999e+02, 1.9473500e+02, 2.2271001e+02, 2.5747000e+02, 2.9992999e+02,
                      3.5067999e+02, 4.0992999e+02, 4.7747000e+02, 5.5271002e+02, 6.3473505e+02,
                      7.2240002e+02, 8.1447009e+02, 9.0974011e+02, 1.0071550e+03, 1.1059050e+03,
                      1.2055350e+03, 1.3062051e+03, 1.4091499e+03, 1.5170950e+03, 1.6341748e+03,
                      1.7651348e+03, 1.9141499e+03, 2.0840349e+03, 2.2762251e+03, 2.4912500e+03,
                      2.7292500e+03, 2.9902500e+03, 3.2742500e+03, 3.5812500e+03, 3.9112500e+03,
                      4.2642500e+03, 4.6402500e+03, 5.0392500e+03, 5.4612500e+03, 5.9062500e+03])
    return(depth)

def great_circle_distance(lon_ref, lat_ref, Lon, Lat):
    earth_radius = 6371000
    lon_ref_radians = np.radians(lon_ref)
    lat_ref_radians = np.radians(lat_ref)
    lons_radians = np.radians(Lon)
    lats_radians = np.radians(Lat)
    lat_diff = lats_radians - lat_ref_radians
    lon_diff = lons_radians - lon_ref_radians
    d = np.sin(lat_diff * 0.5) ** 2 + np.cos(lat_ref_radians) * np.cos(lats_radians) * np.sin(lon_diff * 0.5) ** 2
    h = 2 * earth_radius * np.arcsin(np.sqrt(d))
    return(h)

def read_transect_from_shapefile(project_folder):
    sf = shapefile.Reader(os.path.join(project_folder,'Map','Shapefiles','Greenland Shelf Break'))
    line = np.array(sf.shapes()[0].points)
    line = series_to_N_points(line, N=1000)

    dist = np.zeros((np.shape(line)[0], ))
    for d in range(len(dist)-1):
        dist[d+1] = dist[d] + great_circle_distance(line[d,0], line[d,1], line[d+1,0], line[d+1,1])
    return(line, dist)

def interpolate_model_output_onto_line(results_dir, var_name, year, month, line):

    date_str = str(year)+'{:02d}'.format(month)
    output_grid = np.zeros((50, np.shape(line)[0]))

    file_name_to_read = ''
    for file_name in os.listdir(results_dir):
        if file_name[0]!='.' and file_name[-3:]=='.nc':
            if file_name.split('.')[1]==date_str:
                file_name_to_read = file_name

    if file_name_to_read!='':

        ds = nc4.Dataset(os.path.join(results_dir, file_name_to_read))
        Lon = ds.variables['longitude'][:, :]
        Lat = ds.variables['latitude'][:, :]
        depth = ds.variables['depths'][:]
        var_grid = ds.variables[var_name][:, :, :, :]
        ds.close()

        points = np.column_stack([Lon.ravel(), Lat.ravel()])

        for d in range(len(depth)):
            # print(d)
            interp_values = griddata(points, var_grid[0,d,:,:].ravel(), (line[:,0], line[:,1]), fill_value=0, method='nearest')
            output_grid[d,:] = interp_values

    return(output_grid)

def output_annual_nc_file(project_folder, var_name, year, line, distance, depth, year_output,
                          drF, line_depth, line_hFacC):

    line_3413 = reproject_polygon(line,4326,3413)

    output_folder = os.path.join(project_folder,'Data','Shelf Break Timeseries')

    if var_name not in os.listdir(output_folder):
        os.mkdir(os.path.join(output_folder,var_name))

    months = np.arange(1, 13)

    ds = nc4.Dataset(os.path.join(output_folder,var_name,'Shelf_Break_'+var_name+'_'+str(year)+'.nc'), 'w')

    ds.createDimension('months',len(months))
    ds.createDimension('distance',len(distance))
    ds.createDimension('depth',len(depth))

    var = ds.createVariable('months','i4',('months',))
    var[:] = months

    var = ds.createVariable('depth', 'f4', ('depth',))
    var[:] = depth

    var = ds.createVariable('bathymetry', 'f4', ('distance',))
    var[:] = line_depth

    var = ds.createVariable('drF', 'f4', ('depth',))
    var[:] = drF

    var = ds.createVariable('distance', 'f4', ('distance',))
    var[:] = distance

    var = ds.createVariable('longitude', 'f4', ('distance',))
    var[:] = line[:,0]

    var = ds.createVariable('latitude', 'f4', ('distance',))
    var[:] = line[:,1]

    var = ds.createVariable('x', 'f4', ('distance',))
    var[:] = line_3413[:, 0]

    var = ds.createVariable('y', 'f4', ('distance',))
    var[:] = line_3413[:, 1]

    var = ds.createVariable('hFacC', 'f4', ('depth', 'distance'))
    var[:, :] = line_hFacC

    var = ds.createVariable(var_name, 'f4', ('months', 'depth', 'distance'))
    var[:, :, :] = year_output

    ds.close()


project_folder = '/Users/mhwood/Documents/Research/Projects/Greenland Shelf Transport'

W_Greenland_config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/' \
                          'configurations/downscale_darwin'

N_Greenland_config_dir = '/Volumes/petermann/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
                          'configurations/downscale_greenland'

break_line, distance = read_transect_from_shapefile(project_folder)

print('    - Computing geometry along transect')
drF, line_depth_N, line_hFacC_N = read_model_geometry_onto_line(N_Greenland_config_dir, 'L1_N_Greenland', break_line)
_, line_depth_W, line_hFacC_W = read_model_geometry_onto_line(W_Greenland_config_dir, 'L1_W_Greenland', break_line)

line_hFacC = np.zeros_like(line_hFacC_W)
line_depth = np.zeros_like(line_depth_W)

W_indices = distance > 3000e3
for d in range(np.shape(line_hFacC_W)[0]):
    line_hFacC[d, W_indices] = line_hFacC_W[d, W_indices]
line_depth[W_indices] = line_depth_W[W_indices]

N_indices = distance < 3000e3
for d in range(np.shape(line_hFacC_N)[0]):
    line_hFacC[d, N_indices] = line_hFacC_N[d, N_indices]
line_depth[N_indices] = line_depth_N[N_indices]

# plt.subplot(2,1,1)
# plt.plot(line_depth)
# plt.subplot(2,1,2)
# C = plt.pcolormesh(line_hFacC)
# plt.colorbar(C)
# plt.show()

depth = get_depth()

for year in range(1992, 2022):
    print('    - Interpolating the grid in year ' + str(year))

    for var_name in ['Uvel', 'Vvel']:#,'Theta']:
        if var_name in ['Theta', 'Salt']:
            subset = 'state_3D_mon_mean'
        if var_name in ['Uvel', 'Vvel']:
            subset = 'vel_3D_mon_mean'

        results_dir_W_Greenland = W_Greenland_config_dir+'/L1/L1_W_Greenland/results/' + subset

        results_dir_N_Greenland = N_Greenland_config_dir+'/L1/L1_N_Greenland/results/' + subset

        year_output = np.zeros((12, 50, np.shape(break_line)[0]))

        for month in range(1,13):
            print('        - Working on '+var_name+' in month '+str(month))

            print('            - Interpolating the W Greenland results')
            theta_grid_W_Greenland = interpolate_model_output_onto_line(results_dir_W_Greenland, var_name, year, month, break_line)

            print('            - Interpolating the N Greenland results')
            theta_grid_N_Greenland = interpolate_model_output_onto_line(results_dir_N_Greenland, var_name, year, month, break_line)

            theta_grid = np.zeros_like(theta_grid_W_Greenland)

            W_indices = distance > 3000e3
            for d in range(np.shape(theta_grid_N_Greenland)[0]):
                theta_grid[d, W_indices] = theta_grid_W_Greenland[d, W_indices]

            N_indices = distance < 3000e3
            for d in range(np.shape(theta_grid_N_Greenland)[0]):
                theta_grid[d, N_indices] = theta_grid_N_Greenland[d, N_indices]

            year_output[month-1, :, :] = theta_grid

            # plt.subplot(2,1,1)
            # plt.pcolormesh(theta_grid)
            #
            # plt.subplot(2,1,2)
            # plt.plot(break_line[:,1])
            # plt.show()

        if np.any(year_output!=0):
            output_annual_nc_file(project_folder, var_name, year, break_line, distance, depth, year_output,
                                  drF, line_depth, line_hFacC)




