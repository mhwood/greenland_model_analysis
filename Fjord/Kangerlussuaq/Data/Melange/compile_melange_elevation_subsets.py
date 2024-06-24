

import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import shapefile
from pyproj import Transformer
import xarray as xr

def read_melange_polygon_from_shapefile(project_folder):

    sf = shapefile.Reader(os.path.join(project_folder,'Map','Shapefiles','Kangerlussuaq Melange Area'))
    shapes = sf.shapes()
    polygon = np.array(shapes[0].points)
    sf.close()

    return(polygon)

def read_glistin_fronts_from_shapefile(project_folder):

    sf = shapefile.Reader(os.path.join(project_folder,'Map','Shapefiles','Kangerlussuaq GLISTIN Ice Fronts'))
    shapes = sf.shapes()
    records = sf.records()

    for i in range(len(records)):
        if records[i][0] == '2016':
            front_2016 = np.array(shapes[i].points)
        if records[i][0] == '2017':
            front_2017 = np.array(shapes[i].points)
        if records[i][0] == '2018':
            front_2018 = np.array(shapes[i].points)
        if records[i][0] == '2019':
            front_2019 = np.array(shapes[i].points)

    sf.close()

    return([front_2016, front_2017, front_2018, front_2019])

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

def create_melange_polygons(polygon, fronts):

    polygon = series_to_N_points(polygon, 1000)

    melange_polygons = []

    for front in fronts:

        dist_1 = (polygon[:,0]-front[0,0])**2 + (polygon[:,1]-front[0,1])**2
        index_1 = np.argmin(dist_1)

        dist_2 = (polygon[:, 0] - front[-1, 0]) ** 2 + (polygon[:, 1] - front[-1, 1]) ** 2
        index_2 = np.argmin(dist_2)

        melange_polygon = np.vstack([polygon[:index_1,:],
                                     front,
                                     polygon[index_2:,:]])

        # plt.plot(polygon[:, 0], polygon[:, 1], 'k-')
        # plt.plot(melange_polygon[:, 0], melange_polygon[:, 1], 'b-',linewidth=1)
        # plt.show()

        melange_polygons.append(melange_polygon)

    return(melange_polygons)

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
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon

def read_glistin_elevation_data(elevation_folder, polygon, year):

    if year==2016:
        file_name = 'OMG_Ice_GLISTIN-A_L3_20160323141236_20.nc'
        offset = 0.1
    if year==2017:
        file_name = 'OMG_Ice_GLISTIN-A_L3_20170313114240_20.nc'
        offset = 0.65
    if year==2018:
        file_name = 'OMG_Ice_GLISTIN-A_L3_20180310145744_20.nc'
        offset = -1.2
    if year==2019:
        file_name = 'OMG_Ice_GLISTIN-A_L3_20190314144048_20.nc'
        offset = -0.6

    ds = nc4.Dataset(os.path.join(elevation_folder,file_name))
    lon = ds.variables['longitude'][:, :]
    lat = ds.variables['latitude'][:, :]
    # elevation = np.array(ds.variables['elevation'][:, :]).astype(float)
    # geoid = np.array(ds.variables['geoid'][:, :]).astype(float)
    ds.close()
    # print(np.min(elevation), np.max(elevation), elevation[-1, 325])
    # print(type(elevation[0,0]))

    ds = xr.open_dataset(os.path.join(elevation_folder,file_name))
    elevation = ds.variables['elevation'].values
    geoid = ds.variables['geoid'].values
    ds.close()
    # print(np.min(elevation), np.max(elevation), elevation[-1, 325])
    # print(type(elevation[0, 0]))

    elevation+=offset

    elevation = elevation-geoid

    # elevation = elevation

    # plt.imshow(elevation)
    # plt.show()

    points = reproject_polygon(np.column_stack([lon.ravel(), lat.ravel()]), 4326, 3413)
    X = np.reshape(points[:,0], np.shape(lon))
    Y = np.reshape(points[:,1], np.shape(lat))

    min_col = np.argmin(np.abs(X[-1,:] - np.min(polygon[:,0]))) - 10
    max_col = np.argmin(np.abs(X[0, :] - np.max(polygon[:, 0]))) + 10
    max_row = np.argmin(np.abs(Y[:, -1] - np.min(polygon[:, 1]))) - 10
    min_row = np.argmin(np.abs(Y[:, 0] - np.max(polygon[:, 1]))) + 10

    # print(min_row, max_row)
    # print(min_col, max_col)
    # print(np.shape(elevation))
    X = X[min_row:max_row, :]
    Y = Y[min_row:max_row, :]
    elevation = elevation[min_row:max_row, :]

    X = X[:, min_col:max_col]
    Y = Y[:, min_col:max_col]
    elevation = elevation[:, min_col:max_col]

    elevation[np.abs(elevation)>1e5] = -9999
    elevation[np.isnan(elevation)] = -9999
    elevation[~np.isreal(elevation)] = -9999
    print(np.min(elevation), np.max(elevation), elevation[-1,325])

    return(X,Y,elevation)

def write_elevation_to_nc(project_folder, years, fronts, melange_polygons, Xs, Ys, elevations):

    output_file = os.path.join(project_folder,'Data','Glacier','Kangerlussuaw Melange GLISTIN Elevation.nc')

    ds = nc4.Dataset(output_file,'w')

    for yy in range(len(years)):
        grp = ds.createGroup(str(years[yy]))
        grp.createDimension('rows', np.shape(Xs[yy])[0])
        grp.createDimension('cols', np.shape(Ys[yy])[1])
        grp.createDimension('front_points', np.shape(fronts[yy])[0])
        grp.createDimension('polygon_points', np.shape(melange_polygons[yy])[0])
        v = grp.createVariable('X','f4',('rows','cols'))
        v[:, :] = Xs[yy]
        v = grp.createVariable('Y', 'f4', ('rows', 'cols'))
        v[:, :] = Ys[yy]
        v = grp.createVariable('elevation', 'f4', ('rows', 'cols'), fill_value = -9999)
        # v._FillValue = -9999
        v[:, :] = elevations[yy]
        v = grp.createVariable('front_x', 'f4', ('front_points',))
        v[:] = fronts[yy][:,0]
        v = grp.createVariable('front_y', 'f4', ('front_points', ))
        v[:] = fronts[yy][:,1]
        v = grp.createVariable('polygon_x', 'f4', ('polygon_points',))
        v[:] = melange_polygons[yy][:, 0]
        v = grp.createVariable('polygon_y', 'f4', ('polygon_points',))
        v[:] = melange_polygons[yy][:, 1]



    ds.close()

project_folder = '/Users/mike/Documents/Research/Projects/Kangerlussuaq'

elevation_folder = '/Users/mike/Documents/Research/Data Repository/Greenland/Elevation'

polygon = read_melange_polygon_from_shapefile(project_folder)

fronts = read_glistin_fronts_from_shapefile(project_folder)

melange_polygons = create_melange_polygons(polygon, fronts)

Xs = []
Ys = []
elevations = []
years = [2016, 2017, 2018, 2019]

for year in years:
    print('Reading data in year '+str(year))
    X, Y, elevation = read_glistin_elevation_data(elevation_folder, polygon, year)
    Xs.append(X)
    Ys.append(Y)
    elevations.append(elevation)

# plt.pcolormesh(X,Y,elevation)
# plt.plot(polygon[:,0], polygon[:,1])
# plt.show()

write_elevation_to_nc(project_folder, years, fronts, melange_polygons, Xs, Ys, elevations)






