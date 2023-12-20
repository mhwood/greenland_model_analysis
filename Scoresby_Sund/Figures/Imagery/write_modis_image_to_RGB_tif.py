
import os
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc4
from scipy.interpolate import griddata
from osgeo import gdal
from osgeo import osr
from pyproj import Transformer

def read_extent_from_model_grid_nc(config_dir,model_name):

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids',model_name+'_grid.nc'))
    Lon = ds.variables['XC'][:, :]
    Lat = ds.variables['YC'][:, :]
    ds.close()

    return(Lon, Lat)

def reproject_points(points,inputCRS,outputCRS,x_column=0,y_column=1):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(points[:, x_column], points[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
        run_test = False
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon = np.copy(points)
    output_polygon[:, x_column] = x2
    output_polygon[:, y_column] = y2
    return output_polygon

def read_RGB_bands_from_nc(nc_file_path):

    ds = nc4.Dataset(nc_file_path)
    Lon = ds.variables['Longitude'][:, :]
    Lat = ds.variables['Latitude'][:, :]
    band_1 = ds.variables['sur_refl_b01'][:, :, :]
    band_1 = band_1[0, :, :]
    band_3 = ds.variables['sur_refl_b03'][:, :, :]
    band_3 = band_3[0, :, :]
    band_4 = ds.variables['sur_refl_b04'][:, :, :]
    band_4 = band_4[0, :, :]

    ds.close()

    return(Lon, Lat, band_1, band_3, band_4)

def read_MODIS_points_to_domain(modis_path,min_lon,max_lon,min_lat,max_lat,reproject_to_polar):

    file_names = ['h16v01.20220720.ncml.nc4','h17v01.20220720.ncml.nc4',
                  'h16v02.20220720.ncml.nc4','h17v02.20220720.ncml.nc4']

    resolution_buffer = 0.1

    points_started = False

    for file_name in file_names:
        print('   - Reading '+file_name)
        nc_file_path = modis_path+'/'+file_name

        Lon, Lat, band_1, band_3, band_4 = read_RGB_bands_from_nc(nc_file_path)

        points = np.column_stack([np.ravel(Lon), np.ravel(Lat)])
        values_1 = np.reshape(band_1,(np.size(band_1),1)) # red
        values_3 = np.reshape(band_3,(np.size(band_3),1)) # green
        values_4 = np.reshape(band_4,(np.size(band_4),1)) # blue

        indices_lon = np.logical_and(points[:,0]>=min_lon-resolution_buffer,
                                     points[:,0]<=max_lon+resolution_buffer)
        indices_lat = np.logical_and(points[:,1]>=min_lat-resolution_buffer,
                                     points[:,1]<=max_lat+resolution_buffer)
        indices = np.logical_and(indices_lon, indices_lat)

        if np.any(indices):

            points_subset = points[indices, :]
            band_1_subset = values_1[indices]
            band_3_subset = values_3[indices]
            band_4_subset = values_4[indices]

            if not points_started:
                points_started = True
                all_points = points_subset
                all_band_1_points = band_1_subset
                all_band_3_points = band_3_subset
                all_band_4_points = band_4_subset
            else:
                all_points = np.vstack([all_points,points_subset])
                all_band_1_points = np.vstack([all_band_1_points,band_1_subset])
                all_band_3_points = np.vstack([all_band_3_points,band_3_subset])
                all_band_4_points = np.vstack([all_band_4_points,band_4_subset])

    if reproject_to_polar:
        all_points = reproject_points(all_points, inputCRS=4326, outputCRS=3413)

    return(all_points,all_band_1_points,all_band_3_points,all_band_4_points)

def interpolate_points_to_grid(points, X, Y, band_1_points, band_3_points, band_4_points):

    print(np.shape(X))

    print('    - Interpolating band 1')
    band_1_grid = griddata(points,band_1_points.ravel(),(X,Y))
    if np.any(np.isnan(band_1_grid)):
        band_1_nearest = griddata(points,band_1_points.ravel(),(X,Y), method='nearest')
        band_1_grid[np.isnan(band_1_grid)] = band_1_nearest[np.isnan(band_1_grid)]

    print('    - Interpolating band 3')
    band_3_grid = griddata(points,band_3_points.ravel(),(X,Y))
    if np.any(np.isnan(band_3_grid)):
        band_3_nearest = griddata(points,band_3_points.ravel(),(X,Y), method='nearest')
        band_3_grid[np.isnan(band_3_grid)] = band_3_nearest[np.isnan(band_3_grid)]

    print('    - Interpolating band 4')
    band_4_grid = griddata(points,band_4_points.ravel(),(X,Y))
    if np.any(np.isnan(band_4_grid)):
        band_4_nearest = griddata(points,band_4_points.ravel(),(X,Y), method='nearest')
        band_4_grid[np.isnan(band_4_grid)] = band_4_nearest[np.isnan(band_4_grid)]

    return(band_1_grid, band_3_grid, band_4_grid)

def write_data_to_tif(output_file, epsg, x,y, band_1_grid, band_3_grid, band_4_grid):

    geotransform = (np.min(x), x[1]-x[0], 0, np.max(y), 0, y[0]-y[1])

    output_raster = gdal.GetDriverByName('GTiff').Create(output_file, len(x), len(y), 3,
                                                         gdal.GDT_Float32)  # Open the file
    output_raster.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)

    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(np.flipud(band_3_grid))  # Writes my array to the raster
    output_raster.GetRasterBand(2).WriteArray(np.flipud(band_4_grid))  # Writes my array to the raster
    output_raster.GetRasterBand(3).WriteArray(np.flipud(band_1_grid))  # Writes my array to the raster

    output_raster.FlushCache()

config_dir = '/Volumes/ifenty/Wood/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/downscaled_greenland'

level_name = 'L3'
model_name = 'L3_Scoresby_Sund'

# level_name = 'L2'
# model_name = 'L2_CE_Greenland'
#
# level_name = 'L1'
# model_name = 'L1_CE_Greenland'

buffer = 0.1

reproject_to_polar = True
resolution = 250

# reproject_to_polar = False
# resolution = 0.1

# step 1: read in the geometry
print(' - Reading in the model geometry')
Lon, Lat = read_extent_from_model_grid_nc(config_dir,model_name)

min_lon = np.min(Lon)
max_lon = np.max(Lon)
min_lat = np.min(Lat)
max_lat = np.max(Lat)

# step 2: reproject to polar coorindates
print(' - Creating the output grids')
if reproject_to_polar:
    points = np.column_stack([Lon.ravel(),Lat.ravel()])
    points = reproject_points(points,inputCRS=4326,outputCRS=3413)
    if model_name=='L3_Scoresby_Sund':
        extra_shift = 50000
        x = np.arange(np.min(points[:,0]) - resolution - extra_shift, np.max(points[:,0]) + 2 * resolution, resolution)
    else:
        x = np.arange(np.min(points[:, 0]) - resolution, np.max(points[:, 0]) + 2 * resolution, resolution)
    y = np.arange(np.min(points[:,1]) - resolution, np.max(points[:,1]) + 2 * resolution, resolution)
    X, Y = np.meshgrid(x,y)
    epsg = 3413
else:
    x = np.arange(min_lon - resolution, max_lon + 2 * resolution, resolution)
    y = np.arange(min_lat - resolution, max_lat + 2 * resolution, resolution)
    X, Y = np.meshgrid(x,y)
    epsg = 4326

# step 3: read in the modis points
print(' - Reading in the MODIS data')
modis_path = '/Users/michwood/Documents/Research/Projects/Scoresby Sund/Data/Remote Sensing/Raw/MODIS'
points, band_1_points, band_3_points, band_4_points = \
    read_MODIS_points_to_domain(modis_path,min_lon,max_lon,min_lat,max_lat,reproject_to_polar)

# step 4: interpolate the modis points onto the grid
print(' - Interpolating the points onto the domain')
band_1_grid, band_3_grid, band_4_grid = \
    interpolate_points_to_grid(points, X, Y, band_1_points, band_3_points, band_4_points)
print(np.min(band_1_grid),np.max(band_1_grid))
print(np.min(band_3_grid),np.max(band_3_grid))
print(np.min(band_4_grid),np.max(band_4_grid))

# step 5: output the files to tif
print(' - Outputting the bands to tif')
output_file = '/Users/michwood/Documents/Research/Projects/Scoresby Sund/Data/Remote Sensing/' \
              'Processed/MODIS/'+model_name+'_MODIS_20220720_'+str(epsg)+'.tif'
write_data_to_tif(output_file, epsg, x, y, band_1_grid, band_3_grid, band_4_grid)

