
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import griddata
from osgeo import gdal
import shapefile
from pyproj import Transformer

def read_melange_polygon_from_shapefile(project_folder):

    sf = shapefile.Reader(os.path.join(project_folder,'Map','Shapefiles','Kangerlussuaq Melange Area'))
    shapes = sf.shapes()
    polygon = np.array(shapes[0].points)
    sf.close()

    return(polygon)


def read_landsat_image(file_path):

    ds = gdal.Open(file_path)
    R = np.array(ds.GetRasterBand(1).ReadAsArray())
    rows = np.shape(R)[0]
    cols = np.shape(R)[1]

    transform = ds.GetGeoTransform()
    extents = [transform[0], transform[0] + transform[1] * cols,
               transform[3] + transform[5] * rows, transform[3]]

    ds = None

    return(R, extents)



    # x_resolution = transform[1]
    # y_resolution = transform[5]


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
    elif inputCRS == 3413 and str(outputCRS)[:3] == '326':
        x2, y2 = transformer.transform(polygon_array[:,x_column], polygon_array[:,y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon


def subset_landsat_to_new_grid(grid, extent, dest_X_32625, dest_Y_32625):

    src_x = np.arange(extent[0], extent[1], 30)
    src_y = np.arange(extent[2], extent[3], 30)
    src_y = np.flip(src_y)
    src_X, src_Y = np.meshgrid(src_x, src_y)

    # print(np.shape(src_X), np.shape(src_Y), np.shape(grid))

    ll_dist = (src_X-np.min(dest_X_32625))**2 + (src_Y-np.min(dest_Y_32625))**2
    ll_row, ll_col = np.where(ll_dist==np.min(ll_dist))
    ll_row = ll_row[0]
    ll_col = ll_col[0]

    ul_dist = (src_X - np.max(dest_X_32625)) ** 2 + (src_Y -  np.min(dest_Y_32625)) ** 2
    ul_row, ul_col = np.where(ul_dist == np.min(ul_dist))
    ul_row = ul_row[0]
    ul_col = ul_col[0]

    ur_dist = (src_X - np.min(dest_X_32625)) ** 2 + (src_Y - np.max(dest_Y_32625)) ** 2
    ur_row, ur_col = np.where(ur_dist == np.min(ur_dist))
    ur_row = ur_row[0]
    ur_col = ur_col[0]

    lr_dist = (src_X - np.max(dest_X_32625)) ** 2 + (src_Y - np.max(dest_Y_32625)) ** 2
    lr_row, lr_col = np.where(lr_dist == np.min(lr_dist))
    lr_row = lr_row[0]
    lr_col = lr_col[0]

    min_row = np.min([ll_row, lr_row, ur_row, ul_row])-100
    max_row = np.max([ll_row, lr_row, ur_row, ul_row])+100

    min_col = np.min([ll_col, lr_col, ur_col, ul_col])-100
    max_col = np.max([ll_col, lr_col, ur_col, ul_col])+100

    print(ll_row, lr_row, ur_row, ul_row)
    print(min_row, max_row)
    print(ll_col, lr_col, ur_col, ul_col)
    print(min_col, max_col)

    src_X = src_X[min_row:max_row, :]
    src_Y = src_Y[min_row:max_row, :]
    grid = grid[min_row:max_row, :]

    src_X = src_X[:, min_col:max_col]
    src_Y = src_Y[:, min_col:max_col]
    grid = grid[:, min_col:max_col]

    print(np.shape(src_X), np.shape(src_Y), np.shape(grid))
    print(np.shape(dest_X_32625), np.shape(dest_Y_32625))

    interp_grid = griddata(np.column_stack([src_X.ravel(), src_Y.ravel()]),
                           grid.ravel(), (dest_X_32625, dest_Y_32625),
                           method='nearest')

    return(interp_grid)

def stack_landsat_subsets_to_nc(output_file, X, Y, grids):

    ds = nc4.Dataset(output_file, 'w')

    ds.createDimension('rows',np.shape(X)[0])
    ds.createDimension('cols',np.shape(Y)[1])

    x = ds.createVariable('X','f4',('rows','cols'))
    x[:, :] = X

    x = ds.createVariable('Y', 'f4', ('rows', 'cols'))
    x[:, :] = Y

    x = ds.createVariable('band_2', 'f4', ('rows', 'cols'))
    x[:, :] = grids[0]

    x = ds.createVariable('band_3', 'f4', ('rows', 'cols'))
    x[:, :] = grids[1]

    x = ds.createVariable('band_4', 'f4', ('rows', 'cols'))
    x[:, :] = grids[2]

    ds.close()


project_folder = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

subset = 'Fjord'

if subset=='Melange':
    polygon = read_melange_polygon_from_shapefile(project_folder)
    dest_x = np.arange(np.min(polygon[:, 0]), np.max(polygon[:, 0]), 50)
    dest_y = np.arange(np.min(polygon[:, 1]), np.max(polygon[:, 1]), 50)

if subset=='Fjord':
    dest_x = np.arange(487800, 577500, 100)
    dest_y = np.arange(-2353900, -2283500, 100)

dest_X, dest_Y = np.meshgrid(dest_x, dest_y)
points = reproject_polygon(np.column_stack([dest_X.ravel(), dest_Y.ravel()]), 3413, 32625)
dest_X_32625 = np.reshape(points[:,0], np.shape(dest_X))
dest_Y_32625 = np.reshape(points[:,1], np.shape(dest_Y))

# plt.plot(dest_X_32625[:,0], dest_Y_32625[:,0])
# plt.plot(dest_X_32625[:,-1], dest_Y_32625[:,-1])
# plt.plot(dest_X_32625[0, :], dest_Y_32625[0, :])
# plt.plot(dest_X_32625[-1, :], dest_Y_32625[-1, :])
# plt.show()

file_names = ['LC08_L2SP_230012_20160720_20200906_02_T1_SR_B2.TIF',
              'LC08_L2SP_230012_20160720_20200906_02_T1_SR_B3.TIF',
              'LC08_L2SP_230012_20160720_20200906_02_T1_SR_B4.TIF']

grids = []

for file_name in file_names:

    file_path = '/Users/mike/Documents/Research/Data Repository/Greenland/Landsat/'+file_name

    grid, extent = read_landsat_image(file_path)

    # plt.plot(dest_X_32625[:,0], dest_Y_32625[:,0])
    # plt.plot(dest_X_32625[:,-1], dest_Y_32625[:,-1])
    # plt.plot(dest_X_32625[0, :], dest_Y_32625[0, :])
    # plt.plot(dest_X_32625[-1, :], dest_Y_32625[-1, :])
    #
    # plt.plot([extent[0], extent[1]], [extent[2], extent[2]], 'k-')
    # plt.plot([extent[0], extent[1]], [extent[3], extent[3]], 'k-')
    # plt.plot([extent[0], extent[0]], [extent[2], extent[3]], 'k-')
    # plt.plot([extent[1], extent[1]], [extent[2], extent[3]], 'k-')
    # plt.show()

    interp_grid = subset_landsat_to_new_grid(grid, extent, dest_X_32625, dest_Y_32625)

    grids.append(interp_grid)

    # plt.imshow(interp_grid, origin='lower')
    # plt.show()

output_file = project_folder+'/Data/Glacier/Kangerlussuaq '+subset+' Landsat.nc'
stack_landsat_subsets_to_nc(output_file, dest_X, dest_Y, grids)








