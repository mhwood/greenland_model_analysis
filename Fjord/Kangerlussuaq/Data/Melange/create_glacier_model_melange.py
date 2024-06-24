

import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import shapefile
from pyproj import Transformer
from scipy.interpolate import griddata

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1, run_test = True):

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


def read_melange_polygon_from_shapefile(project_folder):

    sf = shapefile.Reader(os.path.join(project_folder,'Map','Shapefiles','Kangerlussuaq Melange Area'))
    shapes = sf.shapes()
    polygon = np.array(shapes[0].points)
    sf.close()

    return(polygon)


def read_melange_data_from_nc(project_folder):

    years = ['2016','2017','2018','2019']
    file_name = os.path.join(project_folder,'Data','Glacier','Kangerlussuaw Melange GLISTIN Elevation.nc')

    elev_dict = {}

    ds = nc4.Dataset(file_name)
    for year in years:
        grp = ds.groups[year]
        X = grp.variables['X'][:, :]
        Y = grp.variables['Y'][:, :]
        elevation = np.array(grp.variables['elevation'][:, :])
        polygon = np.column_stack([grp.variables['polygon_x'][:], grp.variables['polygon_y'][:]])

        p = mplPath.Path(polygon)
        inside = p.contains_points(np.column_stack([X.ravel(), Y.ravel()]))
        inside = np.reshape(inside, np.shape(elevation))
        elevation[elevation==-9999] = np.nan
        elevation[~inside] = np.nan
        elevation[np.isfinite(elevation)]+=2
        # print(np.min(elevation[~np.isnan(elevation)]), np.max(elevation[~np.isnan(elevation)]))
        # print(np.sum(elevation<0), np.median(elevation[np.logical_and(np.isfinite(elevation),elevation<0)]))

        elev_dict[year] = [X, Y, elevation, polygon]
    ds.close()

    return(elev_dict)


def read_model_geometry(config_dir):

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids','L2_Kanger_grid.nc'))
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    ds.close()

    return(XC,YC,Depth)


def create_model_melange_file(polygon, elev_dict, XC, YC, Depth):

    points = reproject_polygon(np.column_stack([XC.ravel(), YC.ravel()]), 4326, 3413)
    X = np.reshape(points[:,0], np.shape(XC))
    Y = np.reshape(points[:,1], np.shape(YC))

    p = mplPath.Path(polygon)
    inside = p.contains_points(np.column_stack([X.ravel(), Y.ravel()]))
    inside = np.reshape(inside, np.shape(XC))

    sum_grid = np.zeros(np.shape(Depth))
    count_grid = np.zeros(np.shape(Depth))
    for year in ['2016', '2017', '2018', '2019']:
        glistin_X = elev_dict[year][0]
        glistin_Y = elev_dict[year][1]
        glistin_elev = elev_dict[year][2]
        points = np.column_stack([glistin_X.ravel(), glistin_Y.ravel()])
        values = glistin_elev.ravel()
        indices = np.logical_and(np.isfinite(values), values>-5)
        points = points[indices, :]
        values = values[indices]
        interp_glistin_elev = griddata(points, values, (X, Y), method='nearest')

        for row in range(np.shape(sum_grid)[0]):
            for col in range(np.shape(sum_grid)[1]):
                if inside[row, col] and Depth[row,col]>0:
                    sum_grid[row, col] += interp_glistin_elev[row, col]
                    count_grid[row,col] += 1

    elevation = np.copy(sum_grid)
    elevation[count_grid>0] = sum_grid[count_grid>0]/count_grid[count_grid>0]

    plt.subplot(1,2,1)
    plt.imshow(sum_grid[500:600,100:200])
    plt.subplot(1, 2, 2)
    plt.imshow(count_grid[500:600, 100:200])
    plt.show()

    return(elevation)



project_folder = '/Users/mike/Documents/Research/Projects/Kangerlussuaq'


config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
             'configurations/downscaled_greenland'


polygon = read_melange_polygon_from_shapefile(project_folder)

elev_dict = read_melange_data_from_nc(project_folder)

XC, YC, Depth = read_model_geometry(config_dir)

elevation = create_model_melange_file(polygon, elev_dict, XC, YC, Depth)

# plt.imshow(elevation)
# plt.gca().set_xlim([100,200])
# plt.gca().set_ylim([500,600])
# plt.show()

output_file = os.path.join(project_folder,'Data','Modeling','Melange_Elevation.bin')
elevation.ravel('C').astype('>f4').tofile(output_file)







