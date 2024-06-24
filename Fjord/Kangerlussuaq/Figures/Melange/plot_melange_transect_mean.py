



import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from scipy.interpolate import griddata
import shapefile
from pyproj import Transformer

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


def read_melange_polygon_from_shapefile(project_folder):

    sf = shapefile.Reader(os.path.join(project_folder,'Map','Shapefiles','Kangerlussuaq Melange Area'))
    shapes = sf.shapes()
    polygon = np.array(shapes[0].points)
    sf.close()

    return(polygon)


def make_transect_from_polygon(polygon):

    split_index = int(np.shape(polygon)[0]/2)+3

    polygon_bottom = polygon[:split_index,:]
    polygon_top = polygon[split_index:-1, :]

    polygon_bottom = series_to_N_points(polygon_bottom, 500)
    polygon_top = series_to_N_points(polygon_top, 500)
    polygon_top = np.flipud(polygon_top)

    transect = np.zeros_like(polygon_top).astype(float)
    transect[:, 0] = (polygon_top[:, 0] + polygon_bottom[:, 0]) / 2
    transect[:, 1] = (polygon_top[:, 1] + polygon_bottom[:, 1]) / 2

    # plt.plot(polygon_bottom[:,0], polygon_bottom[:,1], 'b-')
    # plt.plot(polygon_top[:, 0], polygon_top[:, 1], 'k-')
    # plt.plot(transect[:, 0], transect[:, 1], 'g-')
    # plt.show()

    return(transect)


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


def read_model_melange(config_dir, project_folder):

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids','L2_Kanger_grid.nc'))
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    # Depth = ds.variables['Depth'][:, :]
    ds.close()

    elevation = np.fromfile(os.path.join(project_folder,'Data','Modeling','Melange_Elevation.bin'), '>f4')
    elevation = np.reshape(elevation, np.shape(XC))
    elevation[elevation==0] = np.nan

    points = reproject_polygon(np.column_stack([XC.ravel(), YC.ravel()]), 4326, 3413)
    X = np.reshape(points[:, 0], np.shape(XC))
    Y = np.reshape(points[:, 1], np.shape(YC))


    return(X, Y, elevation)


def sample_melange_on_transect(transect, elev_dict, model_X, model_Y, model_melange):

    transects = []

    for year in ['2016','2017','2018','2019']:
        X = elev_dict[year][0]
        Y = elev_dict[year][1]
        elev = elev_dict[year][2]

        profile = griddata(np.column_stack([X.ravel(), Y.ravel()]), elev.ravel(),
                           (transect[:,0], transect[:,1]), method='nearest', fill_value=0)
        # print(np.shape(profile))
        transects.append(profile)

    model_profile = griddata(np.column_stack([model_X.ravel(), model_Y.ravel()]), model_melange.ravel(),
                       (transect[:, 0], transect[:, 1]), method='nearest', fill_value=0)
    # print(np.shape(profile))
    transects.append(model_profile)

    return(transects)


def plot_melange_transects(project_folder, transect, profiles):

    distance = np.zeros((np.shape(transect)[0],))
    for d in range(1,np.shape(transect)[0]):
        distance[d]= ((transect[d,0]-transect[d-1,0])**2 + (transect[d,1]-transect[d-1,1])**2)**0.5+distance[d-1]
    distance *= 1e-3
    distance = np.max(distance) - distance

    labels = ['2016','2017','2018','2019','Model']
    colors = ['red','orange','green','blue','black']

    fig = plt.figure(figsize=(10,5))

    for p in range(len(profiles)):
        if p == len(profiles)-1:
            plt.plot(distance, profiles[p], label=labels[p], color = 'k', linewidth=2)
        else:
            plt.plot(distance, profiles[p], label=labels[p], color=colors[p], linewidth=1)

    plt.legend(loc=1)
    plt.ylabel('Melange Elevation\n(relative to GOCO05C Geoid)')
    plt.xlabel('Distance Along Transect (km)')

    plt.grid(linestyle='--', linewidth=0.5, alpha=0.5)

    # plt.gca().set_ylim([-1,1])

    output_file = os.path.join(project_folder,'Figures','Glacier','Kangerlussuaq Melange Transects.png')
    plt.savefig(output_file)
    plt.close(fig)
    a=1


def plot_transect_comparison():
    project_folder = '/Users/mike/Documents/Research/Projects/Kangerlussuaq'

    config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
                 'configurations/downscaled_greenland'

    polygon = read_melange_polygon_from_shapefile(project_folder)

    transect = make_transect_from_polygon(polygon)

    elev_dict = read_melange_data_from_nc(project_folder)

    model_X, model_Y, model_melange = read_model_melange(config_dir, project_folder)

    profiles = sample_melange_on_transect(transect, elev_dict, model_X, model_Y, model_melange)

    plot_melange_transects(project_folder, transect, profiles)
