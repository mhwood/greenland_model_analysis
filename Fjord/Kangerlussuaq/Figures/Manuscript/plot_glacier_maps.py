
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
from scipy.interpolate import griddata
import cmocean.cm as cm
import shapefile as sf
from pyproj import Transformer
from osgeo import gdal

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

def read_model_grid_from_nc(config_dir):

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids','L2_Kanger_grid.nc'))
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    ds.close()

    points = reproject_polygon(np.column_stack([XC.ravel(), YC.ravel()]), 4326, 3413)
    X = np.reshape(points[:, 0], np.shape(XC))
    Y = np.reshape(points[:, 1], np.shape(YC))

    return(XC, YC, X, Y, Depth)

def read_L2_CTD_locations(config_dir, L2_model_name, L2_XC, L2_YC):

    dv_file = os.path.join(config_dir,'L2',L2_model_name,'input','dv','CTD_mask.bin')

    dv_grid = np.fromfile(dv_file,'>f4').reshape(np.shape(L2_XC))

    rows, cols = np.where(dv_grid!=0)

    # locations = np.zeros((len(rows),2))
    # for i in range(len(rows)):
    #     locations[i,0] = L2_XC[rows[i],cols[i]]
    #     locations[i,1] = L2_YC[rows[i],cols[i]]

    return(rows,cols)

def read_background_imagery(file_path):
    ds = gdal.Open(file_path)
    R = np.array(ds.GetRasterBand(1).ReadAsArray())
    G = np.array(ds.GetRasterBand(2).ReadAsArray())
    B = np.array(ds.GetRasterBand(3).ReadAsArray())
    rows = np.shape(R)[0]
    cols = np.shape(R)[1]
    R = R.reshape((rows,cols, 1))
    G = G.reshape((rows, cols, 1))
    B = B.reshape((rows, cols, 1))
    image = np.concatenate([B,G,R],axis=2)
    brightness_factor = 0.1 # 0 - 1
    image = (np.max(image)-np.max(image)*(brightness_factor))/(np.max(image)-np.min(image))*(image-np.min(image))+np.max(image)*(brightness_factor)
    image[image<0]=0
    image[image>1]=1
    transform = ds.GetGeoTransform()
    extents = [transform[0],transform[0]+transform[1]*np.shape(image)[1],transform[3]+ transform[5] * np.shape(image)[0], transform[3]]
    # x_resolution = transform[1]
    # y_resolution = transform[5]
    return(image,extents)

def read_landsat_subset_from_nc(file_path):
    ds = nc4.Dataset(file_path)
    X = ds.variables['X'][:, :]
    Y = ds.variables['Y'][:, :]
    R = ds.variables['band_4'][:, :]
    G = ds.variables['band_3'][:, :]
    B = ds.variables['band_2'][:, :]
    ds.close()
    rows = np.shape(X)[0]
    cols = np.shape(Y)[1]

    R = R.reshape((rows, cols, 1))/65535
    G = G.reshape((rows, cols, 1))/65535
    B = B.reshape((rows, cols, 1))/65535
    image = np.concatenate([R, G, B], axis=2)
    brightness_factor = 0.2  # 0 - 1
    image = (np.max(image) - np.max(image) * (brightness_factor)) / (np.max(image) - np.min(image)) * (
                image - np.min(image)) + np.max(image) * (brightness_factor)

    # image[image < 0] = 0
    # image[image > 1] = 1

    # plt.imshow(image)
    # plt.show()
    return(X, Y, image)

def read_crosssectional_profile_locations(project_dir):

    file_path = os.path.join(project_dir,'Data','Ocean','Modeling','L2_Kanger_fjord_crosssection_profiles_Heat_Flux.nc')

    ds = nc4.Dataset(file_path)

    transect_1_lon = np.array(ds.groups['transect_1'].variables['XC'])
    transect_1_lat = np.array(ds.groups['transect_1'].variables['YC'])
    transect_1 = reproject_polygon(np.column_stack([transect_1_lon, transect_1_lat]),4326,3413)

    transect_2_lon = np.array(ds.groups['transect_2'].variables['XC'])
    transect_2_lat = np.array(ds.groups['transect_2'].variables['YC'])
    transect_2 = reproject_polygon(np.column_stack([transect_2_lon, transect_2_lat]), 4326, 3413)

    transect_3_lon = np.array(ds.groups['transect_3'].variables['XC'])
    transect_3_lat = np.array(ds.groups['transect_3'].variables['YC'])
    transect_3 = reproject_polygon(np.column_stack([transect_3_lon, transect_3_lat]), 4326, 3413)

    ds.close()

    return(transect_1, transect_2, transect_3)

def read_transect_points_from_shp(config_dir):
    shp_file = os.path.join(config_dir,'L2','L2_Kanger','input','dv_shp','Kanger_Trough')
    r = sf.Reader(shp_file)
    shapes = r.shapes()
    points = np.array(shapes[0].points)

    total_dist = 0
    for i in range(len(points)-1):
        total_dist += great_circle_distance(points[i,0], points[i,1],
                                            points[i+1,0], points[i+1, 1])
    N = int(total_dist/100)
    # print('        - Subsampling the transect to '+str(N)+' points')
    points = series_to_N_points(points,N)
    points = np.flipud(points)

    XC_boundary = points[:, 0]
    YC_boundary = points[:, 1]
    distance = np.zeros((len(XC_boundary),))
    for i in range(len(YC_boundary) - 1):
        dist = great_circle_distance(XC_boundary[i], YC_boundary[i], XC_boundary[i + 1], YC_boundary[i + 1])
        distance[i + 1] = distance[i] + dist

    distance_indices = distance<100e3
    points = points[distance_indices]
    distance = distance[distance_indices]

    points = reproject_polygon(points,4326,3413)

    return(points, distance)

def create_nested_figure(project_dir, XC, YC, X, Y, Depth,
                         CTD_rows, CTD_cols,
                         modis_image, modis_extents,
                         landsat_X, landsat_Y, landsat_image,
                         transect_1, transect_2, transect_3,
                         along_fjord_transect):

    dmin = 0
    dmax = 1001

    print('        - Reprojecting bounding box to rows and cols')
    ul_dist = ((X - landsat_X[-1, 0]) ** 2 + (Y - landsat_Y[-1, 0]) ** 2) ** 0.5
    ul_row, ul_col = np.where(ul_dist == np.min(ul_dist))

    ur_dist = ((X - landsat_X[-1, -1]) ** 2 + (Y - landsat_Y[-1, -1]) ** 2) ** 0.5
    ur_row, ur_col = np.where(ur_dist == np.min(ur_dist))

    ll_dist = ((X - landsat_X[0, 0]) ** 2 + (Y - landsat_Y[0, 0]) ** 2) ** 0.5
    ll_row, ll_col = np.where(ll_dist == np.min(ll_dist))

    lr_dist = ((X - landsat_X[0, -1]) ** 2 + (Y - landsat_Y[0, -1]) ** 2) ** 0.5
    lr_row, lr_col = np.where(lr_dist == np.min(lr_dist))

    # plt.imshow(X)
    # plt.colorbar()
    # plt.show()

    fig = plt.figure(figsize=(12,6.5))

    gs = GridSpec(14,20, left=0.05, right=0.93, bottom=0.08, top=0.95)

    ax1 = fig.add_subplot(gs[:,:10])
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax1.imshow(np.flipud(modis_image),  alpha=0.7)#extent=extents,
    plot_depth = np.ma.masked_where(Depth == 0, Depth)
    ax1.imshow(plot_depth, cmap=cm.deep, origin='lower', vmin=dmin, vmax=dmax)
    ax1.plot([ll_col, lr_col], [ll_row, lr_row], 'k-')
    ax1.plot([ll_col, ul_col], [ll_row, ul_row], 'k-')
    ax1.plot([ul_col, ur_col], [ul_row, ur_row], 'k-')
    ax1.plot([lr_col, ur_col], [lr_row, ur_row], 'k-')
    ax1.plot(CTD_cols, CTD_rows, 'ko', markersize=7)
    ax1.plot(CTD_cols, CTD_rows, 'wo', markersize=6)
    ax1.plot(320, 30, 'ko', markersize=7)
    ax1.plot(320, 30, 'wo', markersize=6)
    ax1.text(330, 30, '=OMG Survey Sites', ha='left', va='center', fontsize=14)
    ax1.set_ylabel('Model Rows',fontsize=14)
    ax1.set_xlabel('Model Columns',fontsize=14)
    ax1.text(10,10,'a)',ha='left',va='bottom', fontsize=14,
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
    ax1.text(50, 600, 'Kangerlussuaq Glacier', fontsize=14)
    ax1.arrow(90,590,35,-25,color='k',head_width=14)

    ax2 = fig.add_subplot(gs[:10, 10:])
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.imshow(np.flipud(landsat_image),
                 extent=[np.min(landsat_X)/1000, np.max(landsat_X)/1000, np.min(landsat_Y)/1000, np.max(landsat_Y)/1000])
    ax2.plot(along_fjord_transect[:, 0] / 1000, along_fjord_transect[:, 1] / 1000, '-', color='black',linewidth=4)
    ax2.plot(along_fjord_transect[:, 0] / 1000, along_fjord_transect[:, 1] / 1000, '-', color='orange',linewidth=3)
    ax2.plot(transect_1[:, 0]/1000, transect_1[:, 1]/1000, '-', color='black',linewidth=4)
    ax2.plot(transect_2[:, 0]/1000, transect_2[:, 1]/1000, '-', color='black',linewidth=4)
    ax2.plot(transect_3[:, 0]/1000, transect_3[:, 1]/1000, '-', color='black',linewidth=4)
    ax2.plot(transect_1[:, 0] / 1000, transect_1[:, 1] / 1000, '-', color='white', linewidth=3)
    ax2.plot(transect_2[:, 0] / 1000, transect_2[:, 1] / 1000, '-', color='white', linewidth=3)
    ax2.plot(transect_3[:, 0] / 1000, transect_3[:, 1] / 1000, '-', color='white', linewidth=3)
    ax2.text(489.5, -2352.5, 'b)', ha='left', va='bottom', fontsize=14,
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel('Distance North (km, EPSG: 3413)',fontsize=14)
    ax2.set_xlabel('Distance East (km, EPSG: 3413)',fontsize=14)

    ax3 = fig.add_subplot(gs[-3:-1, 10:])
    ax3.plot([0.5, 1], [1, 1], '-', color='black',linewidth=4)
    ax3.plot([0.5,1],[1,1],'-', color='orange',linewidth=3)
    ax3.text(1.1,1,'Along-Fjord Transect\n(Figure 4)',ha='left',va='center',fontsize=14)
    ax3.plot([5.5, 6], [1, 1], '-', color='black',linewidth=4)
    ax3.plot([5.5, 6], [1, 1], '-', color='white',linewidth=3)
    ax3.text(6.1, 1, 'Cross-Fjord Transects\n(Figure 5)', ha='left', va='center',fontsize=14)
    ax3.set_xlim([0,10])
    ax3.set_ylim([0,2])
    ax3.axis('off')

    ax4 = fig.add_subplot(gs[-1, 11:])
    ax4.tick_params(axis='both', which='major', labelsize=12)
    cx = np.arange(dmin,dmax)
    cy = np.array([0, 1])
    CX, CY = np.meshgrid(cx, cy)
    ax4.pcolormesh(CX, CY, CX, cmap=cm.deep)
    ax4.set_yticks([])
    ax4.set_xlabel('Depth (m)',fontsize=14)

    output_file = os.path.join(project_dir,'Figures','Glacier','Glacier Maps.png')
    plt.savefig(output_file)
    plt.close(fig)

config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/' \
             'MITgcm/configurations/downscaled_greenland/'

project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

print('    - Reading in the model grid')
XC, YC, X, Y, Depth = read_model_grid_from_nc(config_dir)

print('    - Reading in the CTD rows and columns')
CTD_rows, CTD_cols = read_L2_CTD_locations(config_dir, 'L2_Kanger', XC, YC)

print('    - Reading in the modis imagery')
modis_image_path = os.path.join(config_dir,'L2','L2_Kanger','plots','basemap','L2_Kanger_MODIS_20220720_row_col.tif')
modis_image, modis_extents = read_background_imagery(modis_image_path)

print('    - Reading in the landsat imagery')
landsat_image_path = os.path.join(project_dir,'Data','Glacier','Kangerlussuaq Fjord Landsat.nc')
landsat_X, landsat_Y, landsat_image = read_landsat_subset_from_nc(landsat_image_path)

print('    - Reading in the cross-sectional transects')
transect_1, transect_2, transect_3 = read_crosssectional_profile_locations(project_dir)

print('    - Reading the along-fjord transect')
along_fjord_transect, distance = read_transect_points_from_shp(config_dir)

print('    - Making the figure')
create_nested_figure(project_dir, XC, YC, X, Y, Depth,
                     CTD_rows, CTD_cols,
                     modis_image, modis_extents,
                     landsat_X, landsat_Y, landsat_image,
                     transect_1, transect_2, transect_3,
                     along_fjord_transect)

