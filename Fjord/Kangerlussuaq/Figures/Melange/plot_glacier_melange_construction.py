


import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from matplotlib.patches import Rectangle
import shapefile
from matplotlib.gridspec import GridSpec
from pyproj import Transformer


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

def read_model_melange(config_dir, project_folder):

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids','L2_Kanger_grid.nc'))
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    ds.close()

    elevation = np.fromfile(os.path.join(project_folder,'Data','Ocean','Modeling','Melange_Elevation.bin'), '>f4')
    elevation = np.reshape(elevation, np.shape(XC))
    elevation[elevation==0] = np.nan

    points = reproject_polygon(np.column_stack([XC.ravel(), YC.ravel()]), 4326, 3413)
    X = np.reshape(points[:, 0], np.shape(XC))
    Y = np.reshape(points[:, 1], np.shape(YC))


    return(X, Y, elevation)

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
    image = np.concatenate([B, G, R], axis=2)
    brightness_factor = 0.1  # 0 - 1
    image = (np.max(image) - np.max(image) * (brightness_factor)) / (np.max(image) - np.min(image)) * (
                image - np.min(image)) + np.max(image) * (brightness_factor)
    image[image < 0] = 0
    image[image > 1] = 1

    # plt.imshow(image)
    # plt.show()
    return(X, Y, image)

def plot_bar_chart(ax, elev, vmin, vmax):

    # bins1 = np.linspace(0.0, 0.25, 1)
    # bins2 = np.linspace(0.025, 2, 39)
    # bins = np.concatenate([bins1, bins2], axis=0)
    bins = np.linspace(-1, 2, 40)

    elev = elev[np.logical_and(np.isfinite(elev), elev > 0)]
    N, _, _ = plt.hist(np.log10(elev), bins=bins)#, density=True)
    ax.cla()
    ax.bar(bins[:-1], N, width=np.diff(bins), align='edge')

    for b in range(len(bins)-1):
        # print(10**((bins[b+1]+bins[b])/2), (np.clip(10**((bins[b+1]+bins[b])/2), vmin, vmax) - vmin) / (vmax-vmin))
        color= plt.cm.turbo((np.clip(10**((bins[b+1]+bins[b])/2), vmin, vmax) - vmin) / (vmax-vmin))
        rect = Rectangle((bins[b], 0), bins[b+1]-bins[b], N[b], facecolor=color, edgecolor='k')
        ax.add_patch(rect)

    return(N)

def plot_melange_with_distributions(project_folder, polygon, elev_dict, model_X, model_Y, model_elev,
                                    landsat_X, landsat_Y, landsat_image):

    output_file = os.path.join(project_folder,'Figures','Glacier','Kangerlussuaq Melange.png')

    image_rows = 6
    image_cols = 6

    vmin = 0
    vmax = 50

    fig = plt.figure(figsize=(7,10))

    gs1 = GridSpec(image_rows*6 + 5, image_cols*2 + 1,
                   left=0.1, right=0.95, bottom=0.05, top=0.95)

    ########################################################################
    # Observations

    ax11m = fig.add_subplot(gs1[:image_rows, :image_cols])
    ax11m.imshow(np.flipud(landsat_image),
                 extent=[np.min(landsat_X)/1000, np.max(landsat_X)/1000, np.min(landsat_Y)/1000, np.max(landsat_Y)/1000],
                 alpha=0.8)
    ax11m.plot(polygon[:,0]/1000, polygon[:,1]/1000, 'y-', label='Melange\nRegion')
    # ax11m.legend(loc=2, facecolor='white')
    ax11m.set_xlim([np.min(polygon[:, 0]/1000), np.max(polygon[:, 0]/1000)])
    ax11m.set_ylim([np.min(polygon[:, 1]/1000), np.max(polygon[:, 1]/1000)])
    ax11m.text(np.max(polygon[:, 0] / 1000) - 1, np.max(polygon[:, 1] / 1000) - 1, 'Landsat 8\n2016/07/20', ha='right', va='top',
               bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    # ax11m.set_ylabel('Distance North (km)')
    ax11m.set_title('Kangerlussuaq Melange')
    ax11m.set_xticks([])
    ax11m.text(487, -2299, 'a)', ha='left', va='bottom', fontsize=12,
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))

    ########################################################################
    # 2016 Melange

    ax21m = fig.add_subplot(gs1[1*image_rows+1:2*image_rows+1, :image_cols])
    ax21m.imshow(np.flipud(landsat_image),
                 extent=[np.min(landsat_X)/1000, np.max(landsat_X)/1000, np.min(landsat_Y)/1000, np.max(landsat_Y)/1000],
                 alpha=0.8)
    ax21m.pcolormesh(elev_dict['2016'][0]/1000, elev_dict['2016'][1]/1000, elev_dict['2016'][2],
                     cmap='turbo', vmin=0, vmax=50)
    ax21m.plot(polygon[:, 0]/1000, polygon[:, 1]/1000, 'y-')
    ax21m.set_xlim([np.min(polygon[:, 0]/1000), np.max(polygon[:, 0]/1000)])
    ax21m.set_ylim([np.min(polygon[:, 1]/1000), np.max(polygon[:, 1]/1000)])
    ax21m.set_xticks([])
    ax21m.text(np.max(polygon[:, 0]/1000)-1, np.max(polygon[:, 1]/1000)-1, '2016/03/23', ha='right', va='top',
               bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    # ax21m.set_yticks([])
    ax21m.text(487, -2299, 'b)', ha='left', va='bottom', fontsize=12,
               bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))

    ########################################################################
    # 2017 Melange

    ax31m = fig.add_subplot(gs1[2*image_rows+2:3*image_rows+2, :image_cols])
    ax31m.imshow(np.flipud(landsat_image),
                 extent=[np.min(landsat_X)/1000, np.max(landsat_X)/1000, np.min(landsat_Y)/1000, np.max(landsat_Y)/1000],
                 alpha=0.8)
    ax31m.pcolormesh(elev_dict['2017'][0]/1000, elev_dict['2017'][1]/1000, elev_dict['2017'][2], cmap='turbo', vmin=vmin, vmax=vmax)
    ax31m.plot(polygon[:, 0]/1000, polygon[:, 1]/1000, 'y-')
    ax31m.set_xlim([np.min(polygon[:, 0]/1000), np.max(polygon[:, 0]/1000)])
    ax31m.set_ylim([np.min(polygon[:, 1]/1000), np.max(polygon[:, 1]/1000)])
    ax31m.set_xticks([])
    ax31m.text(np.max(polygon[:, 0] / 1000) - 1, np.max(polygon[:, 1] / 1000) - 1, '2017/03/13', ha='right', va='top',
               bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    # ax31m.set_yticks([])
    ax31m.text(487, -2299, 'd)', ha='left', va='bottom', fontsize=12,
               bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))

    ########################################################################
    # 2018 Melange

    ax41m = fig.add_subplot(gs1[3*image_rows+3:4*image_rows+3, :image_cols])
    ax41m.imshow(np.flipud(landsat_image),
                 extent=[np.min(landsat_X)/1000, np.max(landsat_X)/1000, np.min(landsat_Y)/1000, np.max(landsat_Y)/1000],
                 alpha=0.8)
    ax41m.pcolormesh(elev_dict['2018'][0]/1000, elev_dict['2018'][1]/1000, elev_dict['2018'][2], cmap='turbo', vmin=vmin, vmax=vmax)
    ax41m.plot(polygon[:, 0]/1000, polygon[:, 1]/1000, 'y-')
    ax41m.set_xlim([np.min(polygon[:, 0]/1000), np.max(polygon[:, 0]/1000)])
    ax41m.set_ylim([np.min(polygon[:, 1]/1000), np.max(polygon[:, 1]/1000)])
    ax41m.set_xticks([])
    ax41m.text(np.max(polygon[:, 0] / 1000) - 1, np.max(polygon[:, 1] / 1000) - 1, '2018/03/10', ha='right', va='top',
               bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    # ax41m.set_yticks([])
    ax41m.set_ylabel('Distance North (km)')
    ax41m.text(487, -2299, 'f)', ha='left', va='bottom', fontsize=12,
               bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))

    ########################################################################
    # 2019 Melange

    ax51m = fig.add_subplot(gs1[4*image_rows+4:5*image_rows+4, :image_cols])
    ax51m.imshow(np.flipud(landsat_image),
                 extent=[np.min(landsat_X)/1000, np.max(landsat_X)/1000, np.min(landsat_Y)/1000, np.max(landsat_Y)/1000],
                 alpha=0.8)
    ax51m.pcolormesh(elev_dict['2019'][0]/1000, elev_dict['2019'][1]/1000, elev_dict['2019'][2], cmap='turbo', vmin=vmin, vmax=vmax)
    ax51m.plot(polygon[:, 0]/1000, polygon[:, 1]/1000, 'y-')
    ax51m.set_xlim([np.min(polygon[:, 0]/1000), np.max(polygon[:, 0]/1000)])
    ax51m.set_ylim([np.min(polygon[:, 1]/1000), np.max(polygon[:, 1]/1000)])
    ax51m.set_xticks([])
    ax51m.text(np.max(polygon[:, 0] / 1000) - 1, np.max(polygon[:, 1] / 1000) - 1, '2019/03/14', ha='right', va='top',
               bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    # ax51m.set_yticks([])
    ax51m.text(487, -2299, 'h)', ha='left', va='bottom', fontsize=12,
               bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))

    ########################################################################
    # Model Melange

    ax61m = fig.add_subplot(gs1[5*image_rows+5:6*image_rows+5, :image_cols])
    ax61m.imshow(np.flipud(landsat_image),
                 extent=[np.min(landsat_X)/1000, np.max(landsat_X)/1000, np.min(landsat_Y)/1000, np.max(landsat_Y)/1000],
                 alpha=0.8)
    ax61m.pcolormesh(model_X/1000, model_Y/1000, model_elev, cmap='turbo', vmin=vmin, vmax=vmax)
    ax61m.plot(polygon[:, 0]/1000, polygon[:, 1]/1000, 'y-')
    ax61m.set_xlim([np.min(polygon[:, 0]/1000), np.max(polygon[:, 0]/1000)])
    ax61m.set_ylim([np.min(polygon[:, 1]/1000), np.max(polygon[:, 1]/1000)])
    ax61m.set_xlabel('Distance East (km)')
    ax61m.text(np.max(polygon[:, 0] / 1000) - 1, np.max(polygon[:, 1] / 1000) - 1, 'L2 Model', ha='right', va='top',
               bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    # ax61m.set_yticks([])
    ax61m.text(487, -2299, 'j)', ha='left', va='bottom', fontsize=12,
               bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))

    #######################################################################
    # Histograms

    hist_xlim = np.array([np.log10(0.1),np.log10(0.4), np.log10(1),np.log10(4), 1, np.log10(40), 2])
    hist_xlabels = ['0.1', '0.4', '1','4','10','40','100']

    ax12h = fig.add_subplot(gs1[0:int(0.3*image_rows), image_cols + 1:])
    x = np.linspace(vmin,vmax,100)
    y = np.array([0,1])
    X, Y = np.meshgrid(x,y)
    ax12h.pcolormesh(X,Y,X, cmap='turbo')
    # ax12h.set_xticks(hist_xlim)
    # ax12h.set_xticklabels(hist_xlabels)
    ax12h.set_yticks([])
    ax12h.set_xlabel('Melange Freeboard Height (m)')

    ax22h = fig.add_subplot(gs1[1*image_rows+1:2*image_rows+1, image_cols + 1:])
    elev = elev_dict['2016'][2]
    N = plot_bar_chart(ax22h, elev, vmin, vmax)
    ax22h.set_xticks(hist_xlim)
    ax22h.set_xticklabels([])
    ax22h.set_title('Normalized Melange\nSize Distribution')
    ax22h.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax22h.text(np.max(hist_xlim) - 0.005, np.max(N), '2016/03/23', ha='right', va='top')
    ax22h.text(np.min(hist_xlim) + 0.005, np.max(N), 'c)', ha='left', va='top', fontsize=12)

    ax32h = fig.add_subplot(gs1[2*image_rows+2:3*image_rows+2, image_cols + 1:])
    elev = elev_dict['2017'][2]
    N = plot_bar_chart(ax32h, elev, vmin, vmax)
    ax32h.set_xticks(hist_xlim)
    ax32h.set_xticklabels([])
    ax32h.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax32h.text(np.max(hist_xlim) - 0.005, np.max(N), '2017/03/13', ha='right', va='top')
    ax32h.text(np.min(hist_xlim) + 0.005, np.max(N), 'e)', ha='left', va='top', fontsize=12)

    ax42h = fig.add_subplot(gs1[3*image_rows+3:4*image_rows+3, image_cols + 1:])
    elev = elev_dict['2018'][2]
    N = plot_bar_chart(ax42h, elev, vmin, vmax)
    ax42h.set_xticks(hist_xlim)
    ax42h.set_ylabel('Count')
    ax42h.set_xticklabels([])
    ax42h.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax42h.text(np.max(hist_xlim) - 0.005, np.max(N), '2018/03/10', ha='right', va='top')
    ax42h.text(np.min(hist_xlim) + 0.005, np.max(N), 'g)', ha='left', va='top', fontsize=12)

    ax52h = fig.add_subplot(gs1[4*image_rows+4:5*image_rows+4, image_cols + 1:])
    elev = elev_dict['2019'][2]
    N = plot_bar_chart(ax52h, elev, vmin, vmax)
    ax52h.set_xticks(hist_xlim)
    ax52h.set_xticklabels([])
    ax52h.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax52h.text(np.max(hist_xlim) - 0.005, np.max(N), '2019/03/14', ha='right', va='top')
    ax52h.text(np.min(hist_xlim) + 0.005, np.max(N), 'i)', ha='left', va='top', fontsize=12)

    ax62h = fig.add_subplot(gs1[5*image_rows+5:6*image_rows+5, image_cols + 1:])
    elev = model_elev
    N = plot_bar_chart(ax62h, elev, vmin, vmax)
    ax62h.set_xticks(hist_xlim)
    ax62h.set_xticklabels(hist_xlabels)
    ax62h.set_xlabel('Melange Freeboard Height (m)')
    ax62h.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax62h.text(np.max(hist_xlim) - 0.005, np.max(N), 'L2 Model', ha='right', va='top')
    ax62h.text(np.min(hist_xlim) + 0.005, np.max(N), 'k)', ha='left', va='top', fontsize=12)

    plt.savefig(output_file)
    plt.close()


project_folder = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
             'configurations/downscaled_greenland'

polygon = read_melange_polygon_from_shapefile(project_folder)

elev_dict = read_melange_data_from_nc(project_folder)

model_X, model_Y, model_elev = read_model_melange(config_dir, project_folder)

landsat_image_path = os.path.join(project_folder,'Data','Glacier','Kangerlussuaq Melange Landsat.nc')
landsat_X, landsat_Y, landsat_image = read_landsat_subset_from_nc(landsat_image_path)

plot_melange_with_distributions(project_folder, polygon, elev_dict, model_X, model_Y, model_elev,
                                landsat_X, landsat_Y, landsat_image)




