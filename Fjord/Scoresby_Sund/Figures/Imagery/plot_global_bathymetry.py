
import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cmocean.cm as cm

def read_geometry_from_grid_nc(config_dir,level_name,model_name):

    grid_path = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')

    ds = nc4.Dataset(grid_path)
    Depth = ds.variables['Depth'][:, :]
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    ds.close()

    bottom = np.column_stack([XC[0, :], YC[0, :]])
    top = np.column_stack([XC[-1, :], YC[-1, :]])
    left = np.column_stack([XC[:, 0], YC[:, 0]])
    right = np.column_stack([XC[:, -1], YC[:, -1]])
    polygon = np.vstack([bottom, right, np.flipud(top), np.flipud(left)])

    return(XC, YC, Depth, polygon)

def read_ecco_field_to_faces(file_path, llc, dim):

    Nr = 50

    grid_array = np.fromfile(file_path,'>f4')
    N = 13*llc*llc

    field_faces = {}

    if dim == 2:
        points_counted = 0
        for i in range(1, 6):
            if i < 3:
                n_points = 3 * llc * llc
                grid = grid_array[points_counted:points_counted + n_points]
                grid = np.reshape(grid, (3 * llc, llc))
            if i == 3:
                n_points = llc * llc
                grid = grid_array[points_counted:points_counted + n_points]
                grid = np.reshape(grid, (llc, llc))
            if i > 3:
                n_points = 3 * llc * llc
                grid = grid_array[points_counted:points_counted + n_points]
                grid = np.reshape(grid, (llc, 3 * llc))
            field_faces[i] = grid
            points_counted += n_points

    if dim==3:

        for i in range(1, 6):
            if i < 3:
                face_grid = np.zeros((Nr, 3 * llc, llc))
            elif i == 3:
                face_grid = np.zeros((Nr, llc, llc))
            if i > 3:
                face_grid = np.zeros((Nr, llc, 3 * llc))
            field_faces[i]=face_grid

        for nr in range(Nr):
            points_counted = 0
            level_grid = grid_array[nr * N:(nr + 1) * N]
            for i in range(1,6):
                if i < 3:
                    n_points = 3*llc*llc
                    grid = level_grid[points_counted:points_counted+n_points]
                    grid = np.reshape(grid,(3*llc,llc))
                if i == 3:
                    n_points = llc * llc
                    grid = level_grid[points_counted:points_counted + n_points]
                    grid = np.reshape(grid, (llc, llc))
                if i > 3:
                    n_points = 3 * llc * llc
                    grid = level_grid[points_counted:points_counted + n_points]
                    grid = np.reshape(grid, (llc, 3*llc))
                field_faces[i][nr,:,:] = grid

                points_counted += n_points

        # plt.imshow(grid,origin='lower')
        # plt.show()

    return(field_faces)

def read_L0_grid(ecco_dir):
    llc = 270

    bathy_file = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','input_init','bathy_llc'+str(llc))
    bathy_faces = read_ecco_field_to_faces(bathy_file, llc, dim=2)

    grid_file_dir = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles')
    XC_faces = {}
    YC_faces = {}
    for i in [1,2,3,4,5]:
        if i < 3:
            grid = np.fromfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), '>f8')
            grid = np.reshape(grid, (16, 3*llc + 1, llc + 1))
            XC_faces[i] = grid[0, :-1, :-1]
            YC_faces[i] = grid[1, :-1, :-1]
        if i == 3:
            grid = np.fromfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'),'>f8')
            grid = np.reshape(grid, (16, llc+1, llc+1))
            XC_faces[i] = grid[0,:-1,:-1]
            YC_faces[i] = grid[1,:-1,:-1]
        if i > 3:
            grid = np.fromfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), '>f8')
            grid = np.reshape(grid, (16, llc + 1, 3*llc + 1))
            XC_faces[i] = grid[0, :-1, :-1]
            YC_faces[i] = grid[1, :-1, :-1]


    return (XC_faces, YC_faces, bathy_faces)

def generate_global_plot(output_path, XC_faces, YC_faces, bathy_faces,
                  center_lon, rotation_lon, center_lat, rotation_lat,
                  polygon):

    plot_mode = 'dark'

    with_land = True

    fig = plt.figure(figsize=(12,12))
    if plot_mode=='dark':
        plt.style.use("dark_background")

    m = Basemap(projection='ortho', resolution=None,
                lat_0=center_lat + rotation_lat, lon_0=center_lon + rotation_lon)
    if with_land:
        m.bluemarble(scale=0.5)


    for face in [1,2,3,4,5]:
        xc_grid = XC_faces[face]
        yc_grid = YC_faces[face]
        var_grid = bathy_faces[face]
        var_grid[var_grid>0] = 0

        if face==2:
            xc_grid-=360

        if with_land:
            var_grid = np.ma.masked_where(var_grid==0, var_grid)

        if 'XY_270_corrected_face_'+str(face)+'.nc' in os.listdir(os.path.join(config_dir, 'plots','Basemap Plot Reference')):
            ds = nc4.Dataset(os.path.join(config_dir, 'plots', 'Basemap Plot Reference', 'XY_270_corrected_face_' + str(face) + '.nc'))
            X = ds.variables['X'][:,:]
            Y = ds.variables['Y'][:,:]
            ds.close()
        else:
            X, Y = m(xc_grid, yc_grid)
            if np.any(np.abs(X) > 1e10):
                rows = np.arange(np.shape(X)[0])
                cols = np.arange(np.shape(Y)[1])
                Cols, Rows = np.meshgrid(cols, rows)
                Cols = Cols.ravel()
                Rows = Rows.ravel()
                X_ravel = X.ravel()
                Y_ravel = Y.ravel()
                Cols = Cols[np.abs(X_ravel) < 1e10]
                Rows = Rows[np.abs(X_ravel) < 1e10]
                X_ravel = X_ravel[np.abs(X_ravel) < 1e10]
                Y_ravel = Y_ravel[np.abs(Y_ravel) < 1e10]
                nan_rows, nan_cols = np.where(np.abs(X) > 1e10)
                for i in range(len(nan_rows)):
                    if i%100==0:
                        print(i,len(nan_rows),round(100*i/len(nan_rows),2))
                    row = nan_rows[i]
                    col = nan_cols[i]
                    closest_index = np.argmin((Cols - col) ** 2 + (Rows - row) ** 2)
                    X[row, col] = X_ravel[closest_index]
                    Y[row, col] = Y_ravel[closest_index]
            ds = nc4.Dataset(os.path.join(config_dir, 'plots','Basemap Plot Reference','XY_270_corrected_face_'+str(face)+'.nc'),'w')
            ds.createDimension('rows',np.shape(X)[0])
            ds.createDimension('cols', np.shape(X)[1])
            Xvar = ds.createVariable('X','f4',('rows','cols'))
            Yvar = ds.createVariable('Y', 'f4', ('rows', 'cols'))
            Xvar[:] = X
            Yvar[:] = Y
            ds.close()
        # m.pcolormesh(X, Y, var_grid, vmin=-5000, vmax = 5000, cmap = cm.topo)
        print(np.shape(X))
        if face==2:
            m.pcolormesh(X[400:,150:], Y[400:,150:], var_grid[400:,150:], vmin=-1500, vmax=0, cmap=cm.deep_r)
        else:
            m.pcolormesh(X, Y, var_grid, vmin=-1500, vmax=0, cmap=cm.deep_r)

    polygon_lon, polygon_lat = m(polygon[:, 0], polygon[:, 1])
    m.plot(polygon_lon, polygon_lat, 'r-', linewidth=3)

    axicon = fig.add_axes([0.22, 0.9, 0.5, 0.1])
    # plt.text(0, 0, 'Global Source Model (1/3$^{\circ}$)', fontsize=20)
    axicon.axis('off')
    axicon.set_xticks([])
    axicon.set_yticks([])

    plt.savefig(output_path)
    plt.close(fig)

def create_globe_domain_plot(output_dir, ecco_dir, config_dir, level_name, model_name):

    XC, YC, Depth, polygon = read_geometry_from_grid_nc(config_dir, level_name, model_name)

    XC_faces, YC_faces, bathy_faces = read_L0_grid(ecco_dir)

    # plt.subplot(1,3,1)
    # plt.imshow(XC_faces[2])
    #
    # plt.subplot(1, 3, 2)
    # plt.imshow(YC_faces[2])
    #
    # plt.subplot(1, 3, 3)
    # plt.imshow(bathy_faces[2],vmin=-1,vmax=0)
    #
    # plt.show()

    ##################################
    # some metadata

    center_lon = np.mean(XC)
    center_lat = np.mean(YC)

    rotation_lon = -45
    rotation_lat = -20

    output_path = os.path.join(output_dir,'L0_global_bathymetry.png')

    generate_global_plot(output_path, XC_faces, YC_faces, bathy_faces,
                  center_lon, rotation_lon, center_lat, rotation_lat,
                  polygon)

# config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/configurations/' \
#              'downscale_darwin'
config_dir = '/Volumes/petermann/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
             'configurations/downscale_greenland'

# output_dir = '/Users/mhwood/Documents/Research/Projects/Disko Bay/Figures'
output_dir = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund/Figures'

ecco_dir = '/Users/mhwood/Documents/Research/Projects/Ocean_Modelling/ECCO'

level_name = 'L1'
model_name = 'L1_W_Greenland'
model_name = 'L1_CE_Greenland'

create_globe_domain_plot(output_dir, ecco_dir, config_dir, level_name, model_name)