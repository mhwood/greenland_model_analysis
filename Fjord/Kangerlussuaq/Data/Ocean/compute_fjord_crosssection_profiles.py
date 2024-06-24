
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4


def read_grid_from_nc(config_dir):

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids','L2_Kanger_grid.nc'))
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    dxG = ds.variables['dxG'][:, :]
    dyG = ds.variables['dyG'][:, :]
    drF = ds.variables['drF'][:]
    ds.close()

    Z_bottom = np.cumsum(drF)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2

    return(XC, YC, Depth, Z, dxG, dyG, drF)


def define_cross_fjord_transect(Depth):

    mask = np.copy(Depth)
    mask[mask>1]=1

    transect1_x = np.arange(233,250)
    transect1_y = np.ones_like(transect1_x)*482
    transect1 = np.column_stack([transect1_x, transect1_y])

    transect2_y = np.arange(506, 526)
    transect2_x = np.ones_like(transect2_y) * 207
    transect2 = np.column_stack([transect2_x, transect2_y])

    transect3_y = np.arange(536, 552)
    transect3_x = np.ones_like(transect3_y) * 171
    transect3 = np.column_stack([transect3_x, transect3_y])

    # plt.pcolormesh(mask)
    # plt.plot(transect1[:,0], transect1[:,1],'g-')
    # plt.plot(transect2[:, 0], transect2[:, 1], 'g-')
    # plt.plot(transect3[:, 0], transect3[:, 1], 'g-')
    # plt.gca().set_xlim([150,300])
    # plt.gca().set_ylim([400,600])
    # plt.show()

    transects = [transect1, transect2, transect3]
    return(transects)


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


def sample_geometry_on_transect(transects, XC, YC, Depth, dxG, dyG):

    transect_XCs = []
    transect_YCs = []
    transect_Depths = []
    transect_distances = []

    transect_dxGs = []
    transect_dyGs = []

    for transect in transects:
        if transect[0, 0] == transect[-1, 0]:  # if is a vertical profile
            xc = XC[transect[:, 1], transect[0, 0]]
            yc = YC[transect[:, 1], transect[0, 0]]
            dxg = dxG[transect[:, 1], transect[0, 0]]
            dyg = dyG[transect[:, 1], transect[0, 0]]
            dep = Depth[transect[:, 1], transect[0, 0]]
        else:
            xc = XC[transect[0, 1], transect[:, 0]]
            yc = YC[transect[0, 1], transect[:, 0]]
            dxg = dxG[transect[0, 1], transect[:, 0]]
            dyg = dyG[transect[0, 1], transect[:, 0]]
            dep = Depth[transect[0, 1], transect[:, 0]]

        distance = np.zeros((len(xc),))
        for d in range(len(distance)-1):
            distance[d+1] = distance[d] + great_circle_distance(xc[d], yc[d], xc[d+1], yc[d+1])

        transect_XCs.append(xc)
        transect_YCs.append(yc)
        transect_dxGs.append(dxg)
        transect_dyGs.append(dyg)
        transect_Depths.append(dep)
        transect_distances.append(distance)


    return(transect_XCs, transect_YCs, transect_dxGs, transect_dyGs, transect_Depths, transect_distances)


def sample_fields_on_transect(config_dir, transects, results_dir, subset, years, u_variable_name, v_variable_name):

    output_fields = {}

    # make zero'd transects
    for t in range(len(transects)):
        output_fields[t] = np.zeros((len(years)*12, Nr, np.shape(transects[t])[0]))

    # loop through the transects
    for t in range(len(transects)):

        transect = transects[t]

        if transect[0, 0] == transect[-1, 0]:  # if its a vertical profile
            var_name = u_variable_name
            multiplier = -1
            vertical = True
        else:
            var_name = v_variable_name
            multiplier = 1
            vertical = False

        # loop through the files and fill em in
        for year in years:
            file_name = var_name+'_'+str(year)+'.nc'
            if file_name in os.listdir(os.path.join(config_dir,'L2','L2_Kanger',results_dir,'diags',var_name)):
                print('    - Reading from '+file_name)
                ds = nc4.Dataset(os.path.join(config_dir,'L2','L2_Kanger',results_dir,'diags',var_name,file_name))
                var_grid = ds.variables[var_name][:, :, :, :]
                if vertical:
                    profile = multiplier*var_grid[:, :, transect[:, 1], transect[0, 0]]
                else:
                    profile = multiplier*var_grid[:, :, transect[0, 1], transect[:, 0]]
                output_fields[t][(year-2015)*12:(year+1-2015)*12,:,:] = profile
                ds.close()

    return(output_fields)


def output_profiles_to_nc(project_dir, experiments, output_fields_all, var_name, Z, drF, transect_XCs, transect_YCs, transect_dxGs, transect_dyGs,
                          transect_Depths, transect_distances):

    output_file = os.path.join(project_dir,'Data','Ocean','Modeling','L2_Kanger_fjord_crosssection_profiles_'+var_name+'_Flux.nc')

    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('time',np.shape(output_fields_all[0][0])[0])
    ds.createDimension('depth', np.shape(output_fields_all[0][0])[1])

    v = ds.createVariable('depth','f4',('depth'))
    v[:] = Z

    v = ds.createVariable('drF', 'f4', ('depth'))
    v[:] = drF

    for t in range(len(transect_XCs)):
        grp = ds.createGroup('transect_'+str(t+1))
        grp.createDimension('n_points',len(transect_XCs[t]))

        v = grp.createVariable('XC','f4',('n_points'))
        v[:] = transect_XCs[t]

        v = grp.createVariable('YC', 'f4', ('n_points'))
        v[:] = transect_YCs[t]

        v = grp.createVariable('dxG', 'f4', ('n_points'))
        v[:] = transect_dxGs[t]

        v = grp.createVariable('dyG', 'f4', ('n_points'))
        v[:] = transect_dyGs[t]

        v = grp.createVariable('Depth', 'f4', ('n_points'))
        v[:] = transect_Depths[t]

        v = grp.createVariable('distance', 'f4', ('n_points'))
        v[:] = transect_distances[t]

        for exp in range(len(experiments)):
            experiment = experiments[exp]
            output_fields = output_fields_all[exp]
            v = grp.createVariable(var_name+'_'+experiment, 'f4', ('time','depth','n_points'))
            v[:,:,:] = output_fields[t]

    ds.close()



config_dir = '/Volumes/inngia/Research/downscale_greenland'

project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

XC, YC, Depth, Z, dxG, dyG, drF = read_grid_from_nc(config_dir)

transects = define_cross_fjord_transect(Depth)

transect_XCs, transect_YCs, transect_dxGs, transect_dyGs, transect_Depths, transect_distances = sample_geometry_on_transect(transects, XC, YC, Depth, dxG, dyG)

Nr = 51
years = [2015]

flux_variable = 'Volume' # Volume or Heat

subset = 'budg3d_hflux_set2'

if flux_variable == 'Volume':
    u_variable_name = 'UVELMASS'
    v_variable_name = 'VVELMASS'

if flux_variable == 'Heat':
    u_variable_name = 'ADVx_TH'
    v_variable_name = 'ADVy_TH'

output_fields_all = []

experiments = ['control','melange','melange_plume','plume']
for experiment in experiments:
    print(' - Computing '+flux_variable+' on the cross-sectional transects ('+experiment+')')
    if experiment=='control':
        results_dir = 'results'
    else:
        results_dir = 'results_'+experiment

    print('    - Gathering fields on the transect')
    output_fields = sample_fields_on_transect(config_dir, transects, results_dir, subset, years,
                                              u_variable_name, v_variable_name)

    output_fields_all.append(output_fields)

print('    - Outputting to an nc file')
output_profiles_to_nc(project_dir, experiments, output_fields_all, flux_variable, Z, drF, transect_XCs, transect_YCs, transect_dxGs, transect_dyGs, transect_Depths, transect_distances)



