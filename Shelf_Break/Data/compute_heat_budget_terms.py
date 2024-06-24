
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4


def read_grid_param(ecco_dir, param_name):

    three_D_fields = ['hFacC','hFacS','hFacW']
    vertical_fields = ['DRF','DRC','RC','RF']

    if param_name in three_D_fields:
        output_grid = np.zeros((50, 13, 270, 270))
    else:
        output_grid = np.zeros((13, 270, 270))

    if param_name in vertical_fields:
        max_tile = 1
    else:
        max_tile = 13

    for tile in range(1,max_tile+1):

        file_path = os.path.join(ecco_dir, 'GRID', 'GRID.'+'{:04d}'.format(tile)+'.nc')
        ds = nc4.Dataset(file_path)
        if param_name in three_D_fields:
            grid = ds.variables[param_name][:, :, :]
        elif param_name in vertical_fields:
            grid = ds.variables[param_name][:]
        else:
            grid = ds.variables[param_name][:, :]
        ds.close()

        if param_name in three_D_fields:
            output_grid[:, tile-1, :, :] = grid
        elif param_name in vertical_fields:
            output_grid = grid
        else:
            output_grid[tile - 1, :, :] = grid

    return(output_grid)


def extend_vertical_dim_to_3d(vertical_array):

    new_array = np.zeros((50, 13, 270, 270))

    for d in range(len(vertical_array)):
        new_array[d, :, :, :] = vertical_array[d]

    return(new_array)


def extend_2d_array_to_3d(two_D_array):

    new_array = np.zeros((50, 13, 270, 270))

    for d in range(50):
        new_array[d, :, :, :] = two_D_array

    return(new_array)


def read_field(ecco_dir, field_name, year):

    file_path = os.path.join(ecco_dir, field_name, field_name+'_'+str(year)+'.nc')
    ds = nc4.Dataset(file_path)
    if 'k' in ds.variables.keys():
        grid = ds.variables[field_name][:, :, :, :, :]
    else:
        grid = ds.variables[field_name][:, :, : ,:]
    ds.close()

    return(grid)


def compute_advection_tendency(CELL_VOL, ADVr_TH, ADVx_TH, ADVy_TH):

    vertical_advection_small = -1*(np.diff(ADVr_TH,axis=1))
    vertical_advection = np.zeros((12,50,13,270,270))
    vertical_advection[:, :-1, :, :, :] = vertical_advection_small
    vertical_advection[:,-1, :, :, :] = vertical_advection_small[:,-1,:,:,:]

    horizontal_advection_x_small = -1 * (np.diff(ADVx_TH, axis=-1))
    horizontal_advection_x = np.zeros((12, 50, 13, 270, 270))
    horizontal_advection_x[:, :, :, :, :-1] = horizontal_advection_x_small
    horizontal_advection_x[:, :, :, :, -1] = horizontal_advection_x_small[:, :, :, :, -1]

    horizontal_advection_y_small = -1 * (np.diff(ADVy_TH, axis=-2))
    horizontal_advection_y = np.zeros((12, 50, 13, 270, 270))
    horizontal_advection_y[:, :, :, :-1, :] = horizontal_advection_y_small
    horizontal_advection_y[:, :, :, -1, :] = horizontal_advection_y_small[:, :, :, -1, :]

    adv = (vertical_advection + horizontal_advection_x + horizontal_advection_y)/CELL_VOL

    return(adv)


def compute_diffusion_tendency(CELL_VOL, DFrI_TH, DFrE_TH, DFxE_TH, DFyE_TH):

    vertical_diffusion_I_small = -1*(np.diff(DFrI_TH,axis=1))
    vertical_diffusion_I = np.zeros((12,50,13,270,270))
    vertical_diffusion_I[:, :-1, :, :, :] = vertical_diffusion_I_small
    vertical_diffusion_I[:,-1, :, :, :] = vertical_diffusion_I_small[:,-1,:,:,:]

    vertical_diffusion_E_small = -1 * (np.diff(DFrE_TH, axis=1))
    vertical_diffusion_E = np.zeros((12, 50, 13, 270, 270))
    vertical_diffusion_E[:, :-1, :, :, :] = vertical_diffusion_E_small
    vertical_diffusion_E[:, -1, :, :, :] = vertical_diffusion_E_small[:, -1, :, :, :]

    horizontal_diffusion_x_small = -1 * (np.diff(DFxE_TH, axis=-1))
    horizontal_diffusion_x = np.zeros((12, 50, 13, 270, 270))
    horizontal_diffusion_x[:, :, :, :, :-1] = horizontal_diffusion_x_small
    horizontal_diffusion_x[:, :, :, :, -1] = horizontal_diffusion_x_small[:, :, :, :, -1]

    horizontal_diffusion_y_small = -1 * (np.diff(DFyE_TH, axis=-2))
    horizontal_diffusion_y = np.zeros((12, 50, 13, 270, 270))
    horizontal_diffusion_y[:, :, :, :-1, :] = horizontal_diffusion_y_small
    horizontal_diffusion_y[:, :, :, -1, :] = horizontal_diffusion_y_small[:, :, :, -1, :]

    dif = (vertical_diffusion_I + vertical_diffusion_E + horizontal_diffusion_x + horizontal_diffusion_y)/CELL_VOL

    return(dif)


ecco_dir = '/Volumes/petermann/Research/Projects/Ocean_Modeling/ECCO'

year = 1992
month = 1
tile = 10

print('        - Reading RAC')
RAC = read_grid_param(ecco_dir, 'RAC')
RAC = extend_2d_array_to_3d(RAC)
print('        - Reading DRF')
DRF = read_grid_param(ecco_dir, 'DRF')
DRF = extend_vertical_dim_to_3d(DRF)
print('        - Reading hFacC')
hFacC = read_grid_param(ecco_dir, 'hFacC')

CELL_VOL = hFacC * DRF * RAC

# print('    - Computing the advective tendency for '+str(year)+'/'+str(month))
#
# print('        - Reading ADVr_TH')
# ADVr_TH = read_field(ecco_dir, 'ADVr_TH', year)
#
# print('        - Reading ADVx_TH')
# ADVx_TH = read_field(ecco_dir, 'ADVx_TH', year)
#
# print('        - Reading ADVy_TH')
# print('        - Reading ADVr_TH')
# ADVr_TH = read_field(ecco_dir, 'ADVr_TH', year)
#
# print('        - Reading ADVx_TH')
# ADVx_TH = read_field(ecco_dir, 'ADVx_TH', year)
#
# print('        - Reading ADVy_TH')
# ADVy_TH = read_field(ecco_dir, 'ADVy_TH', year)
#
# compute_adve
# ADVy_TH = read_field(ecco_dir, 'ADVy_TH', year)
#
# advection = compute_advection_tendency(CELL_VOL, ADVr_TH, ADVx_TH, ADVy_TH)

print('    - Computing the diffusive tendency for '+str(year)+'/'+str(month))

print('        - Reading DFrE_TH')
DFrE_TH = read_field(ecco_dir, 'DFrE_TH', year)

print('        - Reading DFrI_TH')
DFrI_TH = read_field(ecco_dir, 'DFrI_TH', year)

print('        - Reading DFxE_TH')
DFxE_TH = read_field(ecco_dir, 'DFxE_TH', year)

print('        - Reading DFyE_TH')
DFyE_TH = read_field(ecco_dir, 'DFyE_TH', year)

diffusion = compute_diffusion_tendency(CELL_VOL, DFrI_TH, DFrE_TH, DFxE_TH, DFyE_TH)

plt.imshow(diffusion[0,0,0,:,:])
plt.show()



# THETA = read_field(ecco_dir, 'THETA', year)



# THETA = THETA[month-1, :, 10, :, :]
#
# plt.pcolormesh(THETA[0, :, :])
# plt.show()

