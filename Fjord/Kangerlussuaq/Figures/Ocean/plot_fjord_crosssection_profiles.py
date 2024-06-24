

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon
import cmocean.cm as cm
import shutil

def read_crosssectional_profiles_from_nc(project_dir, var_name, experiment):

    file_path = os.path.join(project_dir,'Data','Ocean','Modeling','L2_Kanger_fjord_crosssection_profiles_'+var_name+'.nc')

    distances = []
    Depths = []
    output_profiles = []

    timestep = 3

    ds = nc4.Dataset(file_path)
    depth = ds.variables['depth'][:]
    for transect in range(1,4):
        distance = ds.groups['transect_'+str(transect)].variables['distance'][:]
        distances.append(distance)
        Depth = ds.groups['transect_' + str(transect)].variables['Depth'][:]
        Depths.append(Depth)
        profile = ds.groups['transect_' + str(transect)].variables[var_name.split('_')[0]+'_'+experiment][:, :, :]
        profile = profile[timestep, :, :]

        for d in range(np.shape(profile)[0]):
            if np.any(profile[d,:]!=0):
                non_zero_indices = np.where(profile[d,:]!=0)[0]
                # print(non_zero_indices)
                profile[d,:non_zero_indices[0]] = profile[d,non_zero_indices[0]]
                profile[d, non_zero_indices[-1]:] = profile[d, non_zero_indices[-1]]

        if var_name in ['Heat_Flux']:
            profile*=1e-3

        output_profiles.append(profile)

    return(depth, distances, Depths, output_profiles)


def plot_crosssectional_profiles(project_dir, var_name, depth, distances, Depths,
                                 profiles_control, profiles_melange, profiles_plume, profiles_melange_plume):

    timestep = 3

    if var_name == 'Theta':
        vmin = -1
        vmax = 2.5
        cmap = 'turbo'

    if var_name == 'Volume_Flux':
        vmin = -0.25
        vmax = 0.25
        cmap = cm.balance
        cbar_label = 'Mean Velocity (m/s)\n$\leftarrow$ Away From Glacier       Toward Glacier $\\rightarrow$ '

    if var_name == 'Heat_Flux':
        vmin = -1e2
        vmax = 1e2
        cmap = cm.balance
        cbar_label = 'Advective Flux of Potential Temperature (10^3 $^{\circ}$Cm$^3$/s)\n$\leftarrow$ Away From Glacier       Toward Glacier $\\rightarrow$ '

    max_depth = 900

    letters = ['a','d','g','j','b','e','h','k','c','f','i','l']

    fig = plt.figure(figsize=(8, 10))

    plot_height = 1
    plot_width = 7

    gs2 = GridSpec(plot_height * 4, plot_width * 3 + 2, left=0.09, right=0.87, bottom=0.08, top=0.95)

    plot_counter = 0
    for transect_number in range(len(distances)):

        for experiment_number in range(4):

            if experiment_number==0:
                profiles = profiles_control
                label = 'Control'
            if experiment_number==1:
                profiles = profiles_melange
                label = 'Melange Only'
            if experiment_number==2:
                profiles = profiles_plume
                label = 'Plume Only'
            if experiment_number==3:
                profiles = profiles_melange_plume
                label = 'Melange and Plume'

            ax = fig.add_subplot(gs2[experiment_number * plot_height:(experiment_number + 1) * plot_height,
                                 transect_number*plot_width:(transect_number+1)*plot_width])
            # ax.pcolormesh(distances[transect_number], depth, profiles[transect_number][:, :],
            #               vmin = vmin, vmax = vmax, cmap = cmap)
            ax.contourf(distances[transect_number], depth, profiles[transect_number][:, :], 50,
                          vmin=vmin, vmax=vmax, cmap=cmap)

            ax.set_ylim([max_depth,0])
            ax.set_xlim([0, np.max(distances[transect_number])])

            Depth_polygon_top = np.column_stack([distances[transect_number], Depths[transect_number]])
            Depth_polygon_bottom = np.copy(np.flipud(Depth_polygon_top))
            Depth_polygon_bottom[:, 1] = max_depth + 10
            Depth_polygon_outline = np.vstack([Depth_polygon_top, Depth_polygon_bottom])
            Depth_polygon = Polygon(Depth_polygon_outline, facecolor='silver', edgecolor='black')
            ax.add_patch(Depth_polygon)

            ax.text(50, max_depth-20, letters[plot_counter]+')',
                    ha='left',va='bottom',color='k',fontsize=12)
            plot_counter+=1

            if experiment_number==0 and transect_number==0:
                ax.set_title('Transect 1\nFjord Mouth')
            if experiment_number==0 and transect_number==1:
                ax.set_title('Transect 2\nMid-Fjord')
            if experiment_number==0 and transect_number==2:
                ax.set_title('Transect 3\n Near-Glacier')
            if transect_number>0:
                ax.set_yticklabels([])
            if experiment_number==3:
                if transect_number in [0]:
                    ax.set_xlabel('Distance\nEastward')
                else:
                    ax.set_xlabel('Distance\nNorthward')
            if experiment_number==1 and transect_number==0:
                ax.set_ylabel('Depth (m)')
            if transect_number==1:
                ax.text(np.max(distances[transect_number])-150,max_depth-20,label,ha='right',va='bottom',color='k',fontsize=12)

    axc = fig.add_subplot(gs2[plot_height:3 * plot_height, -1])
    x = np.array([0, 1])
    y = np.linspace(vmin, vmax, 100)
    X, Y = np.meshgrid(x, y)
    axc.pcolormesh(X, Y, Y, cmap=cmap)
    axc.yaxis.tick_right()
    axc.set_ylabel(cbar_label)
    axc.yaxis.set_label_position("right")
    axc.set_xticks([])

    output_file = os.path.join(project_dir,'Figures','Ocean','Modeling',
                               'Kangerlussuaq Fjord Cross-sectional '+var_name+' Profiles.png')
    plt.savefig(output_file)
    plt.close(fig)


def plot_profiles():

    project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Kangerlussuaq'

    var_name = 'Heat_Flux'

    depth, distances, Depths, profiles_control = read_crosssectional_profiles_from_nc(project_dir, var_name, experiment='control')
    _, _ , _, profiles_melange = read_crosssectional_profiles_from_nc(project_dir, var_name, experiment='melange')
    _, _ , _, profiles_plume = read_crosssectional_profiles_from_nc(project_dir, var_name, experiment='plume')
    _, _ , _, profiles_melange_plume = read_crosssectional_profiles_from_nc(project_dir, var_name, experiment='melange_plume')

    print('control bounds: '+str(np.min(profiles_control[0]))+', '+str(np.max(profiles_control[0])))
    print('melange bounds: '+str(np.min(profiles_melange[0]))+', '+str(np.max(profiles_melange[0])))
    print('plume bounds: '+str(np.min(profiles_plume[0]))+', '+str(np.max(profiles_plume[0])))
    print('melange + plume bounds: '+str(np.min(profiles_melange_plume[0]))+', '+str(np.max(profiles_melange_plume[0])))

    plot_crosssectional_profiles(project_dir, var_name, depth, distances, Depths,
                                 profiles_control, profiles_melange, profiles_plume, profiles_melange_plume)

    shutil.copyfile(os.path.join(project_dir,'Figures','Ocean','Modeling',
                               'Kangerlussuaq Fjord Cross-sectional '+var_name+' Profiles.png'),
                    os.path.join(project_dir, 'Manuscript', 'Draft 1', 'Figure_5.png'))

if __name__ == '__main__':
    plot_profiles()