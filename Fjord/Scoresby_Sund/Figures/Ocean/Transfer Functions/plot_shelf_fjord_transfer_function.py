
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.gridspec import GridSpec
import netCDF4 as nc4


def read_modeled_dv_profiles(project_folder,shelf_location,var_name):

    file_name = os.path.join(project_folder,'Data','Modeling','Downscaled','L3_Scoresby_Sund','L3_Scoresby_Sund_dv_CTD_Profiles.nc')
    ds = nc4.Dataset(file_name)
    model_time = ds.variables['time'][:]
    model_depth = ds.variables['depth'][:]
    if shelf_location == 'outer_shelf':
        shelf_model_theta = ds.groups['shelf'].variables[var_name][:,:]
    if shelf_location == 'fjord_entrance':
        shelf_model_theta = ds.groups['fjord_entrance'].variables[var_name][:, :]
    fjord_model_theta = ds.groups['near_glacier'].variables[var_name][:, :]
    ds.close()

    max_fjord_depth_index = np.sum(fjord_model_theta[:,0]!=0)
    shelf_model_theta = shelf_model_theta[:max_fjord_depth_index,:]
    fjord_model_theta = fjord_model_theta[:max_fjord_depth_index,:]
    model_depth = model_depth[:max_fjord_depth_index]

    max_shelf_depth_index = np.sum(shelf_model_theta[:,0]!=0)
    for i in range(max_shelf_depth_index,np.shape(shelf_model_theta)[0]):
        shelf_model_theta[i,:] = shelf_model_theta[max_shelf_depth_index-1,:]

    return(shelf_model_theta, fjord_model_theta, model_time, model_depth, max_fjord_depth_index, max_shelf_depth_index)


def read_modeled_profiles(project_folder,var_name):

    shelf_model_profiles = []
    fjord_model_profiles = []

    interp_profile = np.zeros((1000, 2))
    interp_profile[:, 0] = np.arange(1, 1001, 1)

    file_name = os.path.join(project_folder,'Data','Modeling','Downscaled','L3_Scoresby_Sund','L3_Scoresby_Sund_OMG_Profiles.nc')

    ds = nc4.Dataset(file_name)
    model_depth = ds.variables['depth'][:]

    grp = ds.groups['sound_entrance']
    times = list(grp.groups.keys())
    for t in range(len(times)):
        profile = np.column_stack([model_depth,grp.groups[times[t]].variables[var_name][:]])
        profile = profile[profile[:,1]!=0,:]
        dense_profile = np.copy(interp_profile)
        set_int = interp1d(profile[:,0],profile[:,1])
        for i in range(np.shape(dense_profile)[0]):
            if dense_profile[i,0]<np.max(profile)-1:
                dense_profile[i,1] = set_int(dense_profile[i,0])
            else:
                dense_profile[i,1]=profile[-1,1]
        shelf_model_profiles.append(dense_profile)

    grp = ds.groups['near_glacier']
    times = list(grp.groups.keys())
    for t in range(len(times)):
        profile = np.column_stack([model_depth, grp.groups[times[t]].variables[var_name][:]])

        profile = profile[profile[:, 1] != 0, :]
        dense_profile = np.copy(interp_profile)
        set_int = interp1d(profile[:, 0], profile[:, 1])
        for i in range(np.shape(dense_profile)[0]):
            if dense_profile[i, 0] < np.max(profile) - 1:
                dense_profile[i, 1] = set_int(dense_profile[i, 0])
            else:
                dense_profile[i, 1] = profile[-1, 1]
        fjord_model_profiles.append(dense_profile)

    ds.close()

    return(model_depth, shelf_model_profiles, fjord_model_profiles)


def collect_omg_ctd_profiles(ctd_dir,shelf_location,var_name):

    near_glacier_list = ['CTD_20210820_131357',
                         'CTD_20200904_183011',
                       #  'CTD_20190819_151838',
                         'CTD_20180827_143507',
                         'CTD_20171016_121232']

    interp_profile = np.zeros((1000,2))
    interp_profile[:,0] = np.arange(1,1001,1)

    fjord_profiles = []
    for ctd_name in near_glacier_list:
        year = ctd_name.split('_')[1][:4]
        ds = nc4.Dataset(os.path.join(ctd_dir,year,ctd_name+'.nc'))
        depth = ds.variables['depth'][:]
        if var_name=='Theta':
            pt = ds.variables['potential_temperature'][:]
            profile = np.column_stack([depth,pt])
        if var_name=='Salt':
            sal = ds.variables['practical_salinity'][:]
            profile = np.column_stack([depth, sal])
        ds.close()
        fjord_profile = np.copy(interp_profile)
        for i in range(np.shape(fjord_profile)[0]):
            index = np.argmin(np.abs(profile[:,0]-fjord_profile[i,0]))
            fjord_profile[i,1] = profile[index,1]
        fjord_profiles.append(fjord_profile)

    if shelf_location =='fjord_entrance':
        fjord_entrance_list = ['CTD_20210825_125131',
                               'CTD_20200904_155014',
                               'CTD_20190819_134458',
                               'CTD_20180825_140502',
                               'CTD_20171015_123129']
    if shelf_location=='outer_shelf':
        shelf_list = ['CTD_20210820_155735',
                       'CTD_20200904_123338',
                       #'CTD_20190825_144746',
                       'CTD_20180825_144746',
                       'CTD_20171015_132717']

    shelf_profiles = []
    for ctd_name in shelf_list:
        year = ctd_name.split('_')[1][:4]
        ds = nc4.Dataset(os.path.join(ctd_dir, year, ctd_name + '.nc'))
        depth = ds.variables['depth'][:]
        if var_name == 'Theta':
            pt = ds.variables['potential_temperature'][:]
            profile = np.column_stack([depth, pt])
        if var_name == 'Salt':
            sal = ds.variables['practical_salinity'][:]
            profile = np.column_stack([depth, sal])
        ds.close()
        shelf_profile = np.copy(interp_profile)
        for i in range(np.shape(shelf_profile)[0]):
            index = np.argmin(np.abs(profile[:, 0] - shelf_profile[i, 0]))
            shelf_profile[i, 1] = profile[index, 1]
        shelf_profiles.append(shelf_profile)

    return(shelf_profiles, fjord_profiles)


def plot_theta_shelf_fjord_transfer_functions(project_folder,
                                        shelf_model_theta, fjord_model_theta, model_time, model_depth,
                                        max_fjord_depth_index, max_shelf_depth_index,
                                        shelf_ctd_profiles, fjord_ctd_profiles,
                                              for_publication=False):

    fig = plt.figure(figsize=(10,8))

    # N = 100
    N = len(model_time)

    min_theta = -2
    max_theta = 2

    min_diff = -2.5
    max_diff = 2.5

    plt.subplot(2,3,1)
    for i in range(N):
        plt.plot(shelf_model_theta[:,i],model_depth,'-',linewidth=0.8,color='silver')
    mean_shelf_profile = np.mean(shelf_model_theta,axis=1)
    plt.plot(mean_shelf_profile,model_depth,'-',color='g')
    plt.title('Shelf Profiles')
    plt.plot([min_theta, max_theta],[model_depth[max_shelf_depth_index-1], model_depth[max_shelf_depth_index-1]], 'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.ylabel('High-res Model\nDepth (m)')
    plt.gca().set_xlim([min_theta, max_theta])
    plt.gca().set_ylim([1000,-10])

    plt.subplot(2, 3, 2)
    for i in range(N):
        plt.plot(fjord_model_theta[:, i], model_depth, '-', linewidth=0.8, color='silver')
    mean_fjord_profile = np.mean(fjord_model_theta, axis=1)
    plt.plot(mean_fjord_profile, model_depth, '-', color='g')
    plt.title('Fjord Profiles')
    plt.plot([min_theta, max_theta], [model_depth[max_fjord_depth_index-1], model_depth[max_fjord_depth_index-1]], 'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_xlim([min_theta, max_theta])
    plt.gca().set_ylim([1000, -10])

    plt.subplot(2, 3, 3)
    for i in range(N):
        plt.plot(fjord_model_theta[:, i] - shelf_model_theta[:, i], model_depth, '-', linewidth=0.8, color='silver')
    mean_difference_profile = np.mean(fjord_model_theta-shelf_model_theta, axis=1)
    plt.plot(mean_difference_profile, model_depth, '-', color='g')
    plt.title('Fjord - Shelf Difference')
    plt.plot([min_diff, max_diff], [model_depth[max_shelf_depth_index - 1], model_depth[max_shelf_depth_index - 1]], 'k--')
    plt.plot([min_diff, max_diff], [model_depth[max_fjord_depth_index - 1], model_depth[max_fjord_depth_index - 1]], 'k--')
    plt.grid(linewidth = 0.5, linestyle = '--', alpha = 0.4)
    plt.gca().set_xlim([min_diff, max_diff])
    plt.gca().set_ylim([1000, -10])

    plt.subplot(2, 3, 4)
    mean_shelf_profile = np.zeros_like(shelf_ctd_profiles[0])
    for i in range(len(shelf_ctd_profiles)):
        plt.plot(shelf_ctd_profiles[i][:, 1], shelf_ctd_profiles[i][:, 0], '-', linewidth=0.8, color='silver')
        mean_shelf_profile[:,1] += shelf_ctd_profiles[i][:,1]
        mean_shelf_profile[:,0] = shelf_ctd_profiles[i][:,0]
    mean_shelf_profile[:,1] *= 1/(len(shelf_ctd_profiles))
    plt.plot(mean_shelf_profile[:,1], mean_shelf_profile[:,0], '-', color='g')
    plt.title('Shelf Profiles')
    plt.plot([min_theta, max_theta], [model_depth[max_shelf_depth_index - 1], model_depth[max_shelf_depth_index - 1]], 'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_xlim([min_theta, max_theta])
    plt.gca().set_ylim([1000, -10])
    plt.xlabel('Potential Temperature ($^{\circ}$C)')
    plt.ylabel('CTDs\nDepth (m)')

    plt.subplot(2, 3, 5)
    mean_fjord_profile = np.zeros_like(fjord_ctd_profiles[0])
    for i in range(len(fjord_ctd_profiles)):
        plt.plot(fjord_ctd_profiles[i][:, 1], fjord_ctd_profiles[i][:, 0], '-', linewidth=0.8, color='silver')
        mean_fjord_profile[:, 1] += fjord_ctd_profiles[i][:, 1]
        mean_fjord_profile[:, 0] = fjord_ctd_profiles[i][:, 0]
    mean_fjord_profile[:, 1] *= 1 / (len(fjord_ctd_profiles))
    plt.plot(mean_fjord_profile[:, 1], mean_fjord_profile[:, 0], '-', color='g')
    plt.title('Fjord Profiles')
    plt.plot([min_theta, max_theta], [model_depth[max_fjord_depth_index - 1], model_depth[max_fjord_depth_index - 1]],
             'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_xlim([min_theta, max_theta])
    plt.gca().set_ylim([1000, -10])
    plt.xlabel('Potential Temperature ($^{\circ}$C)')

    plt.subplot(2, 3, 6)
    mean_difference_profile = np.zeros_like(fjord_ctd_profiles[0])
    for i in range(len(fjord_ctd_profiles)):
        plt.plot(fjord_ctd_profiles[i][:, 1] - shelf_ctd_profiles[i][:, 1], fjord_ctd_profiles[i][:, 0], '-', linewidth=0.8, color='silver')
        mean_difference_profile[:, 1] += fjord_ctd_profiles[i][:, 1] - shelf_ctd_profiles[i][:, 1]
        mean_difference_profile[:, 0] = fjord_ctd_profiles[i][:, 0]
    mean_difference_profile[:, 1] *= 1 / (len(fjord_ctd_profiles))
    plt.plot(mean_difference_profile[:, 1], mean_difference_profile[:, 0], '-', color='g')
    # mean_difference_profile = np.mean(fjord_model_theta - shelf_model_theta, axis=1)
    # plt.plot(mean_difference_profile, model_depth, '-', color='g')
    plt.title('Fjord - Shelf Difference')
    plt.plot([min_diff, max_diff], [model_depth[max_shelf_depth_index - 1], model_depth[max_shelf_depth_index - 1]],'k--')
    plt.plot([min_diff, max_diff], [model_depth[max_fjord_depth_index - 1], model_depth[max_fjord_depth_index - 1]],'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_ylim([1000, -10])
    plt.gca().set_xlim([min_diff, max_diff])
    plt.xlabel('Potential Temperature ($^{\circ}$C)')

    output_file = os.path.join(project_folder, 'Figures', 'Ocean', 'Shelf-Fjord Theta Transfer Function.png')
    plt.savefig(output_file)
    plt.close(fig)

def plot_theta_and_salt_shelf_fjord_transfer_functions(project_folder,  #model_time, model_depth,
                                        shelf_model_theta_profiles, fjord_model_theta_profiles,
                                        shelf_model_salt_profiles, fjord_model_salt_profiles,
                                        #max_fjord_depth_index, max_shelf_depth_index,
                                        shelf_ctd_theta_profiles, fjord_ctd_theta_profiles,
                                        shelf_ctd_salt_profiles, fjord_ctd_salt_profiles,
                                                       for_publication=False):

    dpi=300
    fig = plt.figure(figsize=(10,12), dpi=dpi)

    gap_rows = 2
    plot_rows = 6
    gs = GridSpec(gap_rows+4*plot_rows+3, 3, left=0.12, right=0.95, bottom = 0.05, top = 0.92)

    # N = 100
    # N = len(model_time)

    min_theta = -2
    max_theta = 2

    min_theta_diff = -2.5
    max_theta_diff = 2.5

    min_salt = 34
    max_salt = 35

    min_salt_diff = -2.5
    max_salt_diff = 2.5

    fontsize=16


    #################################################################################################################
    # Row 1: Model Theta

    # plt.subplot(4,3,1)
    ax11 = fig.add_subplot(gs[:plot_rows, 0])
    mean_shelf_profile = np.zeros_like(shelf_model_theta_profiles[0])
    for i in range(len(shelf_model_theta_profiles)):
        if i==0:
            plt.plot(shelf_model_theta_profiles[i][:, 1], shelf_model_theta_profiles[i][:, 0],
                     '-', linewidth=0.8,color='silver', label='Individual\nProfile')
        else:
            plt.plot(shelf_model_theta_profiles[i][:, 1], shelf_model_theta_profiles[i][:, 0],
                     '-', linewidth=0.8, color='silver')
        mean_shelf_profile[:, 1] += shelf_model_theta_profiles[i][:, 1]
        mean_shelf_profile[:, 0] = shelf_model_theta_profiles[i][:, 0]
    mean_shelf_profile[:, 1] *= 1 / (len(shelf_model_theta_profiles))
    plt.plot(mean_shelf_profile[:, 1], mean_shelf_profile[:, 0], '--', color='r',alpha=0.6,label='Interpolated\nDownward')
    shallow_indices = mean_shelf_profile[:,0]<400
    plt.plot(mean_shelf_profile[shallow_indices, 1], mean_shelf_profile[shallow_indices, 0], '-', color='r', label='Mean')
    plt.title('Shelf Profiles',fontsize=fontsize)
    # plt.plot([min_theta, max_theta], [model_depth[max_shelf_depth_index - 1], model_depth[max_shelf_depth_index - 1]], 'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_xlim([min_theta, max_theta])
    plt.gca().set_ylim([1000, -10])
    plt.ylabel('L2 Model\n Depth (m)',fontsize=fontsize)
    plt.gca().set_xticklabels([])
    plt.gca().tick_params(axis='both', which='major', labelsize=fontsize)
    # plt.legend(loc=3,fontsize=fontsize)
    plt.text(min_theta+(max_theta-min_theta)*0.97,20,'a)',ha='right',va='top',fontsize=14)

    ax12 = fig.add_subplot(gs[:plot_rows, 1])
    mean_fjord_profile = np.zeros_like(fjord_model_theta_profiles[0])
    for i in range(len(fjord_model_theta_profiles)):
        plt.plot(fjord_model_theta_profiles[i][:, 1], fjord_model_theta_profiles[i][:, 0], '-', linewidth=0.8,
                 color='silver')
        mean_fjord_profile[:, 1] += fjord_model_theta_profiles[i][:, 1]
        mean_fjord_profile[:, 0] = fjord_model_theta_profiles[i][:, 0]
    mean_fjord_profile[:, 1] *= 1 / (len(fjord_model_theta_profiles))
    plt.plot(mean_fjord_profile[:, 1], mean_fjord_profile[:, 0], '-', color='r')
    plt.title('Fjord Profiles',fontsize=fontsize)
    # plt.plot([min_theta, max_theta], [model_depth[max_fjord_depth_index - 1], model_depth[max_fjord_depth_index - 1]],'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_xlim([min_theta, max_theta])
    plt.gca().set_ylim([1000, -10])
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    plt.text(min_theta + (max_theta - min_theta) * 0.97, 20, 'b)', ha='right', va='top', fontsize=14)


    ax12 = fig.add_subplot(gs[:plot_rows, 2])
    mean_difference_profile = np.zeros_like(fjord_model_theta_profiles[0])
    for i in range(len(fjord_model_theta_profiles)):
        plt.plot(fjord_model_theta_profiles[i][:, 1] - shelf_model_theta_profiles[i][:, 1],
                 fjord_model_theta_profiles[i][:, 0], '-', linewidth=0.8, color='silver')
        mean_difference_profile[:, 1] += fjord_model_theta_profiles[i][:, 1] - shelf_model_theta_profiles[i][:, 1]
        mean_difference_profile[:, 0] = fjord_model_theta_profiles[i][:, 0]
    mean_difference_profile[:, 1] *= 1 / (len(fjord_model_theta_profiles))
    plt.plot(mean_difference_profile[:, 1], mean_difference_profile[:, 0], '--', color='r',alpha=0.6)
    shallow_indices = mean_difference_profile[:, 0] < 400
    plt.plot(mean_difference_profile[shallow_indices, 1], mean_difference_profile[shallow_indices, 0], '-', color='r')
    # mean_difference_profile = np.mean(fjord_model_theta - shelf_model_theta, axis=1)
    # plt.plot(mean_difference_profile, model_depth, '-', color='g')
    plt.title('Fjord - Shelf Difference',fontsize=fontsize)
    # plt.plot([min_diff, max_diff], [model_depth[max_shelf_depth_index - 1], model_depth[max_shelf_depth_index - 1]],'k--')
    # plt.plot([min_diff, max_diff], [model_depth[max_fjord_depth_index - 1], model_depth[max_fjord_depth_index - 1]],'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_ylim([1000, -10])
    plt.gca().set_xlim([min_theta_diff, max_theta_diff])
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    plt.text(min_theta_diff + (max_theta_diff - min_theta_diff) * 0.97, 20, 'c)', ha='right', va='top', fontsize=14)

    #################################################################################################################
    # Row 2: CTD Theta

    ax21 = fig.add_subplot(gs[plot_rows+1:2*plot_rows+1, 0])
    mean_shelf_profile = np.zeros_like(shelf_ctd_theta_profiles[0])
    for i in range(len(shelf_ctd_theta_profiles)):
        if i==0:
            plt.plot(shelf_ctd_theta_profiles[i][:, 1], shelf_ctd_theta_profiles[i][:, 0],
                     '-', linewidth=0.8, color='silver',label='Individual\nProfile')
        else:
            plt.plot(shelf_ctd_theta_profiles[i][:, 1], shelf_ctd_theta_profiles[i][:, 0],
                     '-', linewidth=0.8, color='silver')
        mean_shelf_profile[:,1] += shelf_ctd_theta_profiles[i][:,1]
        mean_shelf_profile[:,0] = shelf_ctd_theta_profiles[i][:,0]
    mean_shelf_profile[:,1] *= 1/(len(shelf_ctd_theta_profiles))
    plt.plot(mean_shelf_profile[:,1], mean_shelf_profile[:,0], '--', color='r',label='Extrap.\nDownward',alpha=0.6)
    shallow_indices = mean_shelf_profile[:, 0] < 400
    plt.plot(mean_shelf_profile[shallow_indices, 1], mean_shelf_profile[shallow_indices, 0], '-', color='r',
             label='Mean')
    # plt.title('Shelf Profiles')
    # plt.plot([min_theta, max_theta], [model_depth[max_shelf_depth_index - 1], model_depth[max_shelf_depth_index - 1]], 'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_xlim([min_theta, max_theta])
    plt.gca().set_ylim([1000, -10])
    plt.xlabel('Pot. Temperature ($^{\circ}$C)',fontsize=fontsize)
    plt.ylabel('CTDs\nDepth (m)',fontsize=fontsize)
    plt.legend(loc=3,fontsize=fontsize-3)
    plt.gca().tick_params(axis='both', which='major', labelsize=fontsize)
    plt.text(min_theta + (max_theta - min_theta) * 0.97, 20, 'd)', ha='right', va='top', fontsize=14)

    ax22 = fig.add_subplot(gs[plot_rows+1:2*plot_rows+1, 1])
    mean_fjord_profile = np.zeros_like(fjord_ctd_theta_profiles[0])
    for i in range(len(fjord_ctd_theta_profiles)):
        plt.plot(fjord_ctd_theta_profiles[i][:, 1], fjord_ctd_theta_profiles[i][:, 0], '-', linewidth=0.8, color='silver')
        mean_fjord_profile[:, 1] += fjord_ctd_theta_profiles[i][:, 1]
        mean_fjord_profile[:, 0] = fjord_ctd_theta_profiles[i][:, 0]
    mean_fjord_profile[:, 1] *= 1 / (len(fjord_ctd_theta_profiles))
    plt.plot(mean_fjord_profile[:, 1], mean_fjord_profile[:, 0], '-', color='r')
    # plt.title('Fjord Profiles')
    # plt.plot([min_theta, max_theta], [model_depth[max_fjord_depth_index - 1], model_depth[max_fjord_depth_index - 1]],'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_xlim([min_theta, max_theta])
    plt.gca().set_ylim([1000, -10])
    plt.gca().tick_params(axis='both', which='major', labelsize=fontsize)
    plt.xlabel('Pot. Temperature ($^{\circ}$C)',fontsize=fontsize)
    plt.gca().set_yticklabels([])
    plt.text(min_theta + (max_theta - min_theta) * 0.97, 20, 'e)', ha='right', va='top', fontsize=14)

    ax23 = fig.add_subplot(gs[plot_rows+1:2*plot_rows+1, 2])
    mean_difference_profile = np.zeros_like(fjord_ctd_theta_profiles[0])
    for i in range(len(fjord_ctd_theta_profiles)):
        plt.plot(fjord_ctd_theta_profiles[i][:, 1] - shelf_ctd_theta_profiles[i][:, 1], fjord_ctd_theta_profiles[i][:, 0], '-', linewidth=0.8, color='silver')
        mean_difference_profile[:, 1] += fjord_ctd_theta_profiles[i][:, 1] - shelf_ctd_theta_profiles[i][:, 1]
        mean_difference_profile[:, 0] = fjord_ctd_theta_profiles[i][:, 0]
    mean_difference_profile[:, 1] *= 1 / (len(fjord_ctd_theta_profiles))
    plt.plot(mean_difference_profile[:, 1], mean_difference_profile[:, 0], '--', color='r')
    shallow_indices = mean_difference_profile[:, 0] < 400
    plt.plot(mean_difference_profile[shallow_indices, 1], mean_difference_profile[shallow_indices, 0], '-', color='r')
    # mean_difference_profile = np.mean(fjord_model_theta - shelf_model_theta, axis=1)
    # plt.plot(mean_difference_profile, model_depth, '-', color='g')
    # plt.title('Fjord - Shelf Difference')
    # plt.plot([min_diff, max_diff], [model_depth[max_shelf_depth_index - 1], model_depth[max_shelf_depth_index - 1]],'k--')
    # plt.plot([min_diff, max_diff], [model_depth[max_fjord_depth_index - 1], model_depth[max_fjord_depth_index - 1]],'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_ylim([1000, -10])
    plt.gca().tick_params(axis='both', which='major', labelsize=fontsize)
    plt.gca().set_xlim([min_theta_diff, max_theta_diff])
    plt.xlabel('Pot. Temperature ($^{\circ}$C)',fontsize=fontsize)
    plt.gca().set_yticklabels([])
    plt.text(min_theta_diff + (max_theta_diff - min_theta_diff) * 0.97, 20, 'f)', ha='right', va='top', fontsize=14)

    #################################################################################################################
    # Row 3: Model Salt

    ax31 = fig.add_subplot(gs[gap_rows+2*plot_rows+2:gap_rows+3*plot_rows+2, 0])
    mean_shelf_profile = np.zeros_like(shelf_model_salt_profiles[0])
    for i in range(len(shelf_model_salt_profiles)):
        plt.plot(shelf_model_salt_profiles[i][:, 1], shelf_model_salt_profiles[i][:, 0], '-', linewidth=0.8,
                 color='silver')
        mean_shelf_profile[:, 1] += shelf_model_salt_profiles[i][:, 1]
        mean_shelf_profile[:, 0] = shelf_model_salt_profiles[i][:, 0]
    mean_shelf_profile[:, 1] *= 1 / (len(shelf_model_salt_profiles))
    plt.plot(mean_shelf_profile[:, 1], mean_shelf_profile[:, 0], '--', color='b')
    shallow_indices = mean_shelf_profile[:, 0] < 400
    plt.plot(mean_shelf_profile[shallow_indices, 1], mean_shelf_profile[shallow_indices, 0], '-', color='b')
    plt.title('Shelf Profiles',fontsize=fontsize)
    # plt.plot([min_salt, max_salt], [model_depth[max_shelf_depth_index - 1], model_depth[max_shelf_depth_index - 1]], 'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_xlim([min_salt, max_salt])
    plt.gca().set_xticks([34, 34.25, 34.5, 34.75, 35])
    plt.gca().set_ylim([1000, -10])
    plt.ylabel('L2 Model\nDepth (m)',fontsize=fontsize)
    plt.gca().set_xticklabels([])
    plt.gca().tick_params(axis='both', which='major', labelsize=fontsize)
    plt.text(min_salt + (max_salt - min_salt) * 0.97, 20, 'g)', ha='right', va='top', fontsize=14)

    ax32 = fig.add_subplot(gs[gap_rows+2*plot_rows+2:gap_rows+3*plot_rows+2, 1])
    mean_fjord_profile = np.zeros_like(fjord_model_salt_profiles[0])
    for i in range(len(fjord_model_salt_profiles)):
        plt.plot(fjord_model_salt_profiles[i][:, 1], fjord_model_salt_profiles[i][:, 0], '-', linewidth=0.8,
                 color='silver')
        mean_fjord_profile[:, 1] += fjord_model_salt_profiles[i][:, 1]
        mean_fjord_profile[:, 0] = fjord_model_salt_profiles[i][:, 0]
    mean_fjord_profile[:, 1] *= 1 / (len(fjord_model_salt_profiles))
    plt.plot(mean_fjord_profile[:, 1], mean_fjord_profile[:, 0], '-', color='b')
    plt.title('Fjord Profiles',fontsize=fontsize)
    # plt.plot([min_salt, max_salt], [model_depth[max_fjord_depth_index - 1], model_depth[max_fjord_depth_index - 1]],'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_xlim([min_salt, max_salt])
    plt.gca().set_xticks([34, 34.25, 34.5, 34.75, 35])
    plt.gca().set_ylim([1000, -10])
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    plt.text(min_salt + (max_salt - min_salt) * 0.97, 20, 'h)', ha='right', va='top', fontsize=14)

    ax33 = fig.add_subplot(gs[gap_rows+2*plot_rows+2:gap_rows+3*plot_rows+2, 2])
    mean_difference_profile = np.zeros_like(fjord_model_salt_profiles[0])
    for i in range(len(fjord_model_salt_profiles)):
        plt.plot(fjord_model_salt_profiles[i][:, 1] - shelf_model_salt_profiles[i][:, 1],
                 fjord_model_salt_profiles[i][:, 0], '-', linewidth=0.8, color='silver')
        mean_difference_profile[:, 1] += fjord_model_salt_profiles[i][:, 1] - shelf_model_salt_profiles[i][:, 1]
        mean_difference_profile[:, 0] = fjord_model_salt_profiles[i][:, 0]
    mean_difference_profile[:, 1] *= 1 / (len(fjord_model_salt_profiles))
    plt.plot(mean_difference_profile[:, 1], mean_difference_profile[:, 0], '--', color='b')
    shallow_indices = mean_difference_profile[:, 0] < 400
    plt.plot(mean_difference_profile[shallow_indices, 1], mean_difference_profile[shallow_indices, 0], '-', color='b')
    # mean_difference_profile = np.mean(fjord_model_salt - shelf_model_salt, axis=1)
    # plt.plot(mean_difference_profile, model_depth, '-', color='g')
    plt.title('Fjord - Shelf Difference',fontsize=fontsize)
    # plt.plot([min_diff, max_diff], [model_depth[max_shelf_depth_index - 1], model_depth[max_shelf_depth_index - 1]],'k--')
    # plt.plot([min_diff, max_diff], [model_depth[max_fjord_depth_index - 1], model_depth[max_fjord_depth_index - 1]],'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_ylim([1000, -10])
    plt.gca().set_xlim([min_salt_diff, max_salt_diff])
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    plt.text(min_salt_diff + (max_salt_diff - min_salt_diff) * 0.97, 20, 'i)', ha='right', va='top', fontsize=14)

    #################################################################################################################
    # Row 2: CTD Salt

    ax41 = fig.add_subplot(gs[gap_rows+3*plot_rows+3:gap_rows+4*plot_rows+3, 0])
    mean_shelf_profile = np.zeros_like(shelf_ctd_salt_profiles[0])
    for i in range(len(shelf_ctd_salt_profiles)):

        short_cast = shelf_ctd_salt_profiles[i]
        short_cast = short_cast[short_cast[:,0]>300,:]
        if np.any(short_cast[:,1]<34.4):
            print(short_cast[np.logical_and(short_cast[:,0]>312, short_cast[:,0]<319),:])
            print(i)

        if i==0:
            shelf_ctd_salt_profiles[i][shelf_ctd_salt_profiles[i][:, 0]==316, 1] = 34.8
            plt.plot(shelf_ctd_salt_profiles[i][:, 1], shelf_ctd_salt_profiles[i][:, 0],
                     '-', linewidth=0.8, color='silver',label='Individual\nProfile')
        else:
            plt.plot(shelf_ctd_salt_profiles[i][:, 1], shelf_ctd_salt_profiles[i][:, 0],
                     '-', linewidth=0.8, color='silver')
        mean_shelf_profile[:, 1] += shelf_ctd_salt_profiles[i][:, 1]
        mean_shelf_profile[:, 0] = shelf_ctd_salt_profiles[i][:, 0]
    mean_shelf_profile[:, 1] *= 1 / (len(shelf_ctd_salt_profiles))
    plt.plot(mean_shelf_profile[:, 1], mean_shelf_profile[:, 0], '--', color='b',label='Extrap.\nDownward')
    shallow_indices = mean_shelf_profile[:, 0] < 400
    plt.plot(mean_shelf_profile[shallow_indices, 1], mean_shelf_profile[shallow_indices, 0], '-', color='b',label='Mean')
    # plt.title('Shelf Profiles')
    # plt.plot([min_saltt, max_saltt], [model_depth[max_shelf_depth_index - 1], model_depth[max_shelf_depth_index - 1]],'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_xlim([min_salt, max_salt])
    plt.gca().set_xticks([34, 34.25, 34.5, 34.75, 35])
    plt.gca().set_xticklabels(['34', '', '34.5', '', '35'])
    plt.gca().set_ylim([1000, -10])
    plt.gca().tick_params(axis='both', which='major', labelsize=fontsize)
    plt.xlabel('Salinity (psu)',fontsize=fontsize)
    plt.ylabel('CTDs\nDepth (m)',fontsize=fontsize)
    plt.legend(loc=3,fontsize=fontsize-3)
    plt.text(min_salt + (max_salt - min_salt) * 0.97, 20, 'j)', ha='right', va='top', fontsize=14)

    ax42 = fig.add_subplot(gs[gap_rows+3*plot_rows+3:gap_rows+4*plot_rows+3, 1])
    mean_fjord_profile = np.zeros_like(fjord_ctd_salt_profiles[0])
    for i in range(len(fjord_ctd_salt_profiles)):
        plt.plot(fjord_ctd_salt_profiles[i][:, 1], fjord_ctd_salt_profiles[i][:, 0], '-', linewidth=0.8, color='silver')
        mean_fjord_profile[:, 1] += fjord_ctd_salt_profiles[i][:, 1]
        mean_fjord_profile[:, 0] = fjord_ctd_salt_profiles[i][:, 0]
    mean_fjord_profile[:, 1] *= 1 / (len(fjord_ctd_salt_profiles))
    plt.plot(mean_fjord_profile[:, 1], mean_fjord_profile[:, 0], '-', color='b')
    # plt.title('Fjord Profiles')
    # plt.plot([min_salt, max_salt], [model_depth[max_fjord_depth_index - 1], model_depth[max_fjord_depth_index - 1]], 'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_xlim([min_salt, max_salt])
    plt.gca().set_xticks([34,34.25,34.5,34.75,35])
    plt.gca().set_xticklabels(['34','','34.5','','35'])

    plt.gca().set_ylim([1000, -10])
    plt.gca().tick_params(axis='both', which='major', labelsize=fontsize)
    plt.xlabel('Salinity (psu)',fontsize=fontsize)
    plt.gca().set_yticklabels([])
    plt.text(min_salt + (max_salt - min_salt) * 0.97, 20, 'k)', ha='right', va='top', fontsize=14)

    ax43 = fig.add_subplot(gs[gap_rows+3*plot_rows+3:gap_rows+4*plot_rows+3, 2])
    mean_difference_profile = np.zeros_like(fjord_ctd_salt_profiles[0])
    for i in range(len(fjord_ctd_salt_profiles)):
        plt.plot(fjord_ctd_salt_profiles[i][:, 1] - shelf_ctd_salt_profiles[i][:, 1], fjord_ctd_salt_profiles[i][:, 0], '-',
                 linewidth=0.8, color='silver')
        mean_difference_profile[:, 1] += fjord_ctd_salt_profiles[i][:, 1] - shelf_ctd_salt_profiles[i][:, 1]
        mean_difference_profile[:, 0] = fjord_ctd_salt_profiles[i][:, 0]
    mean_difference_profile[:, 1] *= 1 / (len(fjord_ctd_salt_profiles))
    plt.plot(mean_difference_profile[:, 1], mean_difference_profile[:, 0], '--', color='b')
    shallow_indices = mean_difference_profile[:, 0] < 400
    plt.plot(mean_difference_profile[shallow_indices, 1], mean_difference_profile[shallow_indices, 0], '-', color='b')
    # mean_difference_profile = np.mean(fjord_model_salt - shelf_model_salt, axis=1)
    # plt.plot(mean_difference_profile, model_depth, '-', color='g')
    # plt.title('Fjord - Shelf Difference')
    # plt.plot([min_diff, max_diff], [model_depth[max_shelf_depth_index - 1], model_depth[max_shelf_depth_index - 1]],'k--')
    # plt.plot([min_diff, max_diff], [model_depth[max_fjord_depth_index - 1], model_depth[max_fjord_depth_index - 1]],'k--')
    plt.grid(linewidth=0.5, linestyle='--', alpha=0.4)
    plt.gca().set_ylim([1000, -10])
    plt.gca().tick_params(axis='both', which='major', labelsize=fontsize)
    plt.gca().set_xlim([min_salt_diff, max_salt_diff])
    plt.xlabel('Salinity (psu)',fontsize=fontsize)
    plt.gca().set_yticklabels([])
    plt.text(min_salt_diff + (max_salt_diff - min_salt_diff) * 0.97, 20, 'l)', ha='right', va='top', fontsize=14)

    #################################################################################################################

    plt.suptitle('Temperature and Salinity Modulation through Scoresby Sund',fontsize=fontsize+1)

    if for_publication:
        output_file = os.path.join(project_folder, 'Figures', 'Ocean', 'Transfer Function',
                                   'Shelf-Fjord Theta Transfer Function.pdf')
    else:
        output_file = os.path.join(project_folder, 'Figures', 'Ocean', 'Transfer Function',
                                   'Shelf-Fjord Theta Transfer Function.png')
    plt.savefig(output_file, dpi=dpi)
    plt.close(fig)





project_folder = '/Users/mhwood/Documents/Research/Projects/Scoresby Sund'
ctd_dir = '/Users/mhwood/Documents/Research/Data Repository/Greenland/Ocean Properties/OMG CTDs/Processed'

shelf_location = 'outer_shelf' # or fjord_entrance

# step 1: read in the profiles for the shelf and for the fjord
# shelf_model_theta, fjord_model_theta, model_time, model_depth, max_fjord_depth_index, max_shelf_depth_index = \
#     read_modeled_dv_profiles(project_folder,shelf_location, var_name = 'THETA')


# shelf_model_salt, fjord_model_salt, model_time, model_depth, _, _ = \
#     read_modeled_ctd_profiles(project_folder,shelf_location, var_name = 'SALT')

model_depth, shelf_model_theta_profiles, fjord_model_theta_profiles = \
    read_modeled_profiles(project_folder,var_name = 'Theta')

_, shelf_model_salt_profiles, fjord_model_salt_profiles = \
    read_modeled_profiles(project_folder,var_name = 'Salt')


# step 2: read in the CTD profiles at these locations
shelf_ctd_theta_profiles, fjord_ctd_theta_profiles = collect_omg_ctd_profiles(ctd_dir,shelf_location, var_name='Theta')

shelf_ctd_salt_profiles, fjord_ctd_salt_profiles = collect_omg_ctd_profiles(ctd_dir,shelf_location, var_name='Salt')


#
# # step 3: make a 3-paneled plot with both profiles and the difference
# plot_shelf_fjord_transfer_functions(project_folder,
#                                     shelf_model_theta, fjord_model_theta, model_time, model_depth,
#                                     max_fjord_depth_index, max_shelf_depth_index,
#                                     shelf_ctd_profiles, fjord_ctd_profiles)


plot_theta_and_salt_shelf_fjord_transfer_functions(project_folder,  #model_time, model_depth,
                                        shelf_model_theta_profiles, fjord_model_theta_profiles,
                                        shelf_model_salt_profiles, fjord_model_salt_profiles,
                                        #max_fjord_depth_index, max_shelf_depth_index,
                                        shelf_ctd_theta_profiles, fjord_ctd_theta_profiles,
                                        shelf_ctd_salt_profiles, fjord_ctd_salt_profiles,
                                                   for_publication=True)





