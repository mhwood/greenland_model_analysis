
import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import datetime

def collect_ctd_profiles(ctd_dir,ctd_list):
    ctd_profiles = []
    for ctd_name in ctd_list:
        year = ctd_name.split('_')[1][:4]
        ds = nc4.Dataset(os.path.join(ctd_dir,year,ctd_name+'.nc'))
        depth = ds.variables['depth'][:]
        pt = ds.variables['potential_temperature'][:]
        profile = np.column_stack([depth,pt])
        ds.close()
        ctd_profiles.append(profile)
    return(ctd_profiles)

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_ecco_profiles_for_ctd_list(project_dir, ecco_region_name, ctd_list):

    ecco_file = os.path.join(project_dir,'Data','Ocean_Properties','ECCO','ECCOv5_THETA_'+ecco_region_name+'_profile_timeseries.nc')
    ds = nc4.Dataset(ecco_file)
    time = ds.variables['time'][:]
    z = ds.variables['RC'][:]
    theta = ds.variables['THETA'][:, :]
    ds.close()

    ecco_profiles = []
    for ctd_id in ctd_list:
        year = int(ctd_id.split('_')[1][:4])
        if year<=2017:
            month = int(ctd_id.split('_')[1][4:6])
            day = int(ctd_id.split('_')[1][6:8])
            dec_yr = YMD_to_DecYr(year,month,day)
            index = np.argmin(np.abs(time-dec_yr))
            ecco_profile = np.column_stack([-z,theta[:,index]])
            ecco_profiles.append(ecco_profile)

    return(ecco_profiles)

def create_all_ctd_plots(project_dir, ctd_dir, plot_titles, ecco_region_names, ctd_lists):

    fig = plt.figure(figsize=(12,8))

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    year_to_color = {2016:'silver',2017:'silver',2018:colors[2],2019:'silver',2020:colors[4],2021:'silver'}

    for m in range(len(plot_titles)):
        ctd_profiles = collect_ctd_profiles(ctd_dir,ctd_lists[m])
        # ecco_region_name = ecco_region_names[m]
        # if ecco_region_name!='':
        #     ecco_profiles = read_ecco_profiles_for_ctd_list(project_dir, ecco_region_name, ctd_lists[m])
        # else:
        #     ecco_profiles = []

        print('    - Plotting the '+plot_titles[m]+' profiles')

        ##################################################################################
        # ctds

        print('        - Plotting '+str(len(ctd_profiles))+' profiles')
        plt.subplot(1,4,m+1)
        # for c in range(len(ecco_profiles)):
        #     profile = ecco_profiles[c]
        #     plt.plot(profile[:, 1], profile[:, 0], alpha=0.5, color='silver')
        for c in range(len(ctd_profiles)):
            profile = ctd_profiles[c]
            ctd_id = ctd_lists[m][c]
            year= ctd_id.split('_')[1][:4]
            plt.plot(profile[:,1],profile[:,0],label=year,color=year_to_color[int(year)])

        plt.gca().set_xlim([-1.9,3])
        # plt.plot([-1.9,3],[200,200],'k-')
        # plt.plot([-1.9, 3], [500, 500], 'k-')
        plt.gca().set_ylim([0,1000])

        plt.grid(linestyle='--',alpha=0.5)

        plt.gca().invert_yaxis()

        if m==1:
            plt.legend()
        plt.title(plot_titles[m])

        # ##################################################################################
        # # ecco
        #
        # plt.subplot(2, 4, m + 5)
        # print('        - Plotting ' + str(len(ecco_profiles)) + ' profiles')
        # for c in range(len(ecco_profiles)):
        #     profile = ecco_profiles[c]
        #     ctd_id = ctd_lists[m][c]
        #     year = ctd_id.split('_')[1][:4]
        #     print(c,ctd_id)
        #     plt.plot(profile[:, 1], profile[:, 0], label=year, color=year_to_color[int(year)])
        #
        # plt.gca().set_xlim([-1.9, 3])
        # plt.plot([-1.9, 3], [200, 200], 'k-')
        # plt.plot([-1.9, 3], [500, 500], 'k-')
        # plt.gca().set_ylim([0, 1000])
        #
        # plt.grid(linestyle='--', alpha=0.5)
        #
        # plt.gca().invert_yaxis()

    output_file = os.path.join(project_dir, 'Figures','Scoresby Sund OMG Fjord Temperature.png')
    plt.savefig(output_file, bbox_inches='tight')
    plt.close(fig)



near_glacier_list = ['CTD_20210820_131357',
                     'CTD_20200904_183011',
                     'CTD_20190819_151838',
                     'CTD_20180827_143507',
                     'CTD_20171016_121232']

mid_fjord_list = ['CTD_20210820_134800',
                  'CTD_20200904_175044',
                  'CTD_20190819_143308',
                  'CTD_20180827_152500',
                  'CTD_20171016_130337',
                  'CTD_20161010_140634']

fjord_entrance_list = ['CTD_20210825_125131',
                       'CTD_20200904_155014',
                       'CTD_20190819_134458',
                       'CTD_20180825_140502',
                       'CTD_20171015_123129',
                       'CTD_20161010_145001']

upstream_list = ['CTD_20171016_110925',
                 'CTD_20160929_150010']#,'CTD_20160929_145648']

project_dir = '/Users/michwood/Documents/Research/Projects/Scoresby Sund'

ctd_dir = '/Users/michwood/Documents/Research/Data Repository/Greenland/Ocean Properties/OMG CTDs/Processed'
# ctd_profiles = collect_ctd_profiles(ctd_dir,near_glacier_list)

plot_titles = ['Near-Glacier','Mid-Fjord','Fjord Entrance','Shelf Break']
output_names =['Near_Glacier','Mid_Fjord','Fjord_Entrance','Shelf_Break']
ecco_region_names = ['','Scoresby_Sund_Fjord','Shelf_Entrance','Shelf_Break']

locations = [(-28.3818,71.9304),(-24.9112,71.0983),(-21.3842,70.1023),(-17.515,71.976)]
ctd_lists = [near_glacier_list, sorted(mid_fjord_list), sorted(fjord_entrance_list), sorted(upstream_list)]

create_all_ctd_plots(project_dir, ctd_dir, plot_titles, ecco_region_names, ctd_lists)


