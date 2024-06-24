
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import netCDF4 as nc4
import shapefile
import datetime

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)


def read_model_grid(config_dir, model_name):
    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids',model_name+'_grid.nc'))
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    ds.close()
    return(XC, YC)


def read_sample_areas_from_shapefile(project_folder, shapefile_name):

    file_path = os.path.join(project_folder, 'Map', 'Shapefiles', shapefile_name)

    sf = shapefile.Reader(file_path)
    shapes = sf.shapes()
    records = sf.records()

    polygons = []
    polygon_names = []

    for i in range(len(shapes)):
        record = records[i]
        polygon_names.append(record[0])
        polygons.append(np.array(shapes[i].points))

    return(polygons, polygon_names)


def compute_inside_indices(XC, YC, bounding_box):
    p = mplPath.Path(bounding_box)
    inside = p.contains_points(np.column_stack([XC.ravel(), YC.ravel()]))
    inside = np.reshape(inside, np.shape(XC))
    return(inside)


def create_timeseries_from_model_files(config_dir, subset_name, polygon_indices):

    results_dir = os.path.join(config_dir,'L1','L1_W_Greenland','results',subset_name)

    all_timeseries = []
    for t in range(len(polygon_indices)):
        all_timeseries.append(np.zeros((30*366, 4)).astype(float))

    counter = 0

    for file_name in sorted(os.listdir(results_dir)):
        if file_name[0]!='.' and file_name[-3:]=='.nc':

            year = int(file_name.split('.')[1][:4])
            month = int(file_name.split('.')[1][4:6])

            print('    - Reading from '+file_name)

            ds = nc4.Dataset(os.path.join(results_dir,file_name))
            if subset_name=='EtaN_day_snap':
                grid = ds.variables['EtaN'][:, :, :]
            if subset_name == 'BGC_daily_Chl':
                grid = ds.variables['Chl01'][:, :, :]
                grid += ds.variables['Chl02'][:, :, :]
                grid += ds.variables['Chl03'][:, :, :]
                grid += ds.variables['Chl04'][:, :, :]
                grid += ds.variables['Chl05'][:, :, :]
            ds.close()

            for d in range(np.shape(grid)[0]):

                dec_yr = YMD_to_DecYr(year,month,d+1)

                grid_subset = grid[d, :, :]

                for t in range(len(polygon_indices)):
                    indices = np.logical_and(polygon_indices[t], grid_subset != 0)

                    all_timeseries[t][counter + d, 0] = dec_yr
                    all_timeseries[t][counter + d, 1] = np.mean(grid_subset[indices])
                    all_timeseries[t][counter + d, 2] = np.std(grid_subset[indices], ddof=1)
                    all_timeseries[t][counter + d, 3] = np.sum(indices)

            counter += np.shape(grid)[0]

    for t in range(len(all_timeseries)):
        all_timeseries[t] = all_timeseries[t][all_timeseries[t][:, 0] != 0]

    return(all_timeseries)


def save_timeseries_as_nc(project_dir, output_name, all_timeseries, polygon_names):
    ds = nc4.Dataset(os.path.join(project_dir, 'Data','Modeling', 'L1', output_name+'.nc'), 'w')
    ds.createDimension('dec_yr', np.shape(all_timeseries[0])[0])

    t = ds.createVariable('dec_yr', 'f4', ('dec_yr',))
    t[:] = all_timeseries[0][:, 0]

    for t in range(len(all_timeseries)):
        region_name = '_'.join(polygon_names[t].split())
        grp = ds.createGroup(region_name)
        v = grp.createVariable('mean', 'f4', ('dec_yr',))
        v[:] = all_timeseries[t][:, 1]
        v = grp.createVariable('std', 'f4', ('dec_yr',))
        v[:] = all_timeseries[t][:, 2]
        v = grp.createVariable('count', 'f4', ('dec_yr',))
        v[:] = all_timeseries[t][:, 3]

    ds.close()


# config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/' \
#              'MITgcm/configurations/downscale_darwin'
config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/' \
             'configurations/downscale_darwin'

# project_folder = '/Users/mike/Documents/Research/Projects/Disko Bay'
project_folder = '/Users/mhwood/Documents/Research/Projects/Disko Bay'

shapefile_name = 'Arrigo et al 2017 Sample Areas'
subset_name = 'BGC_daily_Chl'
output_name = 'Arrigo_et_al_Chl_Timeseries'

# read in the model grid
XC, YC = read_model_grid(config_dir,model_name='L1_W_Greenland')

# read the polygons from the shapefile
polygons, polygon_names = read_sample_areas_from_shapefile(project_folder, shapefile_name)

# # find indices in the regions
polygon_indices = []
for polygon in polygons:
    inside = compute_inside_indices(XC, YC, polygon)
    polygon_indices.append(inside)

# subset the files to create timeseries
all_timeseries = create_timeseries_from_model_files(config_dir, subset_name, polygon_indices)

save_timeseries_as_nc(project_folder, output_name, all_timeseries, polygon_names)

for t in range(len(all_timeseries)):
    plt.subplot(len(all_timeseries),1,t+1)
    plt.plot(all_timeseries[t][:,0], all_timeseries[t][:,1], 'k.')
    plt.ylabel(polygon_names[t])


plt.show()




