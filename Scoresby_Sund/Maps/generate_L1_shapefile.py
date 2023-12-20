

import netCDF4 as nc4
import shapefile
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Transformer

def reproject_points(points,inputCRS,outputCRS,x_column=0,y_column=1):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(points[:, x_column], points[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
        run_test = False
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon = np.copy(points)
    output_polygon[:, x_column] = x2
    output_polygon[:, y_column] = y2
    return output_polygon


grid_file = '/Volumes/mhwood/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/configurations/' \
            'downscaled_greenland/nc_grids/L1_CE_Greenland_grid.nc'

ds = nc4.Dataset(grid_file)
xc = ds.variables['XC'][:,:]
yc = ds.variables['YC'][:,:]
ds.close()

north = np.column_stack([xc[-1,:],yc[-1,:]])
south = np.column_stack([xc[0,:],yc[0,:]])
west = np.column_stack([xc[:,0],yc[:,0]])
east = np.column_stack([xc[:,-1],yc[:,-1]])

boundary = np.vstack([south,east,np.flipud(north),np.flipud(west)])

# boundary = reproject_points(boundary,4326,3413)

output_file = '/Users/michwood/Documents/Research/Projects/Scoresby Sund/Map/Shapefiles/Domains/L1_CE_Greenland_boundary'

sf = shapefile.Writer(output_file)

sf.field('Domain','C')
sf.poly([boundary.tolist()])
sf.record('L1_CE_Greenland')
sf.close()

# plt.plot(boundary[:,0],boundary[:,1])
# plt.show()