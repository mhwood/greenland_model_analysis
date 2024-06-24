
import os
import netCDF4 as nc4
import matplotlib.pyplot as plt
import numpy as np
from datetime import timedelta, datetime
# from ecco_v4_py import llc_array_conversion as ec
from MITgcmutils import mds
import llc_array_conversion as ec





def file_iter_to_date(file_iter):
    total_seconds = file_iter*1200-(15*24*60*60)
    date = datetime(1992,1,1)+timedelta(seconds=total_seconds)
    return(date)


ecco_dir = '/Volumes/helheim/Data_Repository/ECCO'

var_name = 'THETA'
years = [2019]

if var_name=='THETA':
    file_prefix = 'dynDiag'
    file_level = 0

i = np.arange(270)
j = np.arange(270)
k = np.arange(50)
tile = np.arange(13)


for year in years:
    print('Stacking the extension in year '+str(year))
    # find the files for this year
    month_files = 12*['']
    file_dates = 12*['']
    for month in range(1,13):
        for file_name in os.listdir(os.path.join(ecco_dir,'Extension Files',file_prefix)):
            if file_name[:len(file_prefix)]==file_prefix and file_name[-5:]=='.data':
                file_iter = int(file_name.split('.')[-2])
                date = file_iter_to_date(file_iter)
                if date.year == year and date.month==month:
                    month_files[month-1] = file_name
                    file_dates[month-1] = date

    print(month_files)
    # print(file_dates)

    # for each file, read the mds file into tiles
    annual_stack_started = False
    for month in range(1,13):
        file_name=month_files[month-1]
        if file_name!='':
            print('    - Adding file '+str(file_name))
            grid_compact = mds.rdmds(os.path.join(ecco_dir,'Extension Files',file_prefix,file_name[:-5]), returnmeta=False)
            grid_compact = grid_compact[file_level,:,:,:]
            grid_tiles = ec.llc_compact_to_tiles(grid_compact,less_output=True)
            if not annual_stack_started:
                annual_stack_started=True
                annual_stack = grid_tiles[np.newaxis,...]
            else:
                annual_stack = np.concatenate([annual_stack,grid_tiles[np.newaxis,...]])

    time = np.arange(np.shape(annual_stack)[0])


    print('Writing out to nc')
    ds = nc4.Dataset(os.path.join(ecco_dir,var_name,var_name+'_'+str(year)+'.nc'),'w')

    ds.createDimension('i',len(i))
    ds.createDimension('j', len(j))
    ds.createDimension('k', len(k))
    ds.createDimension('tile', len(tile))
    ds.createDimension('time', len(time))

    ivar = ds.createVariable('i','i8',('i',))
    ivar[:] = i

    jvar = ds.createVariable('j', 'i8', ('j',))
    jvar[:] = j

    kvar = ds.createVariable('k', 'i8', ('k',))
    kvar[:] = k

    tilevar = ds.createVariable('tile', 'i8', ('tile',))
    tilevar[:] = tile

    timevar = ds.createVariable('time', 'i8', ('time',))
    timevar[:] = time

    var = ds.createVariable(var_name,'f4',('time','k','tile','j','i'))
    print(np.shape(time),np.shape(k),np.shape(tile),np.shape(j),np.shape(i),np.shape(annual_stack))
    var[:,:,:,:,:] = annual_stack

    ds.close()




