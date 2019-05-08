# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: EL Hachem Abbas,
Institut f�r Wasser- und Umweltsystemmodellierung - IWS
"""
import os
import time
import timeit

import netCDF4 as nc
import numpy as np
import pandas as pd


print('\a\a\a\a Started on %s \a\a\a\a\n' % time.asctime())
START = timeit.default_timer()  # to get the runtime of the program

main_dir = os.path.join(r'X:\hiwi\ElHachem\AdvancedPython\SPI_droughts\data')
os.chdir(main_dir)


def getFiles(data_dir_, file_ext_str):
    ''' create function to get files based on dir '''
    dfs_files = []
    for r, _, f in os.walk(data_dir_):
        for fs in f:
            if fs.endswith(file_ext_str):
                dfs_files.append(os.path.join(r, fs))
    print('done getting all needed files')
    return dfs_files


def resampleDf(data_frame, temp_freq):
    ''' sample DF based on freq and time shift '''

    df_ = data_frame.copy()
    df_res = df_.resample(temp_freq, label='right', closed='right').sum()
    return df_res


in_nc_files = getFiles(
    r'C:\Users\hachem\Desktop\full_data_monthly_v2018_025.nc',
    '.nc')
assert len(in_nc_files) > 0, 'no nc files found'
wanted_date_fmt = '%Y-%m-%d  %H:%M:%S'
in_time_name = 'time'
x_coord = 233  # idx=224 lat=34.131, value for Lebanon
y_coord = 57  # idx=864 lon=35.90, value for Lebanon

xMin, yMin = 25.2356, -26.1941
xMax, yMax = 33.9993, -19.7914
# subset = ppt_vals.isel(0, lat=slice(439, 466), lon=slice(821, 857)).values
idx_ = []
vals = []
# time_idx = pd.date_range(
#     start='1979',
#     end='22-08-2018',
#     freq='D',
#     closed='left')

for nc_f in in_nc_files:
    in_nc = nc.Dataset(nc_f)

    time_var = in_nc.variables[in_time_name]
#     # convert time array
    time_arr = nc.num2date(time_var[:],
                           time_var.units,
                           calendar=time_var.calendar)
#                           calendar='standard')

    # get the various variables as arrays
#     lat_arr = in_nc.variables['lat'][:]
#     lon_arr = in_nc.variables['lon'][:]

    # read and make time for pandas datetime

    #
    # select ppt data
    ppt_vals = in_nc.variables['tmax']

    for idx, date in enumerate(time_arr):

        all_ppt_vals = ppt_vals[idx, x_coord, y_coord]
        vals.append(all_ppt_vals)
        idx_.append(date)
    in_nc.close()

df_ = pd.DataFrame(data=vals, dtype=float)
str_time = [i.strftime(wanted_date_fmt) for i in idx_]
time_idx = pd.DatetimeIndex(pd.to_datetime(str_time, format=wanted_date_fmt))
df_.columns = ['tmax']
df_.loc[:, 'time'] = time_idx
df_.set_index('time', drop=True, inplace=True)
# df_monthly = resampleDf(df_, 'M')
df_.to_csv(os.path.join(main_dir, 'SA_temp_max_var_.csv'), sep=';')

STOP = timeit.default_timer()  # Ending time
print(('\n\a\a\a Done with everything on %s. Total run time was'
       ' about %0.4f seconds \a\a\a' % (time.asctime(), STOP - START)))
