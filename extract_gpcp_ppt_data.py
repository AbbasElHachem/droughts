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

main_dir = os.path.join(r'C:\Users\hachem\Desktop\data_Gaza\tmin')
os.chdir(main_dir)

variable_name = 'tmin'  # tmax


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
    main_dir,
    '.nc')
assert len(in_nc_files) > 0, 'no nc files found'

wanted_date_fmt = '%Y-%m-%d'

in_time_name = 'time'

x_coord_1 = 69
y_coord_1 = 117
x_coord_2 = 70
y_coord_2 = 118


idx_ = []
vals = []

for nc_f in in_nc_files:
    in_nc = nc.Dataset(nc_f)

    time_var = in_nc.variables[in_time_name]
#     # convert time array
    time_arr = nc.num2date(time_var[:],
                           time_var.units,
                           calendar='standard')
#                           calendar='standard')

    # read and make time for pandas datetime
    # select temp data
    ppt_vals = in_nc.variables[variable_name][:].data  # change tmax, tmin

    for idx, date in enumerate(time_arr):

        all_ppt_vals_1 = ppt_vals[idx, y_coord_1, x_coord_1]
        all_ppt_vals_2 = ppt_vals[idx, y_coord_2, x_coord_2]

        all_ppt_vals = (all_ppt_vals_1 + all_ppt_vals_2) / 2
        vals.append(all_ppt_vals)
        idx_.append(date)

    in_nc.close()

df_ = pd.DataFrame(data=vals, dtype=float)
str_time = [i.strftime(wanted_date_fmt) for i in idx_]
time_idx = pd.DatetimeIndex(pd.to_datetime(str_time, format=wanted_date_fmt))
df_.columns = [variable_name]  # change tmax, tmin
df_.loc[:, 'time'] = time_idx
df_.set_index('time', drop=True, inplace=True)
# df_monthly = resampleDf(df_, 'M')
df_.to_csv(os.path.join(main_dir, 'gaza_temp_%s_var_.csv' %
                        variable_name), sep=';')

STOP = timeit.default_timer()  # Ending time
print(('\n\a\a\a Done with everything on %s. Total run time was'
       ' about %0.4f seconds \a\a\a' % (time.asctime(), STOP - START)))
