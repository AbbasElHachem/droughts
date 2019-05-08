# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: EL Hachem Abbas,
Institut f�r Wasser- und Umweltsystemmodellierung - IWS
"""
import os
import time
import timeit
import warnings

from SPI_droughts import class_calculate_SPI
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd


warnings.filterwarnings("ignore")
plt.ioff()

print('\a\a\a\a Started on %s \a\a\a\a\n' % time.asctime())
START = timeit.default_timer()  # to get the runtime of the program

main_dir = os.path.join(r'X:\hiwi\ElHachem\AdvancedPython\SPI_droughts')

spi_class = class_calculate_SPI.Calculate_SPI()
main_dir = spi_class.set_main_dir(main_dir)

seasons_dict = {'DJF': [12, 1, 2], 'MAM': [3, 4, 5],
                'JJA': [6, 7, 8], 'SON': [9, 10, 11],
                # in SA ppt in Oct-Apr
                'ONDJFMA': [10, 11, 12, 1, 2, 3, 4],
                'MJJAS': [5, 6, 7, 8, 9],
                'All Year': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                '1M': '1M', '3M': '3M', '6M': '6M', '9M': '9M',
                '12M': '12M', '24M': '24M', '48M': '48M'}

# in_ppt_file = os.path.join(main_dir,
#                            r'data\out_subset_ppt_data.nc')
in_ppt_file = r'X:\exchange\ElHachem\PPT_data\full_data_monthly_v2018_025.nc'
#
# in_nc_data = xr.open_dataset(in_ppt_file)
# in_time = in_nc_data.__getitem__('time').values
# in_lat = in_nc_data.__getitem__('lat').values
# in_lon = in_nc_data.__getitem__('lon').values
# # in_ppt_data = in_nc_data.__getitem__('__xarray_dataarray_variable__').values
# in_ppt_data = in_nc_data.__getitem__('precip').values
# raise Exception
# time_idx = pd.DatetimeIndex(in_time)


def read_nc_grid_file(nc_file, nc_var_list,
                      time_name, lon_name, lat_name,
                      cut_idx=False):
    '''fct to read nc file and extract ppt data, time, lat, lon'''
    in_nc = nc.Dataset(nc_file)
    # print('extracting lon, lat and data from nc_file: %s' % nc_file)
    for nc_var_name in nc_var_list:
        if nc_var_name in in_nc.variables:
            # print('reading var: %s' % nc_var_name)
            lon = in_nc.variables[lon_name]
            lat = in_nc.variables[lat_name]
            in_var_data = in_nc.variables[nc_var_name]
            time_var = in_nc.variables[time_name]
            try:
                time_arr = nc.num2date(time_var[:],
                                       time_var.units,
                                       calendar='standard')
            except Exception:
                time_arr = nc.num2date(time_var[:],
                                       time_var.units,
                                       calendar=time_var.calendar)
            time_idx = pd.DatetimeIndex(time_arr)
            if cut_idx:  # use it to match time from two NC files
                time_idx_df = pd.DataFrame(index=time_idx)
                time_idx_df = time_idx_df.iloc[709:1512, :]
                in_var_data = in_var_data[709:1512, :, :]
                time_idx = time_idx_df.index
            return lat, lon, in_var_data, time_idx, in_nc


in_lat, in_lon, in_ppt_data, time_idx, in_nc = read_nc_grid_file(in_ppt_file,
                                                                 ['precip'],
                                                                 'time',
                                                                 'lon', 'lat',
                                                                 True)
dr_period = pd.DataFrame()

# points = [(i, j) for i, x in enumerate(in_lon) for j, y in enumerate(in_lat)]
points_x, points_y = np.meshgrid(in_lon, in_lat)
spi_grid = np.empty(shape=(1, in_lat.shape[0], in_lon.shape[0]))
for latx, _ in enumerate(in_lat):
    print('going through lat idx:', latx)
    for lonx, _ in enumerate(in_lon):
        print('going through lon idx:', lonx)
        try:
            df = pd.DataFrame(index=time_idx)
            try:
                ppt_row = np.array(
                    in_ppt_data[:, latx, lonx].data, dtype='float64')
            except Exception:
                ppt_row = np.array(in_ppt_data[:, latx, lonx], dtype='float64')
            df['ppt'] = ppt_row
            xvls, edf_vls = spi_class.build_edf_fr_vals(df['ppt'].values)
            ppt_vls, edf_ppt = spi_class.get_vals_fr_edf_vals(
                df, 'ppt', edf_vls)
            df.loc[:, 'edf'] = edf_ppt

            ppt_v, edf_v, spi_v = spi_class.get_vals_fr_idx_vals(
                df, 'ppt', edf_vls)
            df.loc[:, 'spi'] = spi_v
            df = spi_class.resampleDfbySeason(df, seasons_dict, 'MAM')
#             print(df)
            df.fillna(0, inplace=True)
# #             print(df.loc[:, 'spi'].values.mean())
            print(df.index)
            spi_grid[:,
                     latx,
                     lonx] = df.loc['1950-03-01':'1960-03-01',
                                    'spi'].values.mean()
        except Exception as msg:
            spi_grid[:, latx, lonx] = 0
            print(msg)
            continue
        break
    break
#
msh = plt.pcolormesh(points_x, points_y, spi_grid)
plt.colorbar(msh)
plt.savefig('year1960mam.png')
plt.show()  # http://www.fao.org/docrep/008/y5744e/y5744e05.htm

# for wtd_season in seasons_dict.keys():
#     if isinstance(seasons_dict[wtd_season], list):
#         ppt_season = spi_class.resampleDfbySeason(in_ppt_df, seasons_dict,
#                                                   wtd_season)
#     if isinstance(seasons_dict[wtd_season], str):
#         temp_freq = seasons_dict[wtd_season]
#         ppt_season = spi_class.resampleDf(in_ppt_df, temp_freq, 0)
#     prob_zero_pp_vls = spi_class.calculate_prob_zero_vls(ppt_season)
#     ppt_season.replace(to_replace=0, value=.01, inplace=True)
#
#     spi_class.plot_orig_fitt_dist(ppt_season, 'Gamma', prob_zero_pp_vls,
#                                   'ppt mm/month', 'cdf', wtd_season,
#                                   main_dir, use_fit_fct=True)
#     plt_spi = spi_class.plot_temp_evol_spi(ppt_season, 'Gamma',
#                                            prob_zero_pp_vls, wtd_season,
#                                            main_dir, use_fit_dist=False)
STOP = timeit.default_timer()  # Ending time
print(('\n\a\a\a Done with everything on %s. Total run time was'
       ' about %0.4f seconds \a\a\a' % (time.asctime(), STOP - START)))
