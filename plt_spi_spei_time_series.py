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
# import numpy as np

print('\a\a\a\a Started on %s \a\a\a\a\n' % time.asctime())
START = timeit.default_timer()  # to get the runtime of the program

main_dir = os.path.join(r'X:\hiwi\ElHachem\AdvancedPython\SPI_droughts')

spi_class = class_calculate_SPI.Calculate_SPI()
main_dir = spi_class.set_main_dir(main_dir)

warnings.filterwarnings("ignore")
plt.ioff()


in_pet_file = os.path.join(main_dir,
                           r'data\SA_pet_var.csv')
in_ppt_minus_pet_file = os.path.join(main_dir,
                                     r'data\SA_ppt_min_pet_var.csv')
oni_df_file = os.path.join(main_dir,
                           r'data\ONI_values.csv')

in_pet_df = spi_class.read_df(in_pet_file, 0, ';', '%Y-%m-%d')
in_ppt_minus_pet_df = spi_class.read_df(
    in_ppt_minus_pet_file, 0, ';', '%Y-%m-%d')
in_oni_df = spi_class.read_df(oni_df_file, 0, ',', '%Y')


seasons_dict = {'DJF': [12, 1, 2], 'MAM': [3, 4, 5],
                'JJA': [6, 7, 8], 'SON': [9, 10, 11],
                # in SA ppt in Oct-Apr
                'ONDJFMA': [10, 11, 12, 1, 2, 3, 4],
                'MJJAS': [5, 6, 7, 8, 9],
                'All Year': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]}

# spi_class.plot_correl_ppt_enso(in_pet_df, in_oni_df, main_dir)
for wtd_season in seasons_dict.keys():

    pet_season = spi_class.resampleDfbySeason(in_pet_df, seasons_dict,
                                              wtd_season)

    in_ppt_minus_pet_df_season = spi_class.resampleDfbySeason(in_ppt_minus_pet_df,
                                                              seasons_dict,
                                                              wtd_season)

    idx_cmn = ppt_season.index.intersection(in_ppt_minus_pet_df_season.index)
    ppt_season_cmn = ppt_season.loc[idx_cmn, :]
    pet_season_cmn = in_ppt_minus_pet_df_season.loc[idx_cmn, :]
    spi_ = spi_class.convert_var_to_idx(ppt_season_cmn.values, 'SPI')
    spei_ = spi_class.convert_var_to_idx(pet_season_cmn.values, 'SPEI')
    prob_zero_pp_vls = spi_class.calculate_prob_zero_vls(ppt_season_cmn)

    spi_class.plot_scatter_plot_two_values(
        spi_, spei_, 'SPI', 'SPEI', wtd_season, main_dir)
    break
#     spi_class.plot_orig_fitt_dist(disch_season, 'Gamma',
#                                   prob_zero_pp_vls,
#                                   'disch m3/month', 'cdf',
#                                   wtd_season, main_dir,
#                                   use_fit_fct=False)


#     spi_class.plot_hist(in_ppt_minus_pet_df_season.values,
#                         'PPT-PET', False, 'PPT-PET', wtd_season, 'r', main_dir)
#     spi_class.plot_hist(pet_season.values,
#                         'PET', False, 'PET', wtd_season, 'm', main_dir)
#     spi_class.plot_hist(ppt_season.values,
# #                         'PPT', False, 'PPT', wtd_season, 'b', main_dir)
# for temp_freq in ['1M', '3M', '6M', '9M', '12M']:
#     ppt_resampled = spi_class.resampleDf(in_ppt_df, temp_freq, 0)
#
#     prob_zero_pp_vls = spi_class.calculate_prob_zero_vls(ppt_resampled)
#
#     ppt_resampled.replace(0, .01, True)
#     in_ppt_minus_pet_df_resampled = spi_class.resampleDf(
#         in_ppt_minus_pet_df, temp_freq, 0)
#
#
# prob_zero_pp_vls = spi_class.calculate_prob_zero_ppt(in_ppt_df)
# plt_fit_orig = spi_class.plot_orig_fitt_dist(in_ppt_df, 'Gamma',
#                                              prob_zero_pp_vls,
#                                              'ppt mm/month', 'cdf',
#                                              'Monthly',
#                                              main_dir, use_fit_fct=True)
# plt_fit_orig = spi_class.plot_orig_spei_fitt_dist(in_ppt_minus_pet_df,
#                                                   'ppt-pet mm/month', 'cdf',
#                                                   'Monthly',
#                                                   main_dir, use_fit_fct=True)
# plt_cycle_orig = spi_class.plot_ppt_anunal_cycle(
#     in_ppt_df, in_pet_df, main_dir)
#
# spi_class.plt_ppt_pet_vals(in_ppt_df, in_pet_df, main_dir)
STOP = timeit.default_timer()  # Ending time
print(('\n\a\a\a Done with everything on %s. Total run time was'
       ' about %0.4f seconds \a\a\a' % (time.asctime(), STOP - START)))
