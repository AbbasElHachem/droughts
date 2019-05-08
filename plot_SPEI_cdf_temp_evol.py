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

warnings.filterwarnings("ignore")
plt.ioff()

print('\a\a\a\a Started on %s \a\a\a\a\n' % time.asctime())
START = timeit.default_timer()  # to get the runtime of the program

main_dir = os.path.join(r'X:\hiwi\ElHachem\AdvancedPython\SPI_droughts')

spi_class = class_calculate_SPI.Calculate_SPI()
main_dir = spi_class.set_main_dir(main_dir)

in_ppt_minus_pet_file = os.path.join(main_dir,
                                     r'data\SA_ppt_min_pet_var.csv')
in_ppt_minus_pet_df = spi_class.read_df(
    in_ppt_minus_pet_file, 0, ';', '%Y-%m-%d')


seasons_dict = {'DJF': [12, 1, 2], 'MAM': [3, 4, 5],
                'JJA': [6, 7, 8], 'SON': [9, 10, 11],
                # in SA ppt in Oct-Apr
                'ONDJFMA': [10, 11, 12, 1, 2, 3, 4],
                'MJJAS': [5, 6, 7, 8, 9],
                'All Year': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                '1M': '1M', '3M': '3M', '6M': '6M', '9M': '9M',
                '12M': '12M', '24M': '24M', '48M': '48M'}
for wtd_season in seasons_dict.keys():
    if isinstance(seasons_dict[wtd_season], list):
        in_ppt_minus_pet_df_season = spi_class.resampleDfbySeason(in_ppt_minus_pet_df,
                                                                  seasons_dict,
                                                                  wtd_season)
    if isinstance(seasons_dict[wtd_season], str):
        temp_freq = seasons_dict[wtd_season]
        in_ppt_minus_pet_df_season = spi_class.resampleDf(
            in_ppt_minus_pet_df, temp_freq, 0)

    prob_zero_pp_vls = spi_class.calculate_prob_zero_vls(
        in_ppt_minus_pet_df_season)
    in_ppt_minus_pet_df_season.replace(to_replace=0, value=.01, inplace=True)
#     be, al, gm = spi_class.lmoments_method(in_ppt_minus_pet_df_season)
#     pdf_log = spi_class.log_logistic_pdf(in_ppt_minus_pet_df_season, be, al, gm)
    spi_class.plot_orig_spei_fitt_dist(in_ppt_minus_pet_df_season,
                                       'PPT-PET (mm/month)', 'cdf', wtd_season,
                                       main_dir, use_fit_fct=True)
    spi_class.plot_temp_evol_spei(in_ppt_minus_pet_df_season,
                                  wtd_season, main_dir)
STOP = timeit.default_timer()  # Ending time
print(('\n\a\a\a Done with everything on %s. Total run time was'
       ' about %0.4f seconds \a\a\a' % (time.asctime(), STOP - START)))
