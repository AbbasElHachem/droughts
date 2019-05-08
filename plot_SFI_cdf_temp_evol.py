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
import pandas as pd
warnings.filterwarnings("ignore")
plt.ioff()

print('\a\a\a\a Started on %s \a\a\a\a\n' % time.asctime())
START = timeit.default_timer()  # to get the runtime of the program

main_dir = os.path.join(r'X:\hiwi\ElHachem\AdvancedPython\SPI_droughts')

spi_class = class_calculate_SPI.Calculate_SPI()
main_dir = spi_class.set_main_dir(main_dir)

in_disch_df_file = r'X:\hiwi\ElHachem\code_data_fr_Faizan\Neckar_data\disch_data_with_catch\neckar_daily_discharge_1961_2015.csv'

in_disch_df = spi_class.read_df(in_disch_df_file, 0, ';', '%Y-%m-%d')
in_disch_df = pd.DataFrame(in_disch_df['2446'])
in_disch_df_monthly = spi_class.resampleDf(in_disch_df, 'M', 0)
seasons_dict = {'DJF': [12, 1, 2], 'MAM': [3, 4, 5],
                'JJA': [6, 7, 8], 'SON': [9, 10, 11],
                # in SA ppt in Oct-Apr
                'ONDJFMA': [10, 11, 12, 1, 2, 3, 4],
                'MJJAS': [5, 6, 7, 8, 9],
                'All Year': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                '1M': '1M', '3M': '3M', '6M': '6M', '9M': '9M',
                '12M': '12M', '24M': '24M', '48M': '48M'}

# disch_df_sri = spi_class.calculate_mod_ppt_agg_vals(in_disch_df, 6, '2446')

disch_df_sri = spi_class.calculate_disch_agg_vals(
    in_disch_df, 3, '2446', mod_idx=True)
for i in range(1, 13):
    print(i)

#     spi_class.plot_temp_evol_spi(xw[i], 'Gamma',
#                                  spi_class.calculate_prob_zero_vls(xw[i]),
#                                  'pr', '5M_month_%s' % str(i),
#                                  main_dir, use_fit_fct=False)
    spi_class.plot_temp_evol_spi(disch_df_sri[i], 'Gamma',
                                 spi_class.calculate_prob_zero_vls(
                                     in_disch_df), '2446',
                                 '3M_2446_disch_%s' % str(i),
                                 main_dir, use_fit_fct=False)
raise Exception
# for wtd_season in seasons_dict.keys():
#     if isinstance(seasons_dict[wtd_season], list):
#         ppt_season = spi_class.resampleDfbySeason(in_ppt_df, seasons_dict,
#                                                   wtd_season)
#     if isinstance(seasons_dict[wtd_season], str):
#         temp_freq = seasons_dict[wtd_season]
#         ppt_season = spi_class.resampleDf(in_ppt_df, temp_freq, 0)
#     prob_zero_pp_vls = spi_class.calculate_prob_zero_vls(ppt_season)
#     ppt_season.replace(to_replace=0, value=.01, inplace=True)

#     xw = spi_class.calculate_modified_SPI_idx(in_ppt_df, 5)
#     print(xw)
#     spi_class.plot_orig_fitt_dist(ppt_season, 'Gamma', prob_zero_pp_vls,
#                                   'ppt mm/month', 'cdf', wtd_season,
#                                   main_dir, use_fit_fct=True)
#     plt_spi = spi_class.plot_temp_evol_spi(ppt_season, 'Gamma',
#                                            prob_zero_pp_vls, wtd_season,
#                                            main_dir, use_fit_fct=False)
STOP = timeit.default_timer()  # Ending time
print(('\n\a\a\a Done with everything on %s. Total run time was'
       ' about %0.4f seconds \a\a\a' % (time.asctime(), STOP - START)))
