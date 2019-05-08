# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: EL Hachem Abbas,
Institut fuer Wasser- und Umweltsystemmodellierung - IWS
"""


class Calculate_PET(object):

    def calculate_PET_values(self, mean_temp_df):
        '''Calculate Potential Evapotraspiration based on daily Temperature'''
        import sys
        if 'pandas' not in sys.modules:
            import pandas as pd
        df_etp = mean_temp_df.copy()
        df_etp.replace(df_etp.values < 0, 0, inplace=True)

        monthly_idx = df_etp.resample('M', label='right',
                                      closed='right').mean()
        for i, tmp_val in zip(monthly_idx.index, monthly_idx.values):
            monthly_idx.loc[i, 'Monthly_heat_idx'] = (tmp_val / 5)**1.514

        yearly_etp = monthly_idx['Monthly_heat_idx'].resample('Y', label='right',
                                                              closed='right').sum()
        yearly_etp.index = pd.to_datetime(yearly_etp.index, format='%Y')
        yearly_etp = pd.DataFrame(data=yearly_etp.values,
                                  index=yearly_etp.index,
                                  columns=['yearly_idx'])

        for year, val in zip(yearly_etp.index, yearly_etp.values):
            yearly_etp.loc[year, 'alpha'] = 675 * 10e-9 * \
                val**3 - 771.10e-7 * val**2 + 1792e-5 * val + 0.49239
        sun_hrs = {1: 8.17, 2: 8.02, 3: 7.42, 4: 8.42, 5: 9.19, 6: 9.12,
                   7: 9.23, 8: 9.54, 9: 9.18, 10: 8.54, 11: 9, 12: 9.05}
        for year in yearly_etp.index:
            alpha = yearly_etp.loc[year, 'alpha']
            yearly_idx = yearly_etp.loc[year, 'yearly_idx']
            for month, temp in zip(monthly_idx.index,
                                   monthly_idx['avg'].values):

                monthly_idx.loc[month, 'PET'] = 16 * (sun_hrs[month.month]
                                                      / 12) * (
                    month.day / 30) * (
                        10 * temp / yearly_idx) ** alpha
        return monthly_idx

    def hargreaves_pet(self, d_o_y, lat_, t_min, t_max, t_avg):
        """
        Purpose: To get the potential evapotranspiration at a given latitude \
                    for a given date and temperature.
        Description of the arguments:
            d_o_y (int): day of the year.
            lat (degrees): latitude of the point.
            t_min (celsius): minimum temperature on that day
            t_max (celsius): maximum temperature on that day
            t_avg (celsius): average temperature on that day
        """
        from numpy import radians, sin, pi
        from math import acos, tan, cos, sqrt
        tot_days = 365
        lat = radians(lat_)
        ndec = 0.409 * sin(((2 * pi * d_o_y) / tot_days) - 1.39)

        nws = acos(-tan(lat) * tan(ndec))
        dfr = 1 + (0.033 * cos((2 * pi * d_o_y) / tot_days))
        fac_1 = 15.342618001389575  # ((1440 * 0.082 * 0.4082)/pi)
        ra = fac_1 * dfr * ((nws * sin(lat) * sin(ndec)) +
                            (cos(lat) * cos(ndec) * sin(nws)))
        fac_2 = 0.002295  # (0.0135 * 0.17)
        pet = fac_2 * ra * sqrt(t_max - t_min) * (t_avg + 17.8)

        if pet > 0:
            return pet
        else:
            return 0
