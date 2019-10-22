import os


import pandas as pd


def hargreaves_pet(d_o_y, lat_, t_min, t_max, t_avg):
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


main_dir = os.path.join(r'C:\Users\hachem\Desktop\data_Gaza')
#main_dir = class_calculate_SPI.set_main_dir(main_dir)


in_temp_min_file = os.path.join(main_dir,
                                r'tmin\gaza_temp_min_var_.csv')
in_temp_max_file = os.path.join(main_dir,
                                r'tmax\gaza_temp_max_var_.csv')


# read df min
in_temp_df_min = pd.read_csv(in_temp_min_file, index_col=0,
                             sep=';', parse_dates=True,
                             infer_datetime_format=True)

# read df max
in_temp_df_max = pd.read_csv(in_temp_max_file, index_col=0,
                             sep=';', parse_dates=True,
                             infer_datetime_format=True)
# calculate df avg
mean_temp = (in_temp_df_max.values + in_temp_df_min.values) / 2
# make it a dataframe
mean_temp_Df = pd.DataFrame(index=in_temp_df_max.index, data=mean_temp,
                            columns=['avg'])

pet = [hargreaves_pet(doy, 31.443895, tmins, tmaxs, t_avg)
       for doy, t_avg, tmaxs, tmins in zip(in_temp_df_min.index.dayofyear,
                                           in_temp_df_min.values,
                                           in_temp_df_max.values,
                                           mean_temp_Df.values)]
df_pet = pd.DataFrame(index=in_temp_df_min.index, data=pet)

# df_mtl = df_pet.resample('M', label='right', closed='right').sum()


# ppt_cmn_idx = in_ppt_df.index.intersection(df_mtl.index)
# spdi_vals = in_ppt_df.loc[ppt_cmn_idx, :].values - \
#     df_mtl.loc[ppt_cmn_idx, :].values
# spdi_df = pd.DataFrame(index=ppt_cmn_idx, data=spdi_vals)

df_pet.to_csv(os.path.join(main_dir, 'gaza_pet_var.csv'), sep=';')

# spdi_df.to_csv(os.path.join(main_dir, 'data\SA_ppt_min_pet_var.csv'), sep=';')
