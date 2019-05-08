import os

from SPI_droughts import class_calculate_PET
import pandas as pd


sss = class_calculate_PET.Calculate_PET()

main_dir = os.path.join(r'X:\hiwi\ElHachem\AdvancedPython\SPI_droughts')
main_dir = sss.set_main_dir(main_dir)

in_temp_file = os.path.join(main_dir,
                            r'data\SA_temp_data_1901_2015.csv')
in_temp_min_file = os.path.join(main_dir,
                                r'data\SA_temp_min_var_.csv')
in_temp_max_file = os.path.join(main_dir,
                                r'data\SA_temp_max_var_.csv')
in_ppt_file = os.path.join(main_dir,
                           r'data\SA_pr_1901_2015.csv')

in_ppt_df = sss.read_df(in_ppt_file, 0, ';', '%Y-%m-%d')
in_temp_df = sss.read_df(in_temp_file, 0, ';', '%Y-%m-%d')
in_temp_df_min = sss.read_df(in_temp_min_file, 0, ';', '%Y-%m-%d')
in_temp_df_max = sss.read_df(in_temp_max_file, 0, ';', '%Y-%m-%d')
in_temp_df_min.dropna(inplace=True)
in_temp_df_max.dropna(inplace=True)

mean_temp = (in_temp_df_max.values + in_temp_df_min.values) / 2
mean_temp_Df = pd.DataFrame(index=in_temp_df_max.index, data=mean_temp,
                            columns=['avg'])

pet = [sss.hargreaves_pet(doy, -26.25, tmins, tmaxs, tmean)
       for doy, tmean, tmaxs, tmins in zip(in_temp_df_min.index.dayofyear,
                                           in_temp_df_min.values,
                                           in_temp_df_max.values,
                                           mean_temp_Df.values)]
df_pet = pd.DataFrame(index=in_temp_df_min.index, data=pet)
df_mtl = df_pet.resample('M', label='right', closed='right').sum()


ppt_cmn_idx = in_ppt_df.index.intersection(df_mtl.index)
spdi_vals = in_ppt_df.loc[ppt_cmn_idx, :].values - \
    df_mtl.loc[ppt_cmn_idx, :].values
spdi_df = pd.DataFrame(index=ppt_cmn_idx, data=spdi_vals)

df_mtl.to_csv(os.path.join(main_dir, 'data\SA_pet_var.csv'), sep=';')

spdi_df.to_csv(os.path.join(main_dir, 'data\SA_ppt_min_pet_var.csv'), sep=';')
