# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: EL Hachem Abbas,
Institut für Wasser- und Umweltsystemmodellierung - IWS
"""

import os

from scipy.special import gamma

import FitKernelDensityFunction
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as optimize
import scipy.stats as stats


class Calculate_SPI(object):
    def set_main_dir(self, main_dir):
        ''' used to create and set main working directory '''
        if not os.path.exists(main_dir):
            os.mkdir(main_dir)
        os.chdir(main_dir)
        return main_dir

    def read_df(self, df_file, idx_col, sep, date_fmt=None):
        ''' fct to read adtaframe, define index column,
            column seperator and date format for the index if time index'''
        assert os.path.exists(df_file), 'df_file is not correct'
        in__df = pd.read_csv(df_file, index_col=idx_col, sep=sep)
        if date_fmt is not None:
            in__df.index = pd.to_datetime(in__df.index, format=date_fmt)
        return in__df

    def load_pickle(self, in_file, mode='rb'):
        ''' fct to read ppt, temp,pet data in a dict pickle file'''
        import pickle
        with open(in_file, mode) as _pkl_hdl:
            return pickle.load(_pkl_hdl)
        return

    def resampleDf(self, df_, temp_freq, temp_shift):
        ''' resample DF based on freq and time shift and the sum '''
        assert isinstance(df_, pd.DataFrame), 'data is not a df'
        return df_.resample(temp_freq, label='right', closed='right',
                            base=temp_shift).sum()

    def resampleDfbySeason(self, data_frame, dict_seasons, chosen_season):
        ''' resamaple DF based on a season, defined by months in dict (J: 1)'''
        assert isinstance(data_frame, pd.DataFrame), 'data is not a df'
        return pd.concat([data_frame[data_frame.index.month == month]
                          for season in dict_seasons.keys()
                          for month in dict_seasons[season]
                          if chosen_season == season], axis=0,
                         ignore_index=False, sort=True)

    def calculate_ppt_agg_vals(self, ppt_df_monthly, w_month,
                               df_out_col_name_str, mod_idx=False):
        ''' D(t) represent the rainfall depth measured at time t,
            Xw(t) = sum (fr:i=t-w+1, to t) D(i) is the total precipitation.
            for a given w-month window with respect to t.
        '''
        assert isinstance(ppt_df_monthly, pd.DataFrame), 'Data should be a df'
        assert isinstance(w_month, int or float), 'window size should be a nbr'
        x_w_sr = []
        for month_ in range(w_month - 1, len(ppt_df_monthly.values)):
            low_intv = month_ - w_month + 1
            x_w, x_w_gr = [], []
            try:
                while low_intv <= month_:
                    x_w.append(ppt_df_monthly.iloc[low_intv, :].values)
                    low_intv += 1
                x_w_gr.append(np.sum(x_w))
            except Exception as msg:
                print(msg)
            x_w_sr.append(x_w_gr)
        df_agg_ptt = pd.DataFrame(data=x_w_sr,
                                  index=ppt_df_monthly.index[w_month - 1:],
                                  columns=[df_out_col_name_str])
        if mod_idx:
            x_w_m = {m_idx: df_agg_ptt[df_agg_ptt.index.month == m_idx]
                     for m_idx in range(1, 13)}
            return x_w_m
        return df_agg_ptt

    def calculate_disch_agg_vals(self, disch_df_monthly, w_month,
                                 df_out_col_name_str, mod_idx=False):
        ''' R(t) represent the discharge measured at time t,
            Yw(t) = sum (fr:i=t-w+1, to t) R(i)/W is the mean discharge
            for a given w-month window with respect to t.
        '''
        assert isinstance(disch_df_monthly, pd.DataFrame), 'in data not a df'
        assert isinstance(w_month, int or float), 'window size should be a nbr'
        y_w_sr = []
        for month_ in range(w_month - 1, len(disch_df_monthly.values)):
            low_intv = month_ - w_month + 1
            y_w, y_w_gr = [], []
            try:
                while low_intv <= month_:
                    y_w.append(disch_df_monthly.iloc[low_intv, :].values)
                    low_intv += 1
                y_w_gr.append(np.sum(y_w) / w_month)
            except Exception as msg:
                print(msg)
            y_w_sr.append(y_w_gr)
        df_agg_disch = pd.DataFrame(data=y_w_sr,
                                    index=disch_df_monthly.index[w_month - 1:],
                                    columns=[df_out_col_name_str])
        if mod_idx:
            y_w_m = {m_idx: df_agg_disch[df_agg_disch.index.month == m_idx]
                     for m_idx in range(1, 13)}
            return y_w_m
        return df_agg_disch

    def gamma_of_ppt_data(self, ppt_data, *params):
        '''calculate log_pdf of Gamma distribution for Precipitation data'''
        return stats.gamma(a=params[0], scale=params[1]).logpdf(ppt_data)

    def pearson3_of_ppt_data(self, ppt_data, *params):
        '''calculate log_pdf of Pearson3 distribution for Precipitation data'''
        return stats.pearson3(skew=params[0], scale=params[1]).logpdf(ppt_data)

    def likelihood(self, params, _data, dist_fct):
        '''calculate the sum of the log_pdf of all values'''
        return -np.sum(np.vectorize(dist_fct)(_data, params[0], params[1]))

    def log_logistic_pdf(self, data, beta, alfa, gamma):
        '''calculate loglogistic density function'''
        return (beta / alfa) * (((data - gamma) / alfa)**(beta - 1)) * \
            (1 + ((data - gamma) / alfa)**beta)**-2

    def log_logistic_cdf(self, data, beta, alfa, gamma):
        '''calculate CDF of Log_logisitc function'''
        return (1 + (alfa / (data - gamma))**beta)**-1

    def optimize_likelihood(self, _data, dist_type):
        ''' minimize the sum and find parameters of used distribution'''
        if dist_type == 'Gamma':
            use_ftn = self.gamma_of_ppt_data
        if dist_type == 'Pearson3':
            use_ftn = self.pearson3_of_ppt_data

        opt_prams = optimize.minimize(self.likelihood, args=(_data, use_ftn),
                                      x0=[1, 1], method='SLSQP',
                                      bounds=[(0., None), (0., None)])
        return opt_prams.x

    def lmoments_method(self, data):
        ''' estimate paramters of Loglogistic function using Lmoments method'''
        data = np.squeeze(data)
        w0, w1, w2 = np.mean(data), np.var(data), stats.skew(data)
        _beta = (2 * w1 - w0) / (6 * w1 - w0 - 6 * w2)
        _alfa = ((w0 - 2 * w1) * _beta) / \
            (gamma(1 + 1 / _beta) * gamma(1 - 1 / _beta))
        _gamma = w0 - _alfa * \
            gamma(1 + 1 / _beta) * \
            gamma(1 - 1 / _beta)
        return _beta, _alfa, _gamma

    def build_edf_fr_vals(self, data):
        '''construct empirical distribution function given data values'''
        data_sorted = np.sort(data, axis=0)[::-1]
        x0 = np.squeeze(data_sorted)[::-1]
        y0 = (np.arange(data_sorted.size) / len(data_sorted))
        return x0, y0

    def get_vals_fr_edf_vals(self, in_df, wtd_col_name, edf_vals):
        '''construct edf and get indices of each corresponding value in DF'''

        in_df_edf = in_df.copy()
        assert wtd_col_name in in_df.columns
        in_df_edf['sorted_vls'] = np.sort(in_df[wtd_col_name].values)
        in_df_edf['edf_vls'] = edf_vals
        indices = []
        for p in in_df[wtd_col_name].values:
            for j, ps in enumerate(in_df_edf['sorted_vls'].values):
                if p == ps:
                    if j not in indices:
                        indices.append(j)

        edf_ppt_ = in_df_edf['edf_vls'].values[indices]
        return in_df_edf[wtd_col_name].values, edf_ppt_

    def calculate_prob_zero_vls(self, df_orig):
        '''calcualte the Probability of 0 rainfall values'''
        return (len(df_orig[df_orig.values == 0].values) / len(df_orig.values))

    def build_cdf_fitted_fct(self, dist_type, ppt_data, params, p0):
        '''construct theoretical CDF based on chosen distribution and params'''
        if dist_type == 'Gamma':
            cdf_vls = stats.gamma.cdf(ppt_data.values, a=params[0],
                                      scale=params[1])
        if dist_type == 'Pearson3':
            cdf_vls = stats.pearson3.cdf(ppt_data.values, skew=params[0],
                                         scale=params[1])
        return p0 + (1 - p0) * cdf_vls

    def calculate_SPI_fr_cdf(self, fitted_cdf):
        '''Calculate Standard Precipitation Index from CDF'''
        c0, c1, c2 = 2.515517, 0.802853, 0.010328
        d1, d2, d3 = 1.432788, 0.189269, 0.001308

        SPI_total = np.empty(shape=(len(fitted_cdf)))
        for i, val in enumerate(fitted_cdf):
            try:
                if 0 < val <= 0.5:
                    t = np.sqrt(np.log(1 / val**2))
                    SPI_total[i] = -(t - (c0 + c1 * t + c2 * t**2) /
                                     (1 + d1 * t + d2 * t**2 + d3 * t**3))
                elif 0.5 < val < 1:
                    t = np.sqrt(np.log(1 / (1 - val)**2))
                    SPI_total[i] = +(t - (c0 + c1 * t + c2 * t**2) /
                                     (1 + d1 * t + d2 * t**2 + d3 * t**3))
            except Exception:
                SPI_total[i] = np.nan
            if np.isclose(SPI_total[i], 0, atol=1e-2):
                SPI_total[i] = np.nan
            if SPI_total[i] < - 10 or 10 <= SPI_total[i]:
                SPI_total[i] = np.nan
        return SPI_total

    def calculate_SPEI_fr_cdf(self, cdf_vals):
        '''Calculate Standard Precipitation Evapostranspiration Index from CDF'''
        # https://journals.ametsoc.org/doi/10.1175/2009JCLI2909.1
        SPEI = np.empty(shape=(len(cdf_vals)))
        c0, c1, c2 = 2.515517, 0.802853, 0.010328
        d1, d2, d3 = 1.432788, 0.189269, 0.001308
        for i, pval in enumerate(cdf_vals):
            try:
                if 0 < pval <= 0.5:
                    w = np.sqrt(-2 * np.log(1 - pval))
                    SPEI[i] = (w - ((c0 + c1 * w + c2 * w**2) /
                                    (1 + d1 * w + d2 * w**2 + d3 * w**3)))
                if pval > 0.5:
                    w = np.sqrt(-2 * np.log(pval))
                    SPEI[i] = -(w - ((c0 + c1 * w + c2 * w**2) /
                                     (1 + d1 * w + d2 * w**2 + d3 * w**3)))
            except Exception:
                SPEI[i] = np.nan
            if np.isclose(SPEI[i], 0, atol=1e-2):
                SPEI[i] = np.nan
            if SPEI[i] < - 7 or 7 <= SPEI[i]:
                SPEI[i] = np.nan
        return SPEI

    def spi_cdf_fit_dist(self, ppt_vals, prob_zero_ppt,
                         dist_type, use_fit_ftn):
        ''' transform from cdf to spi using a fit optimize function '''
        x_sim, y_sim = self.build_edf_fr_vals(ppt_vals)
        if use_fit_ftn:
            dist_params = self.optimize_likelihood(ppt_vals, dist_type)
            x_sim, y_sim = ppt_vals, self.build_cdf_fitted_fct(
                dist_type, ppt_vals, dist_params, prob_zero_ppt)
        spi_ = self.calculate_SPI_fr_cdf(y_sim)
        y_sim_spi = np.delete(y_sim, np.where(np.isnan(spi_)))
        spi_ = spi_[~np.isnan(spi_)]
        return x_sim, y_sim, spi_, y_sim_spi

    def spei_cdf_fit_dist(self, ppt_min_pet_vals):
        ''' transform from cdf to spei using a kernel fit'''
        fitkernel = FitKernelDensityFunction.FitKernelDensity()
        ppt_min_pet_sorted_vals = np.sort(ppt_min_pet_vals.values)
        d_opt = fitkernel.kernel_width_optimization(ppt_min_pet_sorted_vals,
                                                    'gauss', log=False)
        xticks, F_x = fitkernel.max_likelyhood_cdf(ppt_min_pet_sorted_vals,
                                                   d_opt, 'gauss',
                                                   steps=len(
                                                       ppt_min_pet_sorted_vals),
                                                   plot=False, log=False)
        spei_ = self.calculate_SPEI_fr_cdf(F_x)
        return xticks, F_x, spei_

    def convert_var_to_idx(self, values, idx_type):
        ''' convert values (ppt or ppt-pet values to SPI or SPEI values'''
        _, y_orig = self.build_edf_fr_vals(values)
        if idx_type == 'SPI':
            idx_ = self.calculate_SPI_fr_cdf(y_orig)
        if idx_type == 'SPEI':
            idx_ = self.calculate_SPEI_fr_cdf(y_orig)

        _ = np.delete(y_orig, np.where(np.isnan(idx_)))
        idx_final = idx_[~np.isnan(idx_)]
        return idx_final

    def get_vals_fr_idx_vals(self, in_df, wtd_col_name, edf_vals):
        '''construct edf and get indices of each corresponding value in DF'''
        in_df_edf = in_df.copy()
        assert wtd_col_name in in_df.columns
        in_df_edf['sorted_vls'] = np.sort(in_df[wtd_col_name].values)
        in_df_edf['edf_vls'] = edf_vals
        in_df_edf['spi_vls'] = self.calculate_SPI_fr_cdf(edf_vals)
        indices = []
        for p in in_df[wtd_col_name].values:
            for j, ps in enumerate(in_df_edf['sorted_vls'].values):
                if p == ps:
                    if j not in indices:
                        indices.append(j)
        edf_ppt_ = in_df_edf['edf_vls'].values[indices]
        spi_ppt_ = in_df_edf['spi_vls'].values[indices]
        return in_df_edf[wtd_col_name].values, edf_ppt_, spi_ppt_

    def adjust_subplots_labels(self, sub_plt1, x_label, y_label):
        ''' adjust labels and axis of a subplot'''
        sub_plt1.legend(loc=0)
        sub_plt1.set_xlabel(x_label)
        sub_plt1.set_ylabel(y_label)
        sub_plt1.grid(alpha=0.6)
        return sub_plt1

    def plot_orig_fitt_dist(self, ppt_orig_vals, dist_type,
                            prob_zero_ppt, x_label, y_label,
                            season, out_dir, use_fit_fct=False):
        '''2 subplots for Fitted CDF and EDF and SPI for Precipitation data'''
        x_orig, y_orig = self.build_edf_fr_vals(ppt_orig_vals)

        spi_edf = self.calculate_SPI_fr_cdf(y_orig)

        y_orig_spi = np.delete(y_orig, np.where(np.isnan(spi_edf)))
        spi_edf = spi_edf[~np.isnan(spi_edf)]

        _, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8),
                                     dpi=300, sharey=True)
        ax1.scatter(x_orig, y_orig, label='PPT_edf', c='b',
                    alpha=0.5, s=15, marker='o')
        ax2.scatter(spi_edf, y_orig_spi, c='b', label='SPI_edf',
                    marker='o', alpha=0.25, s=12)

        if use_fit_fct:
            x_sim, y_sim, spi_, y_sim_spi = self.spi_cdf_fit_dist(ppt_orig_vals,
                                                                  prob_zero_ppt,
                                                                  dist_type,
                                                                  use_fit_fct)
            ax1.scatter(x_sim, y_sim, c='r', alpha=0.6, s=13, marker='+',
                        label='Ftted %s Dist' % dist_type)
            ax2.scatter(spi_, y_sim_spi, c='r', label='SPI',
                        marker='+', alpha=0.25, s=14)

        ax1 = self.adjust_subplots_labels(ax1, x_label, y_label)
        ax2 = self.adjust_subplots_labels(ax2, 'SPI', y_label)

        plt.savefig(os.path.join(out_dir,
                                 'fitted_%s_to_cdf_season_%s_ppt_data_SA.png'
                                 % (dist_type, season)),
                    frameon=True, papertype='a4', bbox_inches='tight')
        return x_orig, y_orig, spi_edf

    def plot_orig_spei_fitt_dist(self, ppt_min_pet_vals, x_label, y_label,
                                 season, out_dir, use_fit_fct=False):
        '''2 subplots for Fitted CDF and EDF and SPEI for ETP data'''
        x_orig, y_orig = self.build_edf_fr_vals(ppt_min_pet_vals)
        spei_orig = self.calculate_SPEI_fr_cdf(y_orig)

        y_orig_spei = np.delete(y_orig, np.where(np.isnan(spei_orig)))
        spei_orig_ = spei_orig[~np.isnan(spei_orig)]

        _, (ax1, ax2) = plt.subplots(
            1, 2, figsize=(12, 8), dpi=300, sharey=True)
        ax1.scatter(x_orig, y_orig, label='PPT-PET edf', c='b',
                    alpha=0.5, s=15, marker='o')
        ax2.scatter(spei_orig_, y_orig_spei, c='m', label='SPEI Orig',
                    marker='+', alpha=0.5, s=14)

        if use_fit_fct:
            xticks, F_x, spei_ = self.spei_cdf_fit_dist(ppt_min_pet_vals)
            F_x_ = np.delete(F_x, np.where(np.isnan(spei_)))
            _spei_ = spei_[~np.isnan(spei_)]

            ax1.scatter(xticks, F_x, c='r', alpha=0.6, s=25, marker='o',
                        label='Ftted Gauss Kernel')
            ax2.scatter(_spei_, F_x_, c='r', label='SPEI Kernel',
                        marker='+', alpha=0.4, s=20)
#             alfa_, beta_, gamma_ = self.lmoments_method(ppt_min_pet_vals)
#             print(alfa_, beta_, gamma_)
#             xticks, F_x = np.sort(
#                 ppt_min_pet_vals), self.log_logistic_cdf(ppt_min_pet_vals,
#                                                          alfa_, beta_, gamma_)
#             spei_ = self.calculate_SPI_fr_cdf(F_x)
#             print(spei_, F_x)
#             F_x_ = np.delete(F_x, np.where(np.isnan(spei_)))
#             _spei_ = spei_[~np.isnan(spei_)]
#
#             ax1.scatter(xticks, F_x, c='r', alpha=0.6, s=25, marker='o',
#                         label='Ftted Gauss Kernel')
#             ax2.scatter(_spei_, F_x_, c='r', label='SPEI Kernel',
#                         marker='+', alpha=0.4, s=20)

        ax1 = self.adjust_subplots_labels(ax1, x_label, y_label)
        ax2 = self.adjust_subplots_labels(ax2, 'SPEI', y_label)
        plt.savefig(os.path.join(out_dir,
                                 'spei_cdf_season_%s_pet_data_SA.png'
                                 % (season)),
                    frameon=True, papertype='a4', bbox_inches='tight')
        return y_orig, spei_orig_

    def plot_temp_evol_spi(self, ppt_orig_vals, dist_type,
                           prob_zero_ppt, df_col_name,
                           season, out_dir,
                           use_fit_fct=False):
        '''Plot of calculated SPI along time'''
        _, y_orig = self.build_edf_fr_vals(ppt_orig_vals)
        _, edf_not_sorted = self.get_vals_fr_edf_vals(
            ppt_orig_vals, df_col_name, y_orig)

        spi_ = self.calculate_SPI_fr_cdf(edf_not_sorted)
        spi_ = spi_[~np.isnan(spi_)]

        pos_inds, neg_inds = np.where(spi_ >= 0.)[0], np.where(spi_ < 0.)[0]
        time = ppt_orig_vals.index

        years, b_width = mdates.YearLocator(), 30
        fig, ax1 = plt.subplots(1, 1, figsize=(32, 8))
        ax1.bar(time[pos_inds], spi_[pos_inds],
                width=b_width, align='center', color='b')
        ax1.bar(time[neg_inds], spi_[neg_inds],
                width=b_width, align='center', color='r')
        if use_fit_fct:
            _, _, spi_fit, _ = self.spi_cdf_fit_dist(ppt_orig_vals,
                                                     prob_zero_ppt,
                                                     dist_type,
                                                     use_fit_fct)
            spi_fit = spi_fit[~np.isnan(spi_fit)]

            pos_inds_fit = np.where(spi_fit >= 0.)[0]
            neg_inds_fit = np.where(spi_fit < 0.)[0]
            ax1.bar(time[pos_inds_fit], spi_fit[pos_inds_fit],
                    width=b_width, align='center', color='g')
            ax1.bar(time[neg_inds_fit], spi_fit[neg_inds_fit],
                    width=b_width, align='center', color='m')
        ax1.grid(True, alpha=0.5)
        ax1.set_xlabel("Time")
        ax1.xaxis.set_major_locator(years)
        ax1.set_xlim(time[0], time[-1])
        ax1.format_xdata = mdates.DateFormatter('%Y')
        ax1.set_yticks([-3, -2, -1.6, -1.3, -0.8,
                        -0.5, 0, 0.5, 0.8, 1.3, 1.6, 2, 3])
        ax1.set_ylim([-5, 5])
        ax1.set_ylabel('SPI')
        fig.autofmt_xdate()
        plt.savefig(os.path.join(out_dir,
                                 'temp_season_%s_spi_data_SA.png'
                                 % (season)),
                    frameon=True, papertype='a4', bbox_inches='tight')
        return spi_

    def plot_temp_evol_spei(self, pet_vals, season, out_dir):
        '''Plot of calculated SPEI along time'''
        _, y_orig = self.build_edf_fr_vals(pet_vals)
        _, edf_not_sorted = self.get_vals_fr_edf_vals(
            pet_vals, '0', y_orig)

        spei_ = self.calculate_SPEI_fr_cdf(edf_not_sorted)
        spei_ = spei_[~np.isnan(spei_)]
        pos_spei, neg_spei = np.where(spei_ >= 0.)[0], np.where(spei_ < 0.)[0]

        time_spei = pet_vals.index
        years, b_width = mdates.YearLocator(), 30

        fig, ax2 = plt.subplots(1, 1, figsize=(32, 8))
        ax2.bar(time_spei[pos_spei], spei_[pos_spei],
                width=b_width, align='center', color='m')
        ax2.bar(time_spei[neg_spei], spei_[neg_spei],
                width=b_width, align='center', color='c')
        ax2.grid(True, alpha=0.5)
        ax2.set_xlabel("Time")
        ax2.xaxis.set_major_locator(years)
        ax2.set_xlim(time_spei[0], time_spei[-1])
        ax2.format_xdata = mdates.DateFormatter('%Y')
        ax2.set_yticks([-3, -2, -1.6, -1.3, -0.8,
                        -0.5, 0, 0.5, 0.8, 1.3, 1.6, 2, 3])
        ax2.set_ylim([-5, 5])
        ax2.set_ylabel('SPEI')
        fig.autofmt_xdate()

        plt.savefig(os.path.join(out_dir,
                                 'temp_season_%s_spei_data_SA.png'
                                 % (season)),
                    frameon=True, papertype='a4', bbox_inches='tight')
        return spei_, time_spei

    def agg_DF_max_mean_min_cycle(self, data_frame):
        ''' group a Dataframe monthly, to find maxn min and avg yearly cycle'''
        df_cycle_max = data_frame.groupby([data_frame.index.month]).max()
        df_cycle_mean = data_frame.groupby([data_frame.index.month]).mean()
        df_cycle_min = data_frame.groupby([data_frame.index.month]).min()
        assert df_cycle_max.index == df_cycle_mean.index == df_cycle_min.index
        idx = df_cycle_max.index
        return idx, df_cycle_max.values, df_cycle_mean.values, df_cycle_min.values

    def do_line_plot_x_y(self, subplot, xvals, yvals, label_, color_):
        ''' plot line plot of x and y and define label and color'''
        subplot.plot(xvals, yvals, label=label_, color=color_, alpha=0.7)
        subplot.grid(True)
        subplot.set_xticks([i for i in xvals])
        subplot.set_xlabel("month")
        subplot.set_ylabel('mm')
        plt.legend(loc=0)
        return subplot

    def plot_ppt_anunal_cycle(self, in_ppt_df, in_pet_df, out_dir):
        '''Plot of PET and PPT annual cycle'''
        ix_ppt, ppt_max, ppt_mean, ppt_min = self.agg_DF_max_mean_min_cycle(
            in_ppt_df)
        ix_pet, pet_max, pet_mean, pet_min = self.agg_DF_max_mean_min_cycle(
            in_pet_df)
        _, ax = plt.subplots(figsize=(20, 10))

        self.do_line_plot_x_y(ax, ix_ppt, ppt_max, 'PPT Max', 'r')
        self.do_line_plot_x_y(ax, ix_ppt, ppt_mean, 'PPT Mean', 'b')
        self.do_line_plot_x_y(ax, ix_ppt, ppt_min, 'PPT Min', 'c')
        self.do_line_plot_x_y(ax, ix_pet, pet_max, 'PET Max', 'r')
        self.do_line_plot_x_y(ax, ix_pet, pet_mean, 'PET Mean', 'b')
        self.do_line_plot_x_y(ax, ix_pet, pet_min, 'PET Min', 'c')

        plt.savefig(os.path.join(out_dir,
                                 'annual_ppt_pet_cycle_SA.png'),
                    frameon=True, papertype='a4', bbox_inches='tight')
        return

    def plt_ppt_pet_vals(self, in_ppt_df, in_pet_df, out_dir):
        '''Plot of Origibal PPT and PET data along time'''
        if len(in_pet_df.index) < len(in_ppt_df.index):
            cmn_idx = in_ppt_df.index.intersection(in_pet_df.index)
        else:
            cmn_idx = in_pet_df.index.intersection(in_ppt_df.index)
        plt.figure(figsize=(20, 8))
        plt.plot(cmn_idx, in_ppt_df.loc[cmn_idx, :].values, c='b',
                 label='ppt', alpha=0.7)
        plt.plot(cmn_idx, in_pet_df.loc[cmn_idx, :].values, c='r',
                 label='pet', alpha=0.7)
        plt.legend(loc=0)
        plt.grid(alpha=0.5)
        plt.xlabel('Time')
        plt.ylabel('mmm / month')
        plt.savefig(os.path.join(out_dir, 'ppt_pet_all.png'),
                    frameon=True, papertype='a4', bbox_inches='tight')

    def plot_hist(self, data, _label, norme_it, var, time, color, out_dir):
        '''plot hist of a data set'''
        plt.figure(figsize=(20, 8))

        bins = np.linspace(data.min(), data.max(), 50)
        center = (bins[:-1] + bins[1:]) / 2
        hist, _ = np.histogram(data, bins=bins,
                               density=norme_it)

        plt.bar(center, hist, align='center',
                width=.5, alpha=0.7,
                linewidth=2, edgecolor='k',
                color=color,
                label=(_label + ' for ' + time))
        plt.title(_label + ' for ' + time)
        plt.legend(loc=0)
        plt.grid(alpha=0.5)
        plt.xlabel('mmm / month')
        plt.ylabel('Frequency')
        plt.savefig(os.path.join(out_dir, 'hist_%s_%s.png' % (var, time)),
                    frameon=True, papertype='a4', bbox_inches='tight')

    def plot_correl_ppt_enso(self, ppt_data_orig, enso_data, out_dir):
        ''' fct to plot the correlation between ppt data and Oni index data'''
        oni_dict = {'DJF': [12, 1, 2], 'JFM': [1, 2, 3],
                    'FMA': [2, 3, 4], 'MAM': [3, 4, 5],
                    'AMJ': [4, 5, 6], 'MJJ': [5, 6, 7],
                    'JJA': [6, 7, 8], 'JAS': [7, 8, 9],
                    'ASO': [8, 9, 10], 'SON': [9, 10, 11],
                    'OND': [10, 11, 12], 'NDJ': [11, 12, 1]}
        assert isinstance(enso_data, pd.DataFrame)
        for period in oni_dict.keys():
            ppt_data = self.resampleDfbySeason(ppt_data_orig, oni_dict, period)
            ppt_data_avg = ppt_data.groupby(ppt_data.index.year, axis=0).mean()
            enso_data.dropna(inplace=True)
            oni_df = enso_data.loc[ppt_data.index[0]:ppt_data.index[-1],
                                   period]
            oni_df.index = oni_df.index.year

            cmn_idx = ppt_data_avg.index.intersection(oni_df.index)

            ppt_data_season = ppt_data_avg.loc[cmn_idx, :].values
            plt.figure(figsize=(20, 8))
            plt.scatter(cmn_idx, ppt_data_season, marker='o')
            plt.scatter(cmn_idx, (oni_df.values), marker='^')
            plt.ylim([(oni_df.values).min(), ppt_data_season.max()])
            corr = stats.spearmanr(ppt_data_season, oni_df.values)
            plt.title(str(corr) + ' for ' + period)
            plt.savefig(os.path.join(out_dir, 'pet_vs_Oni_%s.png' % (period)),
                        frameon=True, papertype='a4', bbox_inches='tight')

    def plot_scatter_plot_two_values(self, index1, index2, xlabel,
                                     ylabel, season, out_dir):
        ''' fct to plot scatter plot between two indices based on season'''
        vals1 = np.squeeze(index1)
        vals2 = np.squeeze(index2)
        assert len(vals1) == len(vals2), 'wrong shape'
        correl = stats.spearmanr(vals1, vals2)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.scatter(vals1, vals2, s=20, c='b', marker='o',
                    label=correl)
        plt.title(str(correl) + 'between %s and %s' % (xlabel, ylabel))
        plt.show()
        plt.savefig(os.path.join('correlation_between_%s_and_%s_for_%s.png'
                                 % (xlabel, ylabel, season), out_dir),
                    frameon=True, papertype='a4', bbox_inches='tight')
