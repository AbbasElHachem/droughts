# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: EL Hachem Abbas,
Institut für Wasser- und Umweltsystemmodellierung - IWS
"""
import os
import time
import timeit
import warnings

import matplotlib.pyplot as plt
# import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy.optimize as optimize


print('\a\a\a\a Started on %s \a\a\a\a\n' % time.asctime())
START = timeit.default_timer()  # to get the runtime of the program

main_dir = os.path.join(
    r'X:\hiwi\ElHachem\AdvancedPython\KernelDensityFunction')
os.chdir(main_dir)


class KernelDensityEstimate(object):

    def __init__(self, *args, **kwargs):
        object.__init__(self, *args, **kwargs)

    def gauss_kernel(self, t, d):
        return (1. / (d * np.sqrt(2. * np.pi))) * \
            np.exp((-t**2.) / (2. * d**2.))

    def epanechnikov_kernel(self, t, d):
        if np.abs(t) <= np.sqrt(5) * d:
            return (3 / 4 * d * np.sqrt(5)) * (1 - (t**2) / (5 * d**2))
        else:
            return 0

    def fill_mtx(self, kernel_width, data):
        out_matx = np.empty((data.shape[0], data.shape[0]))
        for i, punkte in enumerate(data):
            for j, daten_wert in enumerate(data):
                if punkte == daten_wert:
                    out_matx[i, j] = np.nan
                else:
                    out_matx[i, j] = self.gauss_kernel(
                        punkte - daten_wert, kernel_width)
        return out_matx

    def leaveOneOut_likelihood(self, kernel_width, data):
        out_mtx = self.fill_mtx(kernel_width, data)
        return np.sum([-np.log(np.nanmean(out_mtx[r]))
                       for r in range(out_mtx.shape[0])
                       if np.nanmean(out_mtx[r]) != np.empty])

    def optimize_kernel_width(self, data):
        optimal_width = optimize.minimize_scalar(self.leaveOneOut_likelihood,
                                                 args=(data))

        return optimal_width

    def fit_kernel_to_data(self, _data):
        data_points = np.linspace(_data.min(), _data.max(),
                                  len(_data), endpoint=True, dtype=np.float64)

        out_mtx_cal = np.empty((_data.shape[0], _data.shape[0]))
        kernel_width = self.optimize_kernel_width(_data)['x']
        for i, p in enumerate(data_points):
            for j, v in enumerate(_data):
                out_mtx_cal[i, j] = self.gauss_kernel(
                    (p - v), kernel_width)
        norm_vals = [np.mean(out_mtx_cal[r])
                     for r in range(out_mtx_cal.shape[0])]
        return kernel_width, norm_vals, data_points

    # def multiprocess(processes, samples):
    #     with warnings.catch_warnings():
    #         pool = mp.Pool(processes=processes)
    #         results = [pool.apply_async(fit_kernel_to_data, args=(samples,))]
    #         warnings.simplefilter("ignore", category=RuntimeWarning)
    #         results = [p.get() for p in results]
    #         return results


if __name__ == '__main__':
    in_df = os.path.join(main_dir, 'pet_data.csv')
    in_ppt_dt = pd.read_csv(in_df, index_col=0, sep=';')
    # pool = mp.Pool(processes=4)

    # results = [pool.apply_async(fit_kernel_to_data,
    #                             args=(in_ppt_dt.values,))]
    fit_kernel_fct = KernelDensityEstimate()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        d_opt, no, pp = fit_kernel_fct.fit_kernel_to_data(in_ppt_dt.values)

    plt.ioff()
    plt.figure(figsize=(12, 8), dpi=300)
    plt.hist(
        in_ppt_dt.values,
        bins=68,
        density=True,
        alpha=0.75,
        color='b',
        label='monthly ppt data')
    plt.plot(pp, no, color='r', alpha=0.85,
             label='fitted Gauss kernel, d=%0.2f' % d_opt)
    plt.grid(alpha=0.5)
    plt.legend(loc=0)
    plt.xlabel('ppt')
    plt.ylabel('density')
    plt.savefig(os.path.join(main_dir, 'fitted_gauss_kernel_pet_data.png'),
                frameon=True, papertype='a4', bbox_inches='tight')

    STOP = timeit.default_timer()  # Ending time
    print(('\n\a\a\a Done with everything on %s. Total run time was'
           ' about %0.4f seconds \a\a\a' % (time.asctime(), STOP - START)))
