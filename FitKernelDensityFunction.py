# -*- coding: utf-8 -*-
"""
Created on 13-08-2018

@author: EL Hachem Abbas,
Institut f�r Wasser- und Umweltsystemmodellierung - IWS
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as opt
import scipy.stats as stats


class FitKernelDensity(object):
    def __init__(self):
        pass

    def kernel_gauss(self, t, d):
        # Gauss-Kernel
        return 1 / (d * np.sqrt(2 * np.pi)) * np.exp(-t**2. / (2. * d**2.))

    def calc_edf(self, data, threshold=None, decimals=3, plot=False):
        if threshold is not None:
            data[data <= threshold] = 0
        rank_data = stats.rankdata(np.round(data, decimals), method='max')
        edf = np.unique(rank_data) / np.float(data.shape[0])
        unique = np.unique(np.round(data, decimals))

        if edf.shape[0] != unique.shape[0]:
            raise Exception('EDF is not unique!')

        edf, xval = np.array((edf, unique))
        if plot == True:
            plt.plot(xval, edf, 'x', label='EDF')

        return edf, xval

    def likelyhood_function(self, d, data, this_kernel, opt_quant):
        # get the indices of the quantiles
        q_low = int(data.shape[0] * opt_quant[0])
        q_up = int(data.shape[0] * opt_quant[1])
        # Sort the values
        data_s = np.sort(np.copy(data))
        # Take the sample for optimization depending on the given quantiles
        data_s = data_s[q_low:q_up]
        f_x = np.zeros((data_s.shape[0]))
        f_x[:] = np.nan
        for ii, idata in enumerate(data_s):
            # data without the data equal to the ith xtick
            # get the indices
            idx_del = np.where((data_s == idata))[0]
            data_w = np.delete(data_s, idx_del)
            n_data_w = np.float(data_w.shape[0])
            # calculate the probability density
            if this_kernel == 'gauss':
                # * np.exp(idata) / np.sum(np.exp(data)))
                f_x[ii] = 1 / n_data_w * \
                    np.sum((self.kernel_gauss((idata - data_w), d)))
        # Calculate neg. logaritmic likelyhood-function
        # Which is to be minimized
        return np.sum(-np.log(f_x))

    def kernel_width_optimization(self, data, this_kernel, opt_quant=[0, 1.],
                                  log=False, verbose=False):
        """ Optimization of the kernel width

        Parameters
        ----------
        data : 1D array
            Array, which contains the data.
        this_kernel: string
            defines the kernel (gauss, beta_0_1, beta_1_1, beta1)
        opt_quant: list of float [a,b], optional
            defines the interval of the quantile of the data for
             which the kernel
            width is optimized. Default is [0,1], which means all
            values are
            considered.
        log: bool, optional
            if log is true, the kde is performed in the log-space,
            default = False
        Kernel Info
        -----------
        gauss: Gauss-Kernel
        """

        if log:
            data = np.log(data)
            data = data[data > 0]

        output = opt.minimize_scalar(self.likelyhood_function,
                                     args=(data, this_kernel, opt_quant))
        d_opt = output['x']

        if verbose:
            print('Optimized width: ', d_opt)
        return d_opt

    def max_likelyhood_cdf(self, data, d, this_kernel, plot=False, steps=10000,
                           log=False, limits=None, verbose=False):
        """ CDF of a kernel density estimation

        Parameters
        ----------
        data : 1D array
            Array, which contains the data.
        this_kernel: string
            defines the kernel (gauss, beta_0_1, beta_1_1, beta1)
        plot: bool, optional
            if true, the cdf is plotted in the general figure object,
            default = False
        steps: float, optional
            steps of the kde, default is 10000
        log: bool, optional
            if log is true, the kde is performed in the log-space,
            default = False
        """

        if log:
            data = data[data > 0]  # added by abbas
            data = np.log(data)

        if not limits:
            if this_kernel == 'gauss':
                # The min and max boundary of the ticks is the lowest/highest
                # value +/- 1 times the kernel width
                xmin = data.min() - 1 * d
                xmax = data.max() + 1 * d
        else:
            xmin, xmax = limits
            if log:
                xmin = np.log(xmin)
                xmax = np.log(xmax)

        xticks = np.linspace(xmin, xmax, steps)
        # Maximum Likely Hood Estimation
        n_data = np.float(data.shape[0])
        f_x = np.zeros(xticks.shape[0])
        f_x[:] = np.nan
        dx = xticks[1] - xticks[0]
        for ii, ixtick in enumerate(xticks):
            if this_kernel == 'gauss':
                f_x[ii] = 1 / n_data * \
                    np.sum(self.kernel_gauss(ixtick - data, d))

        F_x = np.cumsum(f_x) * dx

        if verbose:
            print('{}: {:6.5f}, {:6.5f}'.format(this_kernel, F_x[0], F_x[-1]))
        # The CDF will not be exactly between 0 and 1, thus
        # stretch the CDF in order to be between [0,1]
        F_x = (F_x[:] - F_x[0]) / (F_x[:] - F_x[0])[-1]

        if log:
            xticks = np.exp(xticks)
        if plot == True:
            plt.plot(xticks, F_x, label='KDE_{}'.format(this_kernel))
        if xmax > 1:
            "Warning, x values are larger than 1"
        if xmin < -1:
            "Warning, x values are smaller than -1"
        return xticks, F_x


if __name__ == '__main__':
    main_dir = os.path.join(
        r'X:\hiwi\ElHachem\AdvancedPython\KernelDensityFunction')
    os.chdir(main_dir)
    # change file's name ### ATTENTION SEPARATOR INDEX!
    ppt_data = pd.read_csv('pet_data.csv', index_col=0, sep=';')
	
    data = np.sort(ppt_data.values)

    FitKernel = FitKernelDensity()
    # Empirical Distribution
    edf, xval = FitKernel.calc_edf(data, plot=True)

    for ikernel in ['gauss']:
        #
        # Find the optimial Kernelwidth
        d_opt = FitKernel.kernel_width_optimization(data, ikernel, log=False)

        # Fitted maximum likely hood cdf
        #d_opt = 0.1
        # In the line en bas, log=False to consider all the data event the
        # zeros in the data /log= True to take into consideration just the data
        # > 0
        F_x = FitKernel.max_likelyhood_cdf(data, d_opt, ikernel,
                                           plot=True, log=False)

    plt.legend(loc='lower right')
    plt.grid()

    plt.savefig('KDE_gauss_pet_africa.pdf')
    print('done with everything')
