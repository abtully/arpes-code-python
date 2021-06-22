"""
Filtering functions
@author: Alexandra Tully
@date: June 2021
"""

import numpy as np
from scipy.ndimage import gaussian_filter

from plotting_functions import plot2D
from arpes_functions.arpes_dataclasses import Data2D
from arpes_functions.analysis_functions import get_data_region


def fourier_2d(data: np.ndarray, xaxis: np.ndarray, yaxis: np.ndarray):
    fft = np.fft.fft2(data)
    lx, ly = np.shape(xaxis)[0], np.shape(yaxis)[0]
    dx = (xaxis[-1] - xaxis[0]) / lx
    dy = (yaxis[-1] - yaxis[0]) / ly
    x_fft = 1 / dx * np.linspace(-0.5, 0.5, lx)
    y_fft = 1 / dy * np.linspace(-0.5, 0.5, ly)
    dat_fft = 2 * np.abs(np.fft.fftshift(fft))
    return dat_fft, x_fft, y_fft


if __name__ == '__main__':

    """ Load Raw Data """

    # d = Data2D.single_load('January', year='2021', light_source='XUV', filename='OMBE_XUV_2D0004_.ibw')
    data = Data2D.single_load(month='January', year='2021', light_source='XUV', scan_number=4)
    plot2D(data.xaxis, data.yaxis, data.data, title='OMBE_XUV_2D0004_.ibw, January 2021', xlabel='Theta', ylabel='KE')

    # zoom in on cone
    d, x, y = get_data_region(data.data, data.xaxis, data.yaxis, xbounds=(-12, 8), ybounds=(17.1, 18.4), EB=False)
    fig = plot2D(x, y, d, title='OMBE_XUV_2D0004_.ibw, January 2021', xlabel='Theta', ylabel='KE')

    """ Apply Filters """

    # gaussian mask
    sigma = 2
    gauss_data = gaussian_filter(d, sigma=sigma)
    fig = plot2D(x, y, gauss_data,
                 title=f'OMBE_XUV_2D0004_.ibw, January 2021, Gaussian mask (sigma: {sigma})',
                 xlabel='Theta', ylabel='KE')