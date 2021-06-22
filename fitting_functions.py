"""
General fitting functions
@author: Alexandra Tully
@date: June 2021

Note: Includes linear, lorentzian, and hyperbolic fits, as well an offset function (constant, linear, or quadratic).
"""

import numpy as np
import os
import plotly.graph_objects as go
import lmfit as lm
from typing import Union, Optional, Tuple
from scipy import optimize
from scipy.ndimage import gaussian_filter

from HDF5_loader import data_from_hdf, avg_array_from_hdfs, avg_data_hdf
from plotting_functions import plot3D, plot2D, transpose_figure
from analysis_functions import get_data_region, get_vertical_slice, get_averaged_slice, get_horizontal_slice
from k_correction import kcorrect_phimotor, get2Dslice, kcorrect2D_general, kcorrect2D, fix_EkatEF
from polygons import gen_polygon, gen_tiled_hexagons, plot_polygons, plot_polygon
from arpes_dataclasses import Data2D


def offset_model(offset_type: str, a: float = None, b: float = None, c: float = None) -> lm.models.QuadraticModel:
    """

    Args:
        offset_type: constant, linear, or quadratic
        a: quadratic coefficient (ax^2 + bx + c)
        b: linear coefficient
        c: constant coefficient

    Returns:
        offset model to be added to lmfit model (e.g. model + offset_model)

    """
    model = lm.models.QuadraticModel()

    model.set_param_hint('a', value=a)
    model.set_param_hint('b', value=b)
    model.set_param_hint('c', value=c)

    if offset_type == 'quadratic':
        pass
    elif offset_type == 'linear':
        model.set_param_hint('a', value=0, vary=False)  # not allowed to vary parameter 'a' while fitting
    elif offset_type == 'constant':
        model.set_param_hint('a', value=0, vary=False)
        model.set_param_hint('b', value=0, vary=False)
    else:
        raise ValueError(f'offset_type: {offset_type} is not quadratic, linear, or constant')

    return model


def make_line(num, a, b, pos_slope=True) -> lm.models.LinearModel:
    """

    Args:
        num: index of line (0, 1, 2)
        a: initial guess (y = ax + b)
        b: initial guess

    Returns:
        linear model for fit

    """
    pref = f'i{num}_'
    model = lm.models.LinearModel(prefix=pref)
    if pos_slope:
        model.set_param_hint(pref + 'a', value=a, min=0, max=5 * a)
    elif not pos_slope:
        model.set_param_hint(pref + 'a', value=a, min=5 * a, max=0)
    model.set_param_hint(pref + 'b', value=b, min=b - 50, max=b + 50)
    return model


def make_n_lines(num, aes: Union[list, float], bes: Union[list, float]) -> lm.Model:
    """

    Args:
        num: number of lines for fit
        aes: initial guess(es)
        bes: initial guesses(es)

    Returns:
        lmfit model

    """
    if not isinstance(aes, list):
        aes = [aes] * num
    if not isinstance(bes, list):
        bes = [bes] * num

    if any([len(arr) != num for arr in [aes, bes]]):
        raise ValueError(f'length of {aes}, {bes}, not all equal to {num}')

    model = None
    for i, a, b in zip(range(num), aes, bes):
        if a >= 0:
            this_model = make_line(i, a, b, pos_slope=True)
        elif a < 0:
            this_model = make_line(i, a, b, pos_slope=False)
        if not model:
            model = this_model
        else:
            model += this_model
    return model


def fit_linear_data(x: np.ndarray, data: np.ndarray,
                        num: int,
                        aes: Union[list, float] = None, bes: Union[list, float] = None,
                        offset_type: Optional[str] = None, a: float = None, b: float = None, c: float = None) \
        -> lm.model.ModelResult:
    """

    Args:
        x: 1D numpy array xdata
        data: 1D numpy array ydata
        num: number of lines
        amplitudes: initial guess(es)
        centers: initial guesses
        sigmas: initial guess(es)
        offset_type: constant, linear, or quadratic
        a: quadratic coefficient (ax^2 + bx + c)
        b: linear coefficient
        c: constant coefficient

    Returns:
        lmfit model result

    """
    if c is None:
        c = np.mean(data)
    if b is None:
        b = (data[-1] - data[0]) / (x[-1] - x[0])
    if a is None:
        a = 0

    if bes is None:
        bes = np.mean(data)
    if aes is None:
        aes = (data[-1] - data[0]) / (x[-1] - x[0])

    lines = make_n_lines(num, aes, bes)
    if offset_type:
        offset = offset_model(offset_type, a, b, c)
        model = lines + offset
    else:
        model = lines

    fit = model.fit(data.astype(np.float32), x=x.astype(np.float32))

    return fit


def make_lorentzian(num, amplitude, center, sigma) -> lm.models.LorentzianModel:
    """

    Args:
        num: index of lorentzian (0, 1, 2)
        amplitude: initial guess
        center: initial guess
        sigma: initial guess

    Returns:
        lorentzian model for fit

    """
    pref = f'i{num}_'
    model = lm.models.LorentzianModel(prefix=pref)
    model.set_param_hint(pref + 'amplitude', value=amplitude, min=0, max=5 * amplitude)
    model.set_param_hint(pref + 'center', value=center, min=center - 10, max=center + 10)
    model.set_param_hint(pref + 'sigma', value=sigma, min=0, max=10)
    return model


def make_n_lorentzians(num, amplitudes: Union[list, float], centers: list, sigmas: Union[list, float]) -> lm.Model:
    """

    Args:
        num: number of lorentzians for fit
        amplitudes: initial guess(es)
        centers: initial guesses
        sigmas: initial guess(es)

    Returns:
        lmfit model

    """
    if not isinstance(amplitudes, list):
        amplitudes = [amplitudes] * num
    if not isinstance(sigmas, list):
        sigmas = [sigmas] * num

    if any([len(arr) != num for arr in [amplitudes, centers, sigmas]]):
        raise ValueError(f'length of {amplitudes}, {centers}, {sigmas} not all equal to {num}')

    model = None
    for i, amp, cent, sig in zip(range(num), amplitudes, centers, sigmas):
        this_model = make_lorentzian(i, amp, cent, sig)
        if not model:
            model = this_model
        else:
            model += this_model
    return model


def fit_lorentzian_data(x: np.ndarray, data: np.ndarray,
                        num_peaks: int,
                        amplitudes: Union[list, float], centers: list, sigmas: Union[list, float] = None,
                        offset_type: Optional[str] = 'linear', a: float = None, b: float = None, c: float = None) \
        -> lm.model.ModelResult:
    """

    Args:
        x: 1D numpy array xdata
        data: 1D numpy array ydata
        num_peaks: number of lorentzian peaks
        amplitudes: initial guess(es) for lorentzian peaks
        centers: initial guesses for lorentzian peaks
        sigmas: initial guess(es) for lorentzian peaks
        offset_type: constant, linear, or quadratic
        a: quadratic coefficient (ax^2 + bx + c)
        b: linear coefficient
        c: constant coefficient

    Returns:
        lmfit model result

    """
    if c is None:
        c = np.mean(data)
    if b is None:
        b = (data[-1] - data[0]) / (x[-1] - x[0])
    if a is None:
        a = 0
    if sigmas is None:
        sigmas = 1

    lorentzians = make_n_lorentzians(num_peaks, amplitudes, centers, sigmas)
    if offset_type:
        offset = offset_model(offset_type, a, b, c)
        model = lorentzians + offset
    else:
        model = lorentzians

    fit = model.fit(data.astype(np.float32), x=x.astype(np.float32))

    return fit


def hyperbola(x, a, b, h, k, pos=False):
    # return -np.sqrt(a**2 * b**2 + a**2 * (x-h)**2) / b + k
    if pos:
        return np.sqrt(a ** 2 * (1 + (x - h) ** 2 / b ** 2)) + k
    else:
        return -np.sqrt(a**2 * (1 + (x-h)**2 / b**2)) + k


def make_hyperbola(num, a, b, h, k) -> lm.Model:
    """

    Args:
        num: index of hyperbola (0, 1, 2)
        a: initial guess
        b: initial guess
        h: initial guess for center (h, k)
        k: initial guess for center (h, k)

    Returns:
        lmfit model of hyperbola

    """
    pref = f'i{num}_'
    model = lm.Model(hyperbola, prefix=pref)
    model.set_param_hint(pref + 'a', value=a)
    model.set_param_hint(pref + 'b', value=b)
    model.set_param_hint(pref + 'h', value=h, min=h - 5, max=h + 5)
    model.set_param_hint(pref + 'k', value=k, min=k-4, max=k+4)
    return model


def make_n_hyperbolas(num,
                      aes: Union[list, float], bes: Union[list, float], hes: Union[list, float], kes: Union[list, float]) \
        -> lm.Model:
    """

    Args:
        num: number of hyperbolas for fit
        aes: initial guess(es)
        bes: initial guess(es)
        hes: initial guess(es) for center (h, k)
        kes: initial guess(es) for center (h, k)

    Returns:
        lmfit model

    """
    if not isinstance(aes, list):
        aes = [aes] * num
    if not isinstance(bes, list):
        bes = [bes] * num
    if not isinstance(hes, list):
        hes = [hes] * num
    if not isinstance(kes, list):
        kes = [kes] * num

    if any([len(arr) != num for arr in [aes, bes, hes, kes]]):
        raise ValueError(f'length of {aes}, {bes}, {hes}, {kes} not all equal to {num}')

    model = None
    for i, a, b, h, k in zip(range(num), aes, bes, hes, kes):
        this_model = make_hyperbola(i, a, b, h, k)
        if not model:
            model = this_model
        else:
            model += this_model
    return model


def fit_hyperbola_data(x: np.ndarray, data: np.ndarray,
                        num: int,
                        aes: Union[list, float], bes: Union[list, float], hes: Union[list, float], kes: Union[list, float],
                        offset_type: Optional[str] = 'linear', a: float = None, b: float = None, c: float = None) \
        -> lm.model.ModelResult:
    """

    Args:
        x: 1D numpy array xdata
        data: 1D numpy array ydata
        num: number of hyperbolas to fit
        aes: initial guess(es)
        bes: initial guess(es)
        hes: initial guess(es) for center (h, k)
        kes: initial guess(es) for center (h, k)
        offset_type: constant, linear, or quadratic
        a: quadratic coefficient (ax^2 + bx + c)
        b: linear coefficient
        c: constant coefficient

    Returns:
        lmfit model result

    """
    if c is None:
        c = np.mean(data)
    if b is None:
        b = (data[-1] - data[0]) / (x[-1] - x[0])
    if a is None:
        a = 0
    # if sigmas is None:
    #     sigmas = 1

    hyperbolas = make_n_hyperbolas(num, aes, bes, hes, kes)
    if offset_type:
        offset = offset_model(offset_type, a, b, c)
        model = hyperbolas + offset
    else:
        model = hyperbolas

    fit = model.fit(data.astype(np.float32), x=x.astype(np.float32))

    return fit


def make_hyperbola_asymptotes(fit, x):
    """

    Args:
        fit: lmfit model fit
        x: 1D numpy array xdata

    Returns:
        y1, y2 lines for hyperbola

    """
    a, b, h, k = fit.best_values['i0_a'], fit.best_values['i0_b'], fit.best_values['i0_h'], fit.best_values['i0_k']
    y1 = a / b * (x - h) + k
    y2 = - a / b * (x - h) + k
    return y1, y2


if __name__ == '__main__':

    data = Data2D.single_load(month='January', year='2021', light_source='XUV', scan_number=4)

    ## take gaussian mask and differentiate data

    # gauss_data = gaussian_filter(d, sigma=5)  # change sigma for different plots
    # diff_data = np.diff(gauss_data, n=2, axis=0)  # this procedure removes 2 indices from whichever axis you're
    #
    # # differentiating over, in this case yaxis.
    # new_yaxis = np.linspace(y[0], y[-1], diff_data.shape[0])  # fix yaxis
    #
    # # adjust region of data and plot
    # nd, nx, ny = get_data_region(d, xaxis=x, yaxis=y, xbounds=(-10, 10))
    # gdd, gdx, gdy = get_data_region(gauss_data, xaxis=x, yaxis=y, xbounds=(-10, 10))
    # dd, dx, dy = get_data_region(diff_data, xaxis=x, yaxis=new_yaxis, xbounds=(-0.9, 0.8))
    #
    # # Final Plots
    # plot2D(nx, ny - EF, nd, title=f'Data', xlabel='', ylabel='E_B')
    # plot2D(gdx, gdy - EF, gdd, title=f'Gaussian Mask (sigma = 5)', xlabel='', ylabel='E_B')
    # fig = plot2D(dx, dy - EF, d, title=f'2nd Derivative of (Masked) Data', xlabel='', ylabel='E_B')
    # fig.update_layout(coloraxis=dict(cmin=-1, cmax=1))
    # fig.show()  # clipped data, everything > 1 = 1 etc.
    #
    # """Fit Gaussians to Vertical Slices of Data"""
    # # testdata = get_averaged_slice(get_vertical_slice(d, x, x[1], 1), axis='x')
    # testdata2 = get_averaged_slice(get_horizontal_slice(d, y, 18, 0.1), axis='y')
    # fig.add_hline(y=18)
    # fig.add_hline(y=18.05, line_dash='dot')
    # fig.add_hline(y=17.95, line_dash='dot')
    # fig.show()
    # # fig.add_trace(go.Scatter(x=x, y=testdata2))
    # fig2 = go.Figure(data=go.Scatter(x=x, y=testdata2))
    # fig2.update_layout(title='Horizontal Chunk')
    # fig2.show(renderer='browser')
    #
    # data = d
    # y = y - 16.8
    # import lmfit as lm
    # from ali_analysis_functions import gaussian
    #
    # # model = lm.models.GaussianModel()
    # model = lm.Model(gaussian)
    # params = model.make_params()
    # params.add('amplitude', 200, min=np.min(data), max=np.max(data))
    # params.add('center', -1.6, min=np.min(y), max=np.max(y))
    # params.add('width', 1, min=0.001, max=np.max(y) - np.min(y))
    # params.add('const', 150, min=0.5 * np.min(data), max=np.max(data))
    # fit = model.fit(testdata2.astype(np.float32), x=y.astype(np.float32), params=params)
    # print(fit.fit_report())
    # fit.plot()
    # fit.best_values
    # fit.best_values['center']
    #
    # fits = []
    # for row in data.T:
    #     fit = model.fit(row.astype(np.float32), x=y.astype(np.float32), params=params)
    #     fits.append(fit)
    #
    # centers = [fit.best_values['center'] for fit in fits]
    # fig.add_trace(go.Scatter(x=x, y=centers))
    # fig.update_layout(title='Gaussian Mask (sigma = 5) with Fit')
    # fig.show(renderer='browser')
    #
    # # k-correct slice in phi FIXME: Not currently working right :(
    # val = slice_val
    # val = -4
    # kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='XUV_FS_averaged.h5', val=val, Eint=2, EF=18.4, slice_dim='x',
    #                                   phi_m=6.5, theta_m=-1, phi0=-13, theta0=-0.6, alpha0=0)
    # fig = plot2D(x=kx, y=np.flip(ky), data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'{val}eV (XUV)')
    #
    # # k-correction + hexagons
    # slice_val = 17.7
    # # slice_val = 6.5
    # val = np.round(slice_val - EF, 1)
    # kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='XUV_FS_averaged.h5', val=val, Eint=0.5, EF=18.4, slice_dim='y',
    #                                   phi_m=6.5, theta_m=-1, phi0=-13, theta0=-0.6, alpha0=0)
    # fig = plot2D(x=kx, y=ky, data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'{val}eV (XUV)')
    # # fig.add_vline(x=0)
    # # fig.add_hline(y=0)
    # # fig.show()
