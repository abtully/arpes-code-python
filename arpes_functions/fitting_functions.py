"""
General fitting functions
@author: Alexandra Tully
@date: June 2021

Note: Includes linear, lorentzian, and hyperbolic fits, as well an offset function (constant, linear, or quadratic).
"""

import numpy as np
import lmfit as lm
from typing import Union, Optional


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
                        offset_type: Optional[str] = 'linear', a: float = None, b: float = None, c: float = None,
                        method: str = 'leastsq',
                        params = None) \
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
        method: minimization methosd (e.g. leastsq, powell, etc.)

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

    fit = model.fit(data.astype(np.float32), x=x.astype(np.float32), method=method, params=params)

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