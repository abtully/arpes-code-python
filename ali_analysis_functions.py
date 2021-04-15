from ali_functions import get_data_index
import numpy as np


# set region of data
def get_data_region(data: np.ndarray, xaxis: np.ndarray, yaxis: np.ndarray, xbounds: tuple = None, ybounds: tuple = None
                    ) -> tuple:
    if xbounds is not None:
        data_x = data[:, np.logical_and(xbounds[0] < xaxis, xaxis < xbounds[1])]
        new_xaxis = xaxis[np.logical_and(xbounds[0] < xaxis, xaxis < xbounds[1])]
    elif xbounds is None:
        data_x = data
        new_xaxis = xaxis
    if ybounds is not None:
        new_region = data_x[np.logical_and(ybounds[0] < yaxis - 16.8, yaxis - 16.8 < ybounds[1])]
        new_yaxis = yaxis[np.logical_and(ybounds[0] < yaxis - 16.8, yaxis - 16.8 < ybounds[1])]
    elif ybounds is None:
        new_region = data_x
        new_yaxis = yaxis
    return new_region, new_xaxis, new_yaxis


def gaussian(x, amplitude, center, width, const):
    return amplitude * np.exp(-(x - center) ** 2 / width) + const


def get_vertical_slice(data, axis, value, interval):
    """
    Gets vertical chunk of data

    Args:
        data:
        axis:
        value:
        interval:

    Returns:
        1D or 2D vertical chunk of data (kx)

    """
    high = value + interval / 2
    low = value - interval / 2
    low_index, high_index = get_data_index(axis, (low, high))
    if low_index == high_index:
        chunk = data[:, low_index]
    else:
        chunk = data[:, low_index:high_index]
    return chunk


def get_averaged_slice(data) -> np.ndarray:
    """
    Averages a slice of data
    Args:
        data:

    Returns:
        1D numpy array

    """
    data = np.atleast_2d(data)
    return np.mean(data, axis=1)


if __name__ == '__main__':
    """Fit Gaussian or Lorentzian to data"""
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    from scipy.stats import norm
    #
    # n = len(c60_homox_b)
    # mean = np.mean()
    # sigma = sum(c60_homoy_b * (c60_homox_b - mean) ** 2) / n  # note this correction
    #
    #
    # def gaus(x, a, x0, sigma):
    #     return a * exp(-(x - x0) ** 2 / (2 * sigma ** 2))
    #
    # curve_fit(gaus, c60_homo_b[0], c60_homo_b[1])
