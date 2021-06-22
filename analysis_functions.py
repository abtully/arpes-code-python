"""
General analysis functions
@author: Alexandra Tully
@date: November 2020
"""

from misc_functions import get_data_index
import numpy as np


# set region of data
def get_data_region(data: np.ndarray, xaxis: np.ndarray, yaxis: np.ndarray,
                    xbounds: tuple = None, ybounds: tuple = None,
                    EB=False
                    ) -> tuple:
    if xbounds is not None:
        data_x = data[:, np.logical_and(xbounds[0] < xaxis, xaxis < xbounds[1])]
        new_xaxis = xaxis[np.logical_and(xbounds[0] < xaxis, xaxis < xbounds[1])]
    elif xbounds is None:
        data_x = data
        new_xaxis = xaxis
    if ybounds is not None:
        if EB:
            new_region = data_x[np.logical_and(ybounds[0] < yaxis - 16.8, yaxis - 16.8 < ybounds[1])]
            new_yaxis = yaxis[np.logical_and(ybounds[0] < yaxis - 16.8, yaxis - 16.8 < ybounds[1])]
        if not EB:
            new_region = data_x[np.logical_and(ybounds[0] < yaxis, yaxis < ybounds[1])]
            new_yaxis = yaxis[np.logical_and(ybounds[0] < yaxis, yaxis < ybounds[1])]
    elif ybounds is None:
        if EB:
            raise ValueError(f'{ybounds} is None, do not specify EB')
        elif not EB:
            pass
        new_region = data_x
        new_yaxis = yaxis
    return new_region, new_xaxis, new_yaxis


def gaussian(x, amplitude, center, width, const):
    return amplitude * np.exp(-(x - center) ** 2 / width) + const


def get_vertical_slice(data, axis, value, interval):
    """
    Gets vertical chunk of data

    Args:
        data: numpy array
        axis: typically x
        value: in the axis unit (or can specify x[500] if you want the 500th index)
        interval: total interval in the axis unit, centered on value

    Returns:
        1D or 2D vertical chunk of data

    """
    high = value + interval / 2
    low = value - interval / 2
    low_index, high_index = get_data_index(axis, (low, high))
    if low_index == high_index:
        chunk = data[:, low_index]
    else:
        chunk = data[:, low_index:high_index]
    return chunk


def get_horizontal_slice(data, axis, value, interval):
    """
    Gets horizontal chunk of data

    Args:
        data: numpy array
        axis: typically y
        value: in the axis unit (or can specify y[500] if you want the 500th index)
        interval: total interval in the axis unit, centered on value

    Returns:
        1D or 2D horizontal chunk of data

    """
    high = value + interval / 2
    low = value - interval / 2
    low_index, high_index = get_data_index(axis, (low, high))
    if low_index == high_index:
        chunk = data[low_index, :]
    else:
        chunk = data[low_index:high_index, :]
    return chunk


def get_averaged_slice(data, axis) -> np.ndarray:
    """
    Averages a slice of data
    Args:
        data: numpy array
        axis: specify x or y

    Returns:
        1D numpy array

    """
    data = np.atleast_2d(data)
    if axis == 'x':
        return np.mean(data, axis=1)
    if axis == 'y':
        return np.mean(data, axis=0)