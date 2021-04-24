from PIL import Image
import numpy as np
from scipy.constants import h, m_e
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from typing import List, Optional

DEFAULT_RENDERER = 'browser'  # this is a constant of this file

PATH = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\October_2020\LEED\Spinview'


def intensity(array):
    def pixel_intensity(pixel):  # define pixel intensity subfunction
        return np.sqrt(pixel[0]**2 + pixel[1]**2 + pixel[2]**2)  # adding RBG values in quadrature
    array = array.astype(np.float32)  # change dtype from 8 byte integer to 32 byte float
    if array.ndim == 1 and len(array) == 4:  # if array is a pixel, utilize pixel_intensity subfunction
        return pixel_intensity(array)
    elif array.shape[-1] == 4:  # [-1] is the last value
        sq_array = array[:, :, 0:3]**2  # all x, all y, and 0-3 not inclusive, squared
        sum_sq = sq_array[:, :, 0] + sq_array[:, :, 1] + sq_array[:, :, 2]
        return np.sqrt(sum_sq)


def get_z_data(array, x, y):
    return array[:][y][x]


def get_line_distance(coord1: tuple, coord2: tuple) -> float:
    x1, y1 = coord1[0], coord1[1]
    x2, y2 = coord2[0], coord2[1]
    x = np.abs(x1 - x2)
    y = np.abs(y1 - y2)
    return np.sqrt(x**2 + y**2)


def atomic_row_spacing(x: float, R: float, E: int, n=1):  # E in eV, R and x in same units
    return (n*R*12.3)/(x*np.sqrt(E))  # angstroms


def distance_to_screen(a: float, x: float, E: int, n=1):  # E in eV, a and x n same units
    return (a*x*np.sqrt(E))/(12.3*n)  # in angstroms


if __name__ == '__main__':
    from ali_leed_import import load_image, show_image

    fn = r'14_10_20_C60_Au111_1654.tiff'

    d_c60 = get_line_distance((699, 398), (534, 358))  # n=1, 17eV
    d_au = get_line_distance((652, 271), (385, 339))  # n=1, 90eV

    r = distance_to_screen(a=2.88, x=d_au, E=90, n=1)
    atomic_row_spacing(x=d_c60, R=r, E=17, n=1)  # a=10.75 (vs 10.07)

    r = distance_to_screen(a=10, x=d_c60, E=17, n=1)
    atomic_row_spacing(x=d_au, R=r, E=90, n=1)  # a=2.68 (vs 2.88)


    # data = load_image(fn, path=PATH)
    # fig = show_image(data)
    # data[:][0][0]  # z=[2,4,0], y=0, x=0
    # np.array(fig.data[0]).shape
    # fig.data

    # fig.layout.hovermode = 'closest'  # default
    # print(fig.layout.hoverlabel)
    # fig.layout.hoverlabel.
    # scatter = fig.data
    # scatter.on_click()
    # print(fig)
    # fig.data.
    # fig.add_annotation(print(hovertext))
    # print(fig.layout.hover)
    #
    # from plotly.callbacks import Points, InputDeviceState
    # points, state = Points(), InputDeviceState()
    # f = go.FigureWidget(fig)
    # f.on_click()
    #
    # test = go.FigureWidget([go.Scatter()])
    # help(test.data[0].on_click)
    #
    # def click_fn(trace, points, state):
    #     inds = points.point_inds
    #     print(inds)
    #     print(scatter.selectedpoints)
    #     return inds
    #
    # fig2=go.FigureWidget(go.Scatter(x=[1, 2], y=[3, 0]))
    # # fig3=fig.add_trace()
    # scatter = fig2.data[0]
    # scatter.on_click(click_fn)
    # fig2.show(renderer=DEFAULT_RENDERER)