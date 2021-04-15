from PIL import Image
import numpy as np
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


if __name__ == '__main__':
    from ali_leed_import import load_image, show_image

    fn = r'\14_10_20_C60_Au111_1654.tiff'

    data = load_image(fn)
    fig = show_image(data)
    data[:][0][0]  # z=[2,4,0], y=0, x=0
    np.array(fig.data[0]).shape
    fig.data
    # fig.layout.hovermode = 'closest'  # default
    # print(fig.layout.hoverlabel)
    # fig.layout.hoverlabel.
    # scatter = fig.data
    # scatter.on_click()
    # print(fig)
    # fig.data.
    # fig.add_annotation(print(hovertext))
    # print(fig.layout.hover)

    from plotly.callbacks import Points, InputDeviceState
    points, state = Points(), InputDeviceState()
    f = go.FigureWidget(fig)
    f.on_click()

    test = go.FigureWidget([go.Scatter()])
    help(test.data[0].on_click)

    def click_fn(trace, points, state):
        inds = points.point_inds
        print(inds)
        print(scatter.selectedpoints)
        return inds

    fig2=go.FigureWidget(go.Scatter(x=[1, 2], y=[3, 0]))
    # fig3=fig.add_trace()
    scatter = fig2.data[0]
    scatter.on_click(click_fn)
    fig2.show(renderer=DEFAULT_RENDERER)