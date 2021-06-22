"""
Plotting functions for angle-resolved photoemission spectroscopy (ARPES) data
@author: Alexandra Tully
@date: November 2020

Notes:
Colorscale options: https://plotly.com/python/builtin-colorscales/
"""

import plotly.graph_objects as go
from arpes_dataclasses import Data2D, Data3D
from typing import List, Optional
import numpy as np
import copy

DEFAULT_RENDERER = 'browser'  # this is a constant of this file

"""Plotting Functions 2D"""


def plot2D(x: np.ndarray, y: np.ndarray, data: np.ndarray, show=True, opacity=None, xlabel=None, ylabel=None,
           title=None, xrange: tuple = None, yrange: tuple = None, dark_mode=False, colorscale='Plasma',
           E_f: float = None):
    fig = go.Figure()
    if E_f is not None and E_f != '':
        y = y-float(E_f)
        ylabel='E - E_f'
    hm = heatmap(x, y, data, opacity=opacity)
    fig.add_trace(hm)
    fig.update_layout(dict(xaxis_title=xlabel, yaxis_title=ylabel, title=title), coloraxis=dict(colorscale=colorscale))
    if dark_mode:
        fig.update_layout(template='plotly_dark')  # dark mode
    if xrange:
        fig.update_layout(xaxis=dict(range=[xrange[0], xrange[1]]))
    if yrange:
        fig.update_layout(yaxis=dict(range=[yrange[0], yrange[1]]))
    if show:
        fig.show(renderer=DEFAULT_RENDERER)
    return fig


def multi_plot2D(heatmaps: List[go.Heatmap], xlabel=None, ylabel=None,
                 title=None, xrange: tuple = None, yrange: tuple = None,
                 dark_mode=False, colorscale=None):  # expecting heatmaps to be a list of go.Heatmaps (class)
    fig = go.Figure()
    for hm in heatmaps:
        fig.add_trace(hm)
    fig.update_layout(xaxis_title=xlabel, yaxis_title=ylabel, title=title, coloraxis=dict(colorscale=colorscale))
    if dark_mode:
        fig.update_layout(template='plotly_dark')  # dark mode
    if xrange:
        fig.update_layout(xaxis=dict(range=[xrange[0], xrange[1]]))
    if yrange:
        fig.update_layout(yaxis=dict(range=[yrange[0], yrange[1]]))
    fig.show(renderer=DEFAULT_RENDERER)
    return fig


def filepath_plot2D(filepath: str, month, year=2020, laser='lamp', cryo_temp='RT', scan_type='HS', k_cut=None,
                    scan_number=1, filename=None, colorscale='Plasma') -> go.Figure:
    data = Data2D.single_load(month, year=year, light_source=laser, cryo_temp=cryo_temp, scan_type=scan_type, k_cut=k_cut,
                              scan_number=scan_number, filename=filename, filepath=filepath)
    return plot2D(data.xaxis, data.yaxis, data.data, show=False, colorscale=colorscale)


"""Generate Heatmaps"""


def heatmap(x, y, data, opacity=None) -> go.Heatmap:  # this function returns a heatmap
    return go.Heatmap(x=x, y=y, z=data, opacity=opacity, coloraxis="coloraxis")


def multi_heatmaps(datas: List[Data2D], E_f: float, opacity: Optional[float] = None) -> List[go.Heatmap]:
    heatmaps = []
    for data in datas:
        hm = heatmap(data.yaxis, data.xaxis - E_f, data.data, opacity=opacity)
        heatmaps.append(hm)
    return heatmaps


"""Plotting Functions 3D"""


# plot 2D slice of 3D data
def plot3D(x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray, slice_dim: str, slice_val: float,
           int_range: float = 0, show=True, opacity=None, xlabel=None, ylabel=None, title=None, xrange: tuple = None,
           yrange: tuple = None, dark_mode=False, colorscale='Plasma'):
    if slice_dim == 'x':
        a = x
        axis_from = 2
        axes_2d = (y, z)
    elif slice_dim == 'y':
        a = y
        axis_from = 1
        axes_2d = (x, z)
    elif slice_dim == 'z':
        a = z
        axis_from = 0
        axes_2d = (x, y)
    else:
        raise ValueError(f'{slice_dim} is not x, y, or z')
    slice_index = np.argmin(np.abs(a - slice_val))  # gives n; z[n] = value closest to slice_val
    int_index_range = np.floor(int_range / (2 * np.mean(np.diff(a)))).astype(
        int)  # gets avg delta between z, rounds down
    data = np.moveaxis(data, axis_from, 0)
    if int_index_range > 0:
        low = slice_index - int_index_range
        low = low if low > 0 else None
        high = slice_index + int_index_range
        high = high if high < data.shape[0] else None
        data2d = data[low: high].mean(axis=0)
    else:
        data2d = data[slice_index]
    return plot2D(axes_2d[0], axes_2d[1], data2d, show=show, opacity=opacity, xlabel=xlabel, ylabel=ylabel, title=title,
                  xrange=xrange, yrange=yrange, dark_mode=dark_mode, colorscale=colorscale)


def filepath_plot3D(filepath: str, month, year=2020, laser='lamp', cryo_temp='RT', scan_type='FS',
                    scan_number=1, slice_dim='z', slice_val=15, int_range=0.02, filename=None, colorscale='Plasma') \
        -> go.Figure:
    # print(type(int_range), type(slice_val), type(slice_dim))
    data = Data3D.single_load(month, year=year, light_source=laser, cryo_temp=cryo_temp, scan_type=scan_type,
                              scan_number=scan_number, filename=filename, filepath=filepath)
    return plot3D(x=data.yaxis, y=data.xaxis, z=data.zaxis, data=np.moveaxis(data.data, 2, 1), slice_dim=slice_dim,
                  slice_val=slice_val, int_range=int_range, show=False, colorscale=colorscale)


def mutliplot3D(x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray, slice_dim: str, slice_val: float,
           int_range: float = 0, show=True, opacity=None, xlabel=None, ylabel=None, title=None, xrange: tuple = None,
           yrange: tuple = None, dark_mode=False, colorscale='Plasma'):
    raise NotImplementedError
    # hm = heatmap3D(x=x, y=y, z=z, data=data, slice_dim=slice_dim, slice_val=slice_val, int_range=int_range, show=show,
    #                opacity=opacity)


def transpose_figure(fig: go.Figure) -> go.Figure:
    fig = copy.copy(fig)
    fig_data = fig.data[0]
    x, y, data = fig_data['x'], fig_data['y'], fig_data['z']
    fig_data['x'], fig_data['y'], fig_data['z'] = y, x, data.T
    fig.data = [fig_data]
    return fig


if __name__ == '__main__':
    from arpes_dataclasses import Data3D, Data2D
    from misc_functions import multi_load

    # d = Data3D.single_load('October', cryo_temp='RT', scan_subtype='Calibrate')
    # plot3D(x=d.yaxis, y=d.xaxis, z=d.zaxis, data=np.moveaxis(d.data, 2, 1), slice_dim='z', slice_val=15,
    # int_range=0.02)
    # d.berend_data.show(mode='CE', val=15, Eint=0.02)

    d = Data2D.single_load('December', year='2020', filename='UPS_20K0001_001.ibw')
    plot2D(d.xaxis, d.yaxis-16.8, d.data)
    d = Data2D.single_load(month='January', year='2021', light_source='XUV', scan_number=4)
    plot2D(d.xaxis, d.yaxis, d.data)

    raster = Data2D.single_load('October', year='2020', cryo_temp='RT', scan_number=15, filename='Raster0015_015.ibw')
    E_f = 16.8
    plot2D(raster.yaxis, raster.xaxis - E_f, raster.data, opacity=0.7, xrange=(-20, 17))

    rasters = multi_load(i for i in range(10, 14))
    heatmaps = multi_heatmaps(rasters, E_f=16.8, opacity=0.5)
    multi_plot2D(heatmaps, title='Test')
