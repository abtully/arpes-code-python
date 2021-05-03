"""
Plotting Functions for ARPES Data
@author: Alexandra Tully
@date: November 2020

Notes:
Colorscale options: https://plotly.com/python/builtin-colorscales/
"""

import plotly.graph_objects as go
from ali_classes import Data2D, Data3D
from typing import List, Optional
import numpy as np

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


def filepath_plot3D(filepath: str, month, year=2020, laser='lamp', cryo_temp='RT', scan_type='FS', scan_subtype=None,
                    scan_number=1, slice_dim='z', slice_val=15, int_range=0.02, filename=None, colorscale='Plasma') -> go.Figure:
    # print(type(int_range), type(slice_val), type(slice_dim))
    data = Data3D.single_load(month, year=year, light_source=laser, cryo_temp=cryo_temp, scan_type=scan_type,
                              scan_subtype=scan_subtype, scan_number=scan_number, filename=filename, filepath=filepath)
    # print('Here')
    return plot3D(x=data.yaxis, y=data.xaxis, z=data.zaxis, data=np.moveaxis(data.data, 2, 1), slice_dim=slice_dim,
                  slice_val=slice_val, int_range=int_range, show=False, colorscale=colorscale)


def mutliplot3D(x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray, slice_dim: str, slice_val: float,
           int_range: float = 0, show=True, opacity=None, xlabel=None, ylabel=None, title=None, xrange: tuple = None,
           yrange: tuple = None, dark_mode=False, colorscale='Plasma'):

    hm = heatmap3D(x=x, y=y, z=z, data=data, slice_dim=slice_dim, slice_val=slice_val, int_range=int_range, show=show,
                   opacity=opacity)


if __name__ == '__main__':
    from ali_classes import Data3D, Data2D

    # d = Data3D.single_load('October', cryo_temp='RT', scan_subtype='Calibrate')
    # plot3D(x=d.yaxis, y=d.xaxis, z=d.zaxis, data=np.moveaxis(d.data, 2, 1), slice_dim='z', slice_val=15,
    # int_range=0.02)
    # d.berend_data.show(mode='CE', val=15, Eint=0.02)

    d = Data2D.single_load('December', year='2020', light_source='XUV', filename='UPS_20K0001_001.ibw')
    d = Data2D.single_load('December', year='2020', filename='UPS_20K0001_001.ibw')
    plot2D(d.xaxis, d.yaxis-16.8, d.data)
    d = Data2D.single_load(month='January', year='2021', light_source='XUV', scan_number=4)
    plot2D(d.xaxis, d.yaxis, d.data)

    # raster = Data2D.single_load('October', scan_type='Raster', scan_number=15)
    # E_f = 16.8
    # plot2D(raster.yaxis, raster.xaxis - E_f, raster.data, opacity=0.7, xrange=(-20, 17))
    #
    # from AliFunctions import multi_load
    # rasters = multi_load(i for i in range(10, 14))
    # heatmaps = multi_heatmaps(rasters, E_f=16.8, opacity=0.5)
    # multi_plot2D(heatmaps, title='Test')

    # filepath_plot2D('C:/Users/atully/Code/ARPES Code Python/data/', 'October', k_cut='GM')
    # filepath_plot3D('C:/Users/atully/Code/ARPES Code Python/data/', month='October', year=2020, scan_type='FS',
    #                 scan_subtype='Fermi_Surface', scan_number=1, slice_dim='z', slice_val=15, int_range=0.02,
    #                 colorscale='Plasma', cryo_temp='140K')

    # [i for i in range(10, 20)]  # list of numbers 10 to 19
    # np.arange(10, 20)  # an array from 10 to 20

    # data = np.random.random((3, 3))
    # data2 = np.random.random((3, 3)) * 2
    # x = np.linspace(0, 1, 3)
    # fig = go.Figure()
    # fig.add_trace(go.Heatmap(x=x, y=x, z=data, coloraxis="coloraxis"))
    # fig.add_trace(go.Heatmap(x=x, y=x, z=data2, coloraxis="coloraxis"))
    # fig.update_layout(xaxis=dict(range=[0.4, 0.8]))
    # fig.show(renderer=DEFAULT_RENDERER)

    # E_f = 16.8
    # HS = Data2D.single_load('October', k_cut='GM')
    # plot2D(HS.yaxis, HS.xaxis - E_f, HS.data, opacity=0.5, xlabel='Theta', ylabel='E-E_F')
    #
    # E_f = 16.8
    # Raster = Data2D.single_load('October', scan_type='Raster', scan_number=12)
    # _, Raster12 = plot2D(Raster.yaxis, Raster.xaxis - E_f, Raster.data, opacity=0.5, xlabel='Theta', ylabel='E-E_F')
    #
    # Raster = Data2D.single_load('October', scan_type='Raster', scan_number=15)
    # _, Raster15 = plot2D(Raster.yaxis, Raster.xaxis - E_f, Raster.data, opacity=0.5, xlabel='Theta', ylabel='E-E_F')
    #
    # multi_plot2D(Raster12, Raster15)
    #
    # Raster = [Data2D.single_load('October', scan_type='Raster', scan_number=i) for i in
    #           np.ndarray.tolist([np.linspace(1, 24, 25)])]
    #
    # [np.linspace(1, 24, 25)]
    # x = np.linspace(0, 100, 101)
    # y = np.linspace(0, 100, 101)
    # xx, yy = np.meshgrid(x, y)
    # data = np.sin(xx) + np.cos(yy)
    # plot2D(x, y, data)
    #
