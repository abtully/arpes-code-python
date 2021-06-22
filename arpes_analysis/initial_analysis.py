"""
General analysis of early 2021 data
@author: Alexandra Tully
@date: November 2020
"""

from arpes_functions.arpes_dataclasses import Data2D
from arpes_functions.plotting_functions import plot2D
from arpes_functions.misc_functions import bin_data_with_axes
from arpes_functions.analysis_functions import get_data_region, get_vertical_slice, get_averaged_slice
import numpy as np
import plotly.graph_objects as go
import os
from scipy.ndimage import gaussian_filter

if __name__ == '__main__':
    """"""
    """Gaussian Mask, 2nd Derivative, and Analysis of HOMO (lmfit)"""
    """"""
    # load data
    # d = Data2D.single_load(month='January', year='2021', filename='OMBE_Lamp_2D0001_.ibw')
    path = os.path.abspath('/Users/alexandratully/Desktop/ARPES Data/')
    path = os.path.abspath('C:/Users/atully\Code\ARPES Code Python/analysis_data')
    d = Data2D.single_load(month='April', year='2021', light_source='XUV', filepath=path, filename='OMBE_XUV_2D0003_.ibw')
    # data, x, y = fourier_2d(d.data, d.xaxis, d.yaxis)

    # apply mask using scipy's built in gaussian filter
    gauss_data = gaussian_filter(d.data, sigma=5)  # change sigma for different plots
    diff_data = np.diff(gauss_data, n=2, axis=0)  # this procedure removes 2 indices from whichever axis you're
    # differentiating over, in this case yaxis.
    new_yaxis = np.linspace(d.yaxis[0], d.yaxis[-1], diff_data.shape[0])  # fix yaxis

    # plot2D(d.xaxis, d.yaxis - 16.8, d.data, colorscale='RdBu')
    # plot2D(d.xaxis, d.yaxis - 16.8, gauss_data)
    # plot2D(d.xaxis, new_yaxis - 16.8, diff_data, colorscale='Greys')

    # adjust region of data and plot
    dd, dx, dy = get_data_region(d.data, xaxis=d.xaxis, yaxis=d.yaxis, xbounds=(-16, 16),
                                 ybounds=(-4.75, 0.2))
    gdd, gdx, gdy = get_data_region(gauss_data, xaxis=d.xaxis, yaxis=d.yaxis, xbounds=(-16, 16),
                                    ybounds=(-4.75, 0.2))
    c60, c60x, c60y = get_data_region(diff_data, xaxis=d.xaxis, yaxis=new_yaxis, xbounds=(-16, 16),
                                      ybounds=(-4.75, 0.2))

    # Final Plots
    plot2D(dx, dy - 16.8, dd, title='Original Data', xlabel='Theta', ylabel='E_B')
    plot2D(gdx, gdy - 16.8, gdd, title='Gaussian Mask (sigma = 5)', xlabel='Theta', ylabel='E_B')
    fig = plot2D(c60x, c60y - 16.8, c60, title='2nd Derivative of (Masked) Data', xlabel='Theta', ylabel='E_B')
    fig.update_layout(coloraxis=dict(cmin=-1, cmax=1))
    fig.show()  # clipped data, everything > 1 = 1 etc.

    """Analysis of HOMO"""
    # get region of raw data
    homo, homox, homoy = get_data_region(d.data, d.xaxis, d.yaxis, xbounds=(-18, 16),
                                         ybounds=(-2.1, -1.38))
    plot2D(homox, homoy - 16.8, homo, title='Original Data', xlabel='Theta', ylabel='E_B')

    # get region of gaussian mask data
    gauss_homo, gauss_homox, gauss_homoy = get_data_region(gauss_data, d.xaxis, d.yaxis, xbounds=(-18, 16),
                                                           ybounds=(-2.1, -1.38))
    fig2 = plot2D(gauss_homox, gauss_homoy - 16.8, gauss_homo, title='Gaussian Mask (sigma = 5)', xlabel='Theta',
                  ylabel='E_B')

    # bin raw data
    homo_b, homox_b, homoy_b = bin_data_with_axes(homo, (3, 3), x=homox, y=homoy)
    plot2D(homox_b, homoy_b - 16.8, homo_b, title='Binned Original Data', xlabel='Theta', ylabel='E_B')

    # get region of 2nd derivative data
    diff_homo, diff_homox, diff_homoy = get_data_region(diff_data, xaxis=d.xaxis, yaxis=new_yaxis, xbounds=(-18, 16),
                                                        ybounds=(-2.1, -1.38))
    plot2D(diff_homox, diff_homoy - 16.8, diff_homo, title='2nd Derivative of (Masked) Data', xlabel='Theta',
           ylabel='E_B')

    """Fit Gaussians to Vertical Slices of Data"""
    testdata = get_averaged_slice(get_vertical_slice(gauss_homo, gauss_homox, gauss_homox[1], 1))
    data = gauss_homo
    y = gauss_homoy - 16.8
    import lmfit as lm
    from arpes_functions.analysis_functions import gaussian

    # model = lm.models.GaussianModel()
    model = lm.Model(gaussian)
    params = model.make_params()
    params.add('amplitude', 200, min=np.min(data), max=np.max(data))
    params.add('center', -1.6, min=np.min(y), max=np.max(y))
    params.add('width', 1, min=0.001, max=np.max(y) - np.min(y))
    params.add('const', 150, min=0.5 * np.min(data), max=np.max(data))
    fit = model.fit(testdata.astype(np.float32), x=y.astype(np.float32), params=params)
    print(fit.fit_report())
    fit.plot()
    fit.best_values
    fit.best_values['center']

    fits = []
    for row in data.T:
        fit = model.fit(row.astype(np.float32), x=y.astype(np.float32), params=params)
        fits.append(fit)

    centers = [fit.best_values['center'] for fit in fits]
    fig2.add_trace(go.Scatter(x=gauss_homox, y=centers))
    fig2.update_layout(title='Gaussian Mask (sigma = 5) with Fit')
    fig2.show(renderer='browser')

    """"""
    """BAND SUBTRACTION"""
    """"""

    """ Subtract Au(111) bands from C60 on Au(111)"""
    # load data
    data_c60 = Data2D.single_load(month='January', scan_number=1)
    data_au = Data2D.single_load(month='December', year='2020', filename='lamp_overnight0001lamp_overnight.ibw')

    # plot raw data (binding energy)
    plot2D(data_c60.xaxis, data_c60.yaxis - 16.8, data_c60.data, ylabel='Binding Energy', xlabel='Theta',
           title='Lamp Measurement: C60 on Au(111)')
    plot2D(data_au.xaxis, data_au.yaxis - 16.8, data_au.data, ylabel='Binding Energy', xlabel='Theta',
           title='Lamp Measurement: Au(111)')

    # OPTIONAL: bin data
    data_c60.data, data_c60.xaxis, data_c60.yaxis = bin_data_with_axes(data_c60.data, (3, 3), x=data_c60.xaxis,
                                                                       y=data_c60.yaxis)
    data_au.data, data_au.xaxis, data_au.yaxis = bin_data_with_axes(data_au.data, (3, 3), x=data_au.xaxis,
                                                                    y=data_au.yaxis)

    # flip au data to match expected band intensity
    data_au_flipped = np.flip(data_au.data, 1)
    plot2D(data_au.xaxis, data_au.yaxis - 16.8, data_au_flipped, ylabel='Binding Energy', xlabel='Theta',
           title='Lamp Measurement (Flipped): Au(111)')

    # take subset of C60 data (ensure axes match with Au data subset)
    c60, c60x, c60y = get_data_region(data_c60.data, xaxis=data_c60.xaxis, yaxis=data_c60.yaxis, xbounds=(-16, 16),
                                      ybounds=(-4.5, -1.38))
    plot2D(c60x, c60y - 16.8, c60, ylabel='Binding Energy', xlabel='Theta',
           title='Lamp Measurement: C60 on Au(111) (Binned)')

    # take subset of Au data (ensure axes match with C60 data subset)
    auf, aufx, aufy = get_data_region(data_au_flipped, xaxis=data_au.xaxis, yaxis=data_au.yaxis,
                                      xbounds=(-16, 16), ybounds=(-4.5, -1.38))  # alternate ybounds=(-4.43, -1.33)
    plot2D(aufx, aufy - 16.8, auf, ylabel='Binding Energy', xlabel='Theta',
           title='Lamp Measurement: Au(111) (Flipped and Binned)')

    # rescale intensity
    a = np.max(c60) / np.max(auf)
    auf_rescaled = a * auf

    # subtract Au from C60
    new_bands = c60 - auf_rescaled

    # plot2D(c60x, c60y - 16.8, new_bands, colorscale='RdBu', ylabel='Binding Energy', xlabel='Theta',
    #        title='Au(111) Bands Subtracted from C60 on Au(111) Bands')
    plot2D(c60x, c60y - 16.8, new_bands, ylabel='Binding Energy',
           xlabel='Theta', title='Au(111) (Flipped) Bands Subtracted from C60 on Au(111) Bands (Binned)')

    """Set all negative values to zero"""
    data_new_pos = np.where(new_bands < 0, 0, new_bands)
    plot2D(c60x, c60y - 16.8, data_new_pos, ylabel='Binding Energy', xlabel='Theta',
           title='Subtracted Band Plot; All Negative Intensities Set to Zero (Binned)')

    """Set all intensities < 200 to zero"""
    data_200 = np.where(new_bands < 200, 0, new_bands)
    plot2D(c60x, c60y - 16.8, data_200, ylabel='Binding Energy', xlabel='Theta',
           title='Subtracted Band Plot; All Intensities < 200 Set to Zero (Binned)')

    """Analysis of HOMO"""
    c60_homo, c60_homox, c60_homoy = get_data_region(data_c60.data, data_c60.xaxis, data_c60.yaxis, xbounds=(-18, 16),
                                                     ybounds=(-2.1, -1.38))
    plot2D(c60_homox, c60_homoy - 16.8, c60_homo, ylabel='Binding Energy', xlabel='Theta',
           title='C60 HOMO (Binned)')

    test = get_vertical_slice(c60_homo, c60_homox, -5, 3)
    test = get_averaged_slice(test)
    fig = go.Figure()
    test2 = go.Scatter(x=c60_homoy - 16.8, y=test, mode='lines')
    fig.add_trace(test2)
    fig.show(renderer='browser')

    # OPTIONAL: bin data
    c60_homo_b, c60_homox_b, c60_homoy_b = bin_data_with_axes(c60_homo, (20, 10), x=c60_homox, y=c60_homoy)
    fig = go.Figure()
    # binned = go.Scatter(x=c60_homoy_b - 16.8, y=c60_homo_b[7], mode='lines')
    traces = [go.Scatter3d(x=[x] * len(row), y=c60_homoy_b - 16.8, z=row, mode='lines') for row, x in
              zip(c60_homo_b.T, c60_homox_b)]  # waterfall plot
    fig.add_traces(traces)
    fig.show(renderer='browser')
