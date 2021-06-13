import numpy as np
import os
import plotly.graph_objects as go

from ali_HDF5_loader import data_from_hdf, avg_array_from_hdfs, avg_data_hdf
from ali_plotting_functions import plot3D, plot2D, transpose_figure
from ali_analysis_functions import get_data_region, get_vertical_slice, get_averaged_slice, get_horizontal_slice
from ali_k_correction import kcorrect_phimotor, get2Dslice, kcorrect2D_general, kcorrect2D, fix_EkatEF
from ali_polygons import gen_polygon, gen_tiled_hexagons, plot_polygons, plot_polygon
from ali_classes import Data2D

from scipy.ndimage import gaussian_filter

if __name__ == '__main__':
    # XUV_FS_gamma0 (averaged)
    fp = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\April_2021\LT\XUV\3D\phi_motor_scan\XUV_FS_gamma0'

    # load averaged data from hdf5 file
    data, ss, cs, p = data_from_hdf(fp, 'XUV_FS_averaged.h5')  # XUV

    EF = 18.4

    # plot averaged data
    slice_val = 18.1
    int_range = 2
    fig = plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='y', slice_val=slice_val, int_range=0.5,
                 title=f'FS at gamma XUV {slice_val}eV, int_range:{int_range}', xlabel='', show=False)
    fig.add_hline(y=-4.5)
    fig.add_vline(x=-4)
    fig.show()

    # get slice in phi
    slice_val = -4.5
    int_range = 2
    # val = np.round(slice_val - EF, 2)
    plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='z', slice_val=slice_val, int_range=int_range,
           title=f'FS at gamma XUV {slice_val}eV, int_range:{int_range}', xlabel='')

    # get slice in theta
    slice_val = -4
    int_range = 1
    fig = plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='x', slice_val=slice_val, int_range=int_range,
                 title=f'FS at gamma XUV {slice_val}eV, int_range:{int_range}', xlabel='', show=False)
    transpose_figure(fig).show()
    fig = transpose_figure(fig)

    # get horizontal slice from "slice in theta" dataset
    fig_data = fig.data[0]
    x, y, d = fig_data['x'], fig_data['y'], fig_data['z']
    # # take gaussian mask and differentiate data
    #
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

    """Fit Gaussians to Vertical Slices of Data"""
    # testdata = get_averaged_slice(get_vertical_slice(d, x, x[1], 1), axis='x')
    testdata2 = get_averaged_slice(get_horizontal_slice(d, y, 18, 0.1), axis='y')
    fig.add_hline(y=18)
    fig.add_hline(y=18.05, line_dash='dot')
    fig.add_hline(y=17.95, line_dash='dot')
    fig.show()
    # fig.add_trace(go.Scatter(x=x, y=testdata2))
    fig2 = go.Figure(data=go.Scatter(x=x, y=testdata2))
    fig2.update_layout(title='Horizontal Chunk')
    fig2.show(renderer='browser')

    data = d
    y = y - 16.8
    import lmfit as lm
    from ali_analysis_functions import gaussian

    # model = lm.models.GaussianModel()
    model = lm.Model(gaussian)
    params = model.make_params()
    params.add('amplitude', 200, min=np.min(data), max=np.max(data))
    params.add('center', -1.6, min=np.min(y), max=np.max(y))
    params.add('width', 1, min=0.001, max=np.max(y) - np.min(y))
    params.add('const', 150, min=0.5 * np.min(data), max=np.max(data))
    fit = model.fit(testdata2.astype(np.float32), x=y.astype(np.float32), params=params)
    print(fit.fit_report())
    fit.plot()
    fit.best_values
    fit.best_values['center']

    fits = []
    for row in data.T:
        fit = model.fit(row.astype(np.float32), x=y.astype(np.float32), params=params)
        fits.append(fit)

    centers = [fit.best_values['center'] for fit in fits]
    fig.add_trace(go.Scatter(x=x, y=centers))
    fig.update_layout(title='Gaussian Mask (sigma = 5) with Fit')
    fig.show(renderer='browser')

    # k-correct slice in phi FIXME: Not currently working right :(
    val = slice_val
    val = -4
    kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='XUV_FS_averaged.h5', val=val, Eint=2, EF=18.4, slice_dim='x',
                                      phi_m=6.5, theta_m=-1, phi0=-13, theta0=-0.6, alpha0=0)
    fig = plot2D(x=kx, y=np.flip(ky), data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'{val}eV (XUV)')

    # k-correction + hexagons
    slice_val = 17.7
    # slice_val = 6.5
    val = np.round(slice_val - EF, 1)
    kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='XUV_FS_averaged.h5', val=val, Eint=0.5, EF=18.4, slice_dim='y',
                                      phi_m=6.5, theta_m=-1, phi0=-13, theta0=-0.6, alpha0=0)
    fig = plot2D(x=kx, y=ky, data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'{val}eV (XUV)')
    # fig.add_vline(x=0)
    # fig.add_hline(y=0)
    # fig.show()