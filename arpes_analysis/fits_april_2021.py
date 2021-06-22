"""
Fitting ARPES data from April 2021
@author: Alexandra Tully
@date: June 2021
"""

import numpy as np
import plotly.graph_objects as go

from arpes_functions.HDF5_loader import data_from_hdf
from arpes_functions.plotting_functions import plot3D, plot2D, transpose_figure
from arpes_functions.analysis_functions import get_data_region, get_averaged_slice, get_horizontal_slice
from arpes_functions import fitting_functions as ff

""" April Analysis """
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

""" Fitting"""
# get new region for fitting
fig_data = fig.data[0]
x, y, d = fig_data['x'], fig_data['y'], fig_data['z']

d, x, y = get_data_region(d, x, y, xbounds=(-35, 8))
fig = plot2D(x, y, d,
             title=f'FS at gamma XUV {slice_val}eV, int_range:{int_range}', xlabel='Theta', ylabel='Energy')

""" Fit Lorentzians to Row of Data """
# add lines to figure to show averaged region
data_for_fit = get_averaged_slice(get_horizontal_slice(d, y, 18, 0.1), axis='y')
fig.add_hline(y=18)
fig.add_hline(y=18.05, line_dash='dot')
fig.add_hline(y=17.95, line_dash='dot')
fig.show()

fig2 = go.Figure(data=go.Scatter(x=x, y=data_for_fit, name='raw data'))
fig2.update_layout(title='Horizontal Chunk')
fig2.show(renderer='browser')

fit = ff.fit_lorentzian_data(x=x, data=data_for_fit,
                             num_peaks=3, amplitudes=9, centers=[-23, -13, -3], sigmas=1,
                             offset_type='linear')
print(fit.fit_report())

fig2.add_trace(go.Scatter(x=x, y=fit.eval(x=x), mode='lines', name='fit'))
fig2.update_layout(title='Lorentzian Fit, 3 Peaks', xaxis_title='Theta', yaxis_title='Amplitude')
fig2.show()

# fig2.add_trace(go.Scatter(x=x, y=fit.eval(x=x, params=fit.init_params), mode='lines'))
# fig2.show()
#
# print(fit.init_params)

""" Fit Lorentzians to Multiple Slices of Data """

fig = plot2D(x, y, d,
             title=f'FS at gamma XUV {slice_val}eV, int_range:{int_range}', xlabel='Theta', ylabel='Energy')
fits = []
coords = []
num_peaks = 3
# for row, y_ in zip(d, y):
#     if 17.8 < y_ < 18.25:
for yval in np.linspace(17.8, 18.25, 20):
    row = get_averaged_slice(get_horizontal_slice(d, y, yval, 0.1), axis='y')
    fit = ff.fit_lorentzian_data(x=x, data=row,
                                 num_peaks=num_peaks, amplitudes=9, centers=[-23, -13, -3], sigmas=1,
                                 offset_type='linear')
    fits.append(fit)
    coords.extend([(fit.best_values[f'i{i}_center'], yval) for i in range(num_peaks)])

fig.add_trace(go.Scatter(x=[c[0] for c in coords], y=[c[1] for c in coords],
                         mode='markers', marker=dict(color='black'), name='peak centers'))
fig.show()

# fig.data = fig.data[:-1]  # removes last trace addition to figure

fig3 = go.Figure()
spacing = 2
interval = 2  # takes every xth line to plot
for i, fit in enumerate(fits[::interval]):
    x = fit.userkws['x']
    data = fit.data
    fig3.add_trace(go.Scatter(x=x, y=data + (spacing * i), mode='lines'))
    fig3.add_trace(go.Scatter(x=x, y=fit.eval(x=x) + (spacing * i), mode='lines'))

fig3.show()

""" Fit Hyperbolas on Lorentzian Center Coords """

# get new region for fitting
x, y, d = fig_data['x'], fig_data['y'], fig_data['z']

d, x, y = get_data_region(d, x, y, xbounds=(-8, 2))
fig = plot2D(x, y, d,
             title=f'FS at gamma XUV {slice_val}eV, int_range:{int_range}', xlabel='Theta', ylabel='Energy')

# fit lorentzians
fits = []
coords = []
num_peaks = 2
for yval in np.linspace(y[0], 18.25, 10):
    row = get_averaged_slice(get_horizontal_slice(d, y, yval, 0.1), axis='y')
    fit = ff.fit_lorentzian_data(x=x, data=row,
                                 num_peaks=num_peaks, amplitudes=9, centers=[-6, -1], sigmas=1,
                                 offset_type='linear')
    fits.append(fit)
    coords.extend([(fit.best_values[f'i{i}_center'], yval) for i in range(num_peaks)])

fig.add_trace(go.Scatter(x=[c[0] for c in coords], y=[c[1] for c in coords],
                         mode='markers', marker=dict(color='black'), name='peak centers'))
fig.show()

# fig.data = fig.data[:-1]  # removes last trace addition to figure

fig3 = go.Figure()
spacing = 2
interval = 1  # takes every xth line to plot
for i, fit in enumerate(fits[::interval]):
    x = fit.userkws['x']
    data = fit.data
    fig3.add_trace(go.Scatter(x=x, y=data + (spacing * i), mode='lines'))
    fig3.add_trace(go.Scatter(x=x, y=fit.eval(x=x) + (spacing * i), mode='lines'))

fig3.show()

# fit hyperbola
hyper_x = np.asarray([c[0] for c in coords])
hyper_y = np.asarray([c[1] for c in coords])

fit = ff.fit_hyperbola_data(x=hyper_x[np.where(hyper_y < 18.05)], data=hyper_y[np.where(hyper_y < 18.05)],
                            num=1,
                            aes=29, bes=20, hes=-4, kes=18.5,
                            offset_type=None)

# print(fit.fit_report())
# print(fit.ci_report())
fig.add_trace(go.Scatter(x=x, y=fit.eval(x=x), mode='lines', name='fit'))
hyper_pos = ff.hyperbola(x=x,
                         a=fit.best_values['i0_a'], b=fit.best_values['i0_b'],
                         h=fit.best_values['i0_h'], k=fit.best_values['i0_k'],
                         pos=True)
fig.add_trace(go.Scatter(x=x, y=hyper_pos, mode='lines', name='fit'))
# fig.add_trace(go.Scatter(x=hyper_x, y=fit.eval(x=hyper_x), mode='markers', name='fit'))

# using scipy optimize rather than lmfit
# popt, pcov = optimize.curve_fit(hyperbola, hyper_x, hyper_y)
# hyper_fit_data = hyperbola(x, popt[0], popt[1], popt[2], popt[3])

# fig.add_trace(go.Scatter(x=x, y=hyper_fit_data, mode='lines', name='fit'))

fig.update_layout(title='Hyperbola Fit', xaxis_title='Theta', yaxis_title='Amplitude')
fig.show()

# TODO: try gaussian mask and way more / way fewer points to fit?


# # k-correct slice in phi FIXME: Not currently working right :(
# val = slice_val
# val = -4
# kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='XUV_FS_averaged.h5', val=val, Eint=2, EF=18.4, slice_dim='x',
#                                   phi_m=6.5, theta_m=-1, phi0=-13, theta0=-0.6, alpha0=0)
# fig = plot2D(x=kx, y=np.flip(ky), data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'{val}eV (XUV)')