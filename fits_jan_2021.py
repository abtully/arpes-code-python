"""
Fitting ARPES data from January 2021
@author: Alexandra Tully
@date: June 2021
"""

import numpy as np
import plotly.graph_objects as go

from HDF5_loader import data_from_hdf
from plotting_functions import plot3D, plot2D, transpose_figure
from analysis_functions import get_data_region, get_vertical_slice, get_averaged_slice, get_horizontal_slice
from k_correction import kcorrect_phimotor, get2Dslice, kcorrect2D_general, kcorrect2D, fix_EkatEF
from arpes_dataclasses import Data2D
import fitting_functions as ff
from misc_functions import line_intersection
from leed_analysis import get_line_distance

""" January Analysis """

# Load Data #

# raw data
# d = Data2D.single_load('January', year='2021', light_source='XUV', filename='OMBE_XUV_2D0004_.ibw')
data = Data2D.single_load(month='January', year='2021', light_source='XUV', scan_number=4)
plot2D(data.xaxis, data.yaxis, data.data, title='OMBE_XUV_2D0004_.ibw, January 2021', xlabel='Theta', ylabel='KE')

# zoom in on cone
d, x, y = get_data_region(data.data, data.xaxis, data.yaxis, xbounds=(-12, 8), ybounds=(17.1, 18.4), EB=False)
fig = plot2D(x, y, d, title='OMBE_XUV_2D0004_.ibw, January 2021', xlabel='Theta', ylabel='KE')

# Apply Filters #

# # apply gaussian mask
# sigma = 2
# gauss_data = gaussian_filter(d, sigma=sigma)
# fig = plot2D(x, y, gauss_data,
#              title=f'OMBE_XUV_2D0004_.ibw, January 2021, Gaussian mask (sigma: {sigma})',
#              xlabel='Theta', ylabel='KE')

# # increase signal:noise ratio
# low_bound = np.asarray(np.where(gauss_data > 20))
# up_bound = np.asarray(np.where(low_bound < 50))
# scaled_data = 5 * gauss_data[up_bound]

# Run Fits #

# fit lorentzians
fits = []
coords = []
num_peaks = 2
for yval in np.linspace(17.3, 18.35, 10):
    row = get_averaged_slice(get_horizontal_slice(d, y, yval, 0.1), axis='y')
    fit = ff.fit_lorentzian_data(x=x, data=row,
                                 num_peaks=num_peaks, amplitudes=30, centers=[-6, -1], sigmas=1,
                                 offset_type='constant')
    fits.append(fit)
    coords.extend([(fit.best_values[f'i{i}_center'], yval) for i in range(num_peaks)])

fig.add_trace(go.Scatter(x=[c[0] for c in coords], y=[c[1] for c in coords],
                         mode='markers', marker=dict(color='black'), name='lorentzian peak centers'))
fig.show()

# plot figure of waterfall lorentzian fits
fig3 = go.Figure()
spacing = 20
interval = 1  # takes every xth line to plot
for i, fit in enumerate(fits[::interval]):
    x = fit.userkws['x']
    fitdata = fit.data
    fig3.add_trace(go.Scatter(x=x, y=fitdata + (spacing * i), mode='lines'))
    fig3.add_trace(go.Scatter(x=x, y=fit.eval(x=x) + (spacing * i), mode='lines'))

fig3.show()

# fit hyperbola
hyper_x = np.asarray([c[0] for c in coords])
hyper_y = np.asarray([c[1] for c in coords])

fit = ff.fit_hyperbola_data(x=hyper_x[np.where(hyper_y < 18.1)], data=hyper_y[np.where(hyper_y < 18.1)],
                            num=1,
                            aes=1, bes=1, hes=-4, kes=18.5,
                            offset_type=None)
# fit = ff.fit_hyperbola_data(x=hyper_x, data=hyper_y,
#                             num=1,
#                             aes=1, bes=1, hes=-4, kes=18.5,
#                             offset_type=None)

# plot hyperbola
fig.add_trace(go.Scatter(x=x[np.where(fit.eval(x=x) > y[0])],
                         y=fit.eval(x=x)[np.where(fit.eval(x=x) > y[0])],
                         mode='lines', name='hyperbolic fit', line=dict(color='red')))
hyper_pos = ff.hyperbola(x=x,
                         a=fit.best_values['i0_a'], b=fit.best_values['i0_b'],
                         h=fit.best_values['i0_h'], k=fit.best_values['i0_k'],
                         pos=True)
fig.add_trace(go.Scatter(x=x[np.where(hyper_pos < 19.5)], y=hyper_pos[np.where(hyper_pos < 19.5)],
                         mode='lines', name='hyperbolic fit', line=dict(color='red')))

# plot hyperbola line
y1, y2 = ff.make_hyperbola_asymptotes(fit, x=x)
fig.add_trace(go.Scatter(x=x[np.where(np.logical_and(y1 > y[0], y1 < 19.5))],
                         y=y1[np.where(np.logical_and(y1 > y[0], y1 < 19.5))],
                         mode='lines', line=dict(color='red', dash='dash'), name='hyperbola asymptote'))
fig.add_trace(go.Scatter(x=x[np.where(np.logical_and(y2 > y[0], y2 < 19.5))],
                         y=y2[np.where(np.logical_and(y2 > y[0], y2 < 19.5))],
                         mode='lines', line=dict(color='red', dash='dash'), name='hyperbola asymptote'))

# find intersection of asymptotes
fit1 = ff.fit_linear_data(x, y1, num=1)
fit2 = ff.fit_linear_data(x, y2, num=1)
a1, b1 = line_intersection(fit1, fit2)
fig.add_trace(go.Scatter(x=a1, y=b1,
                         mode='markers', marker=dict(color='green', symbol='circle'), name='intersection'))

# fit lines
threshold = -3.7
fit1 = ff.fit_linear_data(x=hyper_x[np.where(hyper_x < threshold)], data=hyper_y[np.where(hyper_x < threshold)],
                          num=1,
                          aes=1, bes=1,
                          offset_type=None)
fig.add_trace(go.Scatter(x=x[np.where(np.logical_and(fit1.eval(x=x) > y[0], fit1.eval(x=x) < 19.5))],
                         y=fit1.eval(x=x)[np.where(np.logical_and(fit1.eval(x=x) > y[0], fit1.eval(x=x) < 19.5))],
                         mode='lines', name='linear fit', line=dict(color='turquoise', dash='dash')))

fit2 = ff.fit_linear_data(x=hyper_x[np.intersect1d(np.where(hyper_x > threshold), np.where(hyper_x < -0.5))],
                          data=hyper_y[np.intersect1d(np.where(hyper_x > threshold), np.where(hyper_x < -0.5))],
                          num=1,
                          offset_type=None)
fig.add_trace(go.Scatter(x=x[np.where(np.logical_and(fit2.eval(x=x) > y[0], fit2.eval(x=x) < 19.5))],
                         y=fit2.eval(x=x)[np.where(np.logical_and(fit2.eval(x=x) > y[0], fit2.eval(x=x) < 19.5))],
                         mode='lines', name='linear fit', line=dict(color='turquoise', dash='dash')))

# find intersection of linear fits
a2, b2 = line_intersection(fit1, fit2)
fig.add_trace(go.Scatter(x=a2, y=b2,
                         mode='markers', marker=dict(color='lightgreen', symbol='circle'), name='intersection'))

# plot figure
fig.update_layout(title='Fitted Data', xaxis_title='Theta', yaxis_title='KE')
fig.update_layout(coloraxis_colorbar=dict(len=0.5))
fig.show()

# calculate energy gap between intersections
gap = get_line_distance((a1, b1), (a2, b2))
# TODO: GAP = 0.1855eV according to this fit


# fit.params['i0_a'].value
# print(fit.fit_report())
# print(fit.ci_report())
# fig.data = fig.data[:-1]  # removes last trace addition to figure