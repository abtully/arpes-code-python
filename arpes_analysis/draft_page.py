

"""Fitting Functions"""
## take gaussian mask and differentiate data

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
#
# """Fit Gaussians to Vertical Slices of Data"""
# # testdata = get_averaged_slice(get_vertical_slice(d, x, x[1], 1), axis='x')
# testdata2 = get_averaged_slice(get_horizontal_slice(d, y, 18, 0.1), axis='y')
# fig.add_hline(y=18)
# fig.add_hline(y=18.05, line_dash='dot')
# fig.add_hline(y=17.95, line_dash='dot')
# fig.show()
# # fig.add_trace(go.Scatter(x=x, y=testdata2))
# fig2 = go.Figure(data=go.Scatter(x=x, y=testdata2))
# fig2.update_layout(title='Horizontal Chunk')
# fig2.show(renderer='browser')
#
# data = d
# y = y - 16.8
# import lmfit as lm
# from ali_analysis_functions import gaussian
#
# # model = lm.models.GaussianModel()
# model = lm.Model(gaussian)
# params = model.make_params()
# params.add('amplitude', 200, min=np.min(data), max=np.max(data))
# params.add('center', -1.6, min=np.min(y), max=np.max(y))
# params.add('width', 1, min=0.001, max=np.max(y) - np.min(y))
# params.add('const', 150, min=0.5 * np.min(data), max=np.max(data))
# fit = model.fit(testdata2.astype(np.float32), x=y.astype(np.float32), params=params)
# print(fit.fit_report())
# fit.plot()
# fit.best_values
# fit.best_values['center']
#
# fits = []
# for row in data.T:
#     fit = model.fit(row.astype(np.float32), x=y.astype(np.float32), params=params)
#     fits.append(fit)
#
# centers = [fit.best_values['center'] for fit in fits]
# fig.add_trace(go.Scatter(x=x, y=centers))
# fig.update_layout(title='Gaussian Mask (sigma = 5) with Fit')
# fig.show(renderer='browser')
#
# # k-correct slice in phi FIXME: Not currently working right :(
# val = slice_val
# val = -4
# kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='XUV_FS_averaged.h5', val=val, Eint=2, EF=18.4, slice_dim='x',
#                                   phi_m=6.5, theta_m=-1, phi0=-13, theta0=-0.6, alpha0=0)
# fig = plot2D(x=kx, y=np.flip(ky), data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'{val}eV (XUV)')
#
# # k-correction + hexagons
# slice_val = 17.7
# # slice_val = 6.5
# val = np.round(slice_val - EF, 1)
# kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='XUV_FS_averaged.h5', val=val, Eint=0.5, EF=18.4, slice_dim='y',
#                                   phi_m=6.5, theta_m=-1, phi0=-13, theta0=-0.6, alpha0=0)
# fig = plot2D(x=kx, y=ky, data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'{val}eV (XUV)')
# # fig.add_vline(x=0)
# # fig.add_hline(y=0)
# # fig.show()


"""Miscellaneous Functions"""
# data = np.random.random((10, 10))
# # data2 = np.random.random((3, 3)) * 2
# x = np.linspace(0, 1, 10)
# E_f = 16.8
# HS = Data2D.single_load('October', k_cut='GM')
# # x, y, data = bin_data2D(HS, 5)
# x, y, data2 = bin_2D(HS, 5, 5)
# # x = np.linspace(HS.xaxis[0], HS.xaxis[-1], data2.shape[0])
#
# fig = plot2D(x, y, data2, colorscale='plasma')
# plot2D(HS.yaxis, HS.xaxis, HS.data)

# data = np.random.random((100, 100))
# fig = go.Figure()

# plot2D(HS.yaxis, HS.xaxis - E_f, HS.data, opacity=0.5, xlabel='Theta', ylabel='E-E_F')

# import os
# from arpes_dataclasses import Data3D
# from plotting_functions import plot3D
#
# path = os.path.abspath('/Users/alexandratully/Desktop/ARPES Data/')
#
# data = Data3D.single_load(month='April', year='2021', light_source='XUV', cryo_temp='LT',
#                           scan_type='deflector_scan', scan_number=1,
#                           filepath=path)
# x, y, z, d = data.xaxis, data.yaxis, data.zaxis, data.data
# slice_val = -1.6
# E_f = 18.2
# fig = plot3D(y, x, z, data=np.moveaxis(d, 2, 1),
#              slice_dim='z', slice_val=slice_val + E_f, int_range=0.02, show=True,
#              colorscale='Plasma', xlabel='Theta (Deg)', ylabel='Phi (Deg)')

"""Plotting Functions"""
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