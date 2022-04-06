"""
April 2022 ARPES measurement data loading and preliminary analysis
@author: Alexandra Tully
@date: April 2022
"""

import numpy as np
from scipy.ndimage import gaussian_filter

from arpes_functions.HDF5_loader import data_from_hdf_2022, avg_array_from_hdfs, data_from_hdf_2D_2022
from arpes_functions.plotting_functions import plot3D, plot2D, transpose_figure, heatmap
from arpes_functions.analysis_functions import get_2Dslice, get_data_region
from arpes_functions.k_correction import kcorrect_phimotor
from arpes_functions.polygons import gen_polygon, gen_tiled_hexagons, plot_polygons, plot_polygon

"""APRIL 2022"""

"""HOMO, Gamma0, XUV @ 22.8eV(ish), EF = 18.5 (April 2022)"""

fp = r'E:\atully\raw arpes data\April_2022\FS1'
fn = r'FS1_000_1.h5'

slice_val = 16.8  # XUV @ 19H, 22.4eV

# get and average data from multiple hdf5 files
fns = ['FS1_000_1.h5', 'FS1_001_2.h5']
data, ss, cs, p = data_from_hdf_2022(fp, fn)
p = p.flatten()
data_avg = avg_array_from_hdfs(fp, fns)
datat = np.swapaxes(data_avg, 2, 0)
plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)

# zoom in on data
x_new, y_new, new_data = get_2Dslice(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
data_region, x, y = get_data_region(new_data, xaxis=x_new, yaxis=y_new, xbounds=(-20, 17.1))
plot2D(x=x, y=y, data=data_region)

# gaussian mask
gauss_data = gaussian_filter(data_region, sigma=2)
plot2D(x=x, y=y, data=gauss_data)


fp = r'E:\atully\raw arpes data\April_2022\KE Scans'
fn = r'KE0.h5'

#TODO: FIXED IN JUPYTER NOTEBOOK
data, ss, cs = data_from_hdf_2D_2022(fp, fn)
hm = heatmap(ss, cs, data)
import plotly.graph_objects as go
fig = go.Figure()
hm = heatmap(ss, cs, data)
fig.add_trace(hm)
plot2D(x=ss, y=cs, data=data)

# diff_data = np.diff(gauss_data, n=2, axis=0)  # this procedure removes 2 indices from whichever axis you're
# # differentiating over, in this case yaxis.
# new_yaxis = np.linspace(d.yaxis[0], d.yaxis[-1], diff_data.shape[0])  # fix yaxis


"""More useful commands..."""
# # # export averaged data to hdf5 file
# avg_data_hdf(fp, 'XUV_HOMO_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)
#
# # load averaged data from hdf5 file
# data, ss, cs, p = data_from_hdf(fp, 'XUV_HOMO_averaged.h5')
#
# # plot averaged data
# plot3D(x=ss, y=cs, z=p, data=data, slice_dim='y', slice_val=slice_val, int_range=0.1,
#        title=f'HOMO at gamma XUV {slice_val}eV')
#
# # k-correction + hexagons
# slice_val = 16.6
# EF = 18.3
# val = np.round(slice_val - EF, 1)
# kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='XUV_HOMO_averaged.h5', val=val, Eint=0.1, EF=EF)
# fig = plot2D(x=kx, y=ky, data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'HOMO at gamma XUV {val}eV')
#
# # add hexagons
# coords = gen_polygon(6, radius=1.1, rotation=30, translation=(0, 0.25))  # Au at -1.2 eV
# coords = gen_polygon(6, radius=0.42, rotation=30, translation=(0, -0.5))  # C60 at HOMO
# new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(coords, radius=0.42,
#                                                                                      rotation=30,
#                                                                                      translation=(0, -0.5))
# plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r], fig=fig, color='black')
# plot_polygon(coords, fig=fig, color='orange')
