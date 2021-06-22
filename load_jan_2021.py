from plotting_functions import plot3D, plot2D
from analysis_functions import get_data_region
from k_correction import get2Dslice, fix_EkatEF, kcorrect3D
from scipy.ndimage import gaussian_filter
from k_correction import kcorrect_phimotor, get2Dslice
from polygons import gen_polygon, gen_tiled_hexagons, plot_polygons, plot_polygon

"""JANUARY 2021"""

"""HOMO, Theta=15"""

# fp = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\January_2021\LT\Lamp\3D\phi_motor_scan\HOMO15_Y2021M01D27'
# fn = r'HOMO15_Y2021M01D27dT21h18m41s.h5'
# fns = ['HOMO15_Y2021M01D27dT21h18m41s.h5', 'HOMO15_Y2021M01D28dT00h08m04s.h5', 'HOMO15_Y2021M01D28dT02h57m28s.h5',
#        'HOMO15_Y2021M01D28dT05h46m53s.h5', 'HOMO15_Y2021M01D28dT08h36m19s.h5']
#
# # get and average data from multiple hdf5 files
# data, ss, cs, p = data_from_hdf(fp, fn)
# data_avg = avg_array_from_hdfs(fp, fns)
# datat = np.swapaxes(data_avg, 2, 0)
# plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=15.7)
#
# # export averaged data to hdf5 file
# avg_data_hdf(fp, 'HOMO15_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)
#
# # load averaged data from hdf5 file
# data, ss, cs, p = data_from_hdf(fp, 'HOMO15_averaged.h5')
# # plot averaged data
# plot3D(x=ss, y=cs, z=p, data=data, slice_dim='y', slice_val=15.7)

"""HOMO, Theta=0"""

# fp = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\January_2021\LT\Lamp\3D\phi_motor_scan\HOMO0_Y2021M01D27'
# fn = r'HOMO0_Y2021M01D27dT18h17m03s.h5'
#
# # load data
# data, ss, cs, p = data_from_hdf(fp, fn)
# # fix data
# datat = np.swapaxes(data, 2, 0)
# test = np.flip(datat, 0)
# # plot data
# plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=15.7)
# plot3D(x=ss, y=cs, z=p, data=test, slice_dim='y', slice_val=15.7)

"""FS, Theta=?"""

# fp = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\January_2021\LT\Lamp\3D\phi_motor_scan\FS_Y2021M01D25'
# fn = r'FS_Y2021M01D25dT22h39m08s.h5'
# fns = ['FS_Y2021M01D25dT22h39m08s.h5', 'FS_Y2021M01D26dT01h28m27s.h5', 'FS_Y2021M01D26dT04h17m48s.h5']
#
# # get and average data from multiple hdf5 files
# data, ss, cs, p = data_from_hdf(fp, fn)
# data_avg = avg_array_from_hdfs(fp, fns)
# datat = np.swapaxes(data_avg, 2, 0)
# plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=15.7)
#
# # export averaged data to hdf5 file
# avg_data_hdf(fp, 'FS_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)

# # load averaged data from hdf5 file
# data, ss, cs, p = data_from_hdf(fp, 'FS_averaged.h5')
# # plot averaged data
# fig = plot3D(x=ss, y=cs, z=p, data=data, slice_dim='y', slice_val=15.7)

"""K-corrected FS, Theta=?"""

# # correct EF
# fix_EkatEF(data, 16.8)
# # generate and plot k corrected data
# kx, ky, kdata = kcorrect3D(self=d, val=15)  # mine
# plot2D(kx, ky, kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title='Constant Energy Slice')  # mine

"""Add Hexagons"""

# coords = gen_polygon(6, radius=0.42, rotation=30, translation=(-0.07, 0.1))
# new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(coords, radius=0.42,
#                                                                                      rotation=30,
#                                                                                      translation=(-0.07, 0.1))
# plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r], fig=fig)

# # datat = np.moveaxis(data, 2, 0)
# datat = np.swapaxes(data_avg, 2, 0)
