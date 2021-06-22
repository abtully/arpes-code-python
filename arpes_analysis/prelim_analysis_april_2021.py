"""
April 2021 ARPES measurement data loading and preliminary analysis
@author: Alexandra Tully
@date: April 2021
"""

import numpy as np

from arpes_functions.HDF5_loader import data_from_hdf
from arpes_functions.plotting_functions import plot3D, plot2D, transpose_figure
from arpes_functions.k_correction import kcorrect_phimotor
from arpes_functions.polygons import gen_polygon, gen_tiled_hexagons, plot_polygons, plot_polygon

"""APRIL 2021"""



if __name__ == '__main__':
    """HOMO short, Gamma, Lamp (April 2021)"""
    # fp = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\April_2021\LT\Lamp\3D\phi_motor_scan\HOMO_short_Gamma0_Y2021M04D27dT20h45m46s'
    # fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/HOMO_short_Gamma0'
    #
    # slice_val = 15.2  # Lamp
    #
    # # fns = ['HOMO15_Y2021M01D27dT21h18m41s.h5', 'HOMO15_Y2021M01D28dT00h08m04s.h5', 'HOMO15_Y2021M01D28dT02h57m28s.h5',
    # #        'HOMO15_Y2021M01D28dT05h46m53s.h5', 'HOMO15_Y2021M01D28dT08h36m19s.h5']
    # # get and average data from multiple hdf5 files
    # data, ss, cs, p = data_from_hdf(fp, fn)
    # # data_avg = avg_array_from_hdfs(fp, fns)
    # datat = np.swapaxes(data, 2, 0)
    # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # x_new, y_new, new_data = get2Dslice(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # data_region, x, y = get_data_region(new_data, xaxis=x_new, yaxis=y_new, xbounds=(-18, 15.5))
    # plot2D(x=x, y=y, data=data_region)
    # gauss_data = gaussian_filter(data_region, sigma=2)
    # plot2D(x=x, y=y, data=gauss_data)
    #
    # # # export averaged data to hdf5 file
    # # avg_data_hdf(fp, 'HOMO_short_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)
    # #
    # # # load averaged data from hdf5 file
    # # data, ss, cs, p = data_from_hdf(fp, 'HOMO_short_averaged.h5')
    # # # plot averaged data
    # # plot3D(x=ss, y=cs, z=p, data=data, slice_dim='y', slice_val=15.7)

    """HOMO overnight, Gamma, Lamp (April 2021)"""

    # fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/phi_motor_scan_overnight'
    # fn = r'HOMO_FS_overnight_Gamma0_Y2021M04D27dT21h57m10s.h5'
    #
    # # slice_val = 15.2  # Lamp: 14.4eV to 16eV
    # #
    # # fns = ['HOMO_FS_overnight_Gamma0_Y2021M04D27dT21h57m10s.h5', 'HOMO_FS_overnight_Gamma0_Y2021M04D27dT22h54m55s.h5',
    # #        'HOMO_FS_overnight_Gamma0_Y2021M04D27dT23h52m41s.h5', 'HOMO_FS_overnight_Gamma0_Y2021M04D28dT00h50m28s.h5',
    # #        'HOMO_FS_overnight_Gamma0_Y2021M04D28dT01h48m14s.h5', 'HOMO_FS_overnight_Gamma0_Y2021M04D28dT02h46m00s.h5',
    # #        'HOMO_FS_overnight_Gamma0_Y2021M04D28dT03h43m46s.h5']
    # #
    # # # get and average data from multiple hdf5 files
    # # data, ss, cs, p = data_from_hdf(fp, fn)
    # # data_avg = avg_array_from_hdfs(fp, fns)
    # # datat = np.swapaxes(data_avg, 2, 0)
    #
    # # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # # x_new, y_new, new_data = get2Dslice(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # # data_region, x, y = get_data_region(new_data, xaxis=x_new, yaxis=y_new, xbounds=(-18, 15.5))
    # # plot2D(x=x, y=y, data=data_region)
    # # gauss_data = gaussian_filter(data_region, sigma=2)
    # # plot2D(x=x, y=y, data=gauss_data)
    #
    # # # export averaged data to hdf5 file
    # # avg_data_hdf(fp, 'HOMO_overnight_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)
    #
    # # load averaged data from hdf5 file
    # data, ss, cs, p = data_from_hdf(fp, 'HOMO_overnight_averaged.h5')
    #
    # # plot averaged data
    # slice_val = 15.6
    # plot3D(x=ss, y=cs, z=p, data=data, slice_dim='y', slice_val=slice_val, int_range=0.1,
    #        title=f'HOMO overnight lamp {slice_val}eV')
    #
    # # k-correction + hexagons
    # val = np.round(slice_val - 16.8, 1)
    # kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='HOMO_overnight_averaged.h5', val=val, Eint=0.1)
    # fig = plot2D(x=kx, y=ky, data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'HOMO at gamma lamp {val}eV')
    #
    # # add hexagons
    # coords = gen_polygon(6, radius=1.1, rotation=30, translation=(0, 0.2))  # Au at -1.2 eV
    # coords = gen_polygon(6, radius=0.42, rotation=30, translation=(0, -0.5))  # C60 at HOMO
    # new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(coords, radius=0.42,
    #                                                                                      rotation=30,
    #                                                                                      translation=(0, -0.5))
    # plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r], fig=fig, color='black')
    # plot_polygon(coords, fig=fig, color='orange')

    """Fermi Surface, Gamma, Lamp (April 2021)"""

    # fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/FermiEdge_gamma0'
    # fn = r'FS_Gamma0_Y2021M04D28dT06h21m31s.h5'
    #
    # # slice_val = 16.8  # Lamp: 15.9eV to 17.5eV
    # #
    # # fns = ['FS_Gamma0_Y2021M04D28dT06h21m31s.h5', 'FS_Gamma0_Y2021M04D28dT07h19m16s.h5']
    # #
    # # # get and average data from multiple hdf5 files
    # # data, ss, cs, p = data_from_hdf(fp, fn)
    # # data_avg = avg_array_from_hdfs(fp, fns)
    # # datat = np.swapaxes(data, 2, 0)
    # #
    # # # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # # # x_new, y_new, new_data = get2Dslice(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # # # data_region, x, y = get_data_region(new_data, xaxis=x_new, yaxis=y_new, xbounds=(-18, 15.5))
    # # # plot2D(x=x, y=y, data=data_region)
    # # # gauss_data = gaussian_filter(data_region, sigma=2)
    # # # plot2D(x=x, y=y, data=gauss_data)
    # #
    # # # export averaged data to hdf5 file
    # # avg_data_hdf(fp, 'FS_gamma0', datat, p=p, slice_scale=ss, channel_scale=cs)
    # #
    # # # load averaged data from hdf5 file
    # data, ss, cs, p = data_from_hdf(fp, 'FS_gamma0.h5')
    #
    # slice_val = 16.6
    #
    # # plot averaged data
    # # plot3D(x=ss, y=cs, z=p, data=data, slice_dim='y', slice_val=slice_val, int_range=0.1,
    # #        title=f'FS at gamma lamp {slice_val}eV')
    #
    # # k-correction + hexagons
    # val = np.round(slice_val - 16.8, 1)
    # kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='FS_gamma0.h5', val=val, Eint=0.1)
    # fig = plot2D(x=kx, y=ky, data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'FS at gamma lamp {val}eV')
    #
    # # add hexagons
    # coords = gen_polygon(6, radius=0.42, rotation=30, translation=(0, -0.5))
    # new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(coords, radius=0.42,
    #                                                                                      rotation=30,
    #                                                                                      translation=(0, -0.5))
    # plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r], fig=fig, color='black')

    """HOMO, Gamma0, XUV (April 2021)"""

    # fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/XUV/3D/phi_motor_scan/XUV_HOMO_gamma0_2'
    # fn = r'XUV_HOMO_gamma0_2_Y2021M04D28dT12h37m51s.h5'
    #
    # slice_val = 17.1  # XUV @ 19H, 22.4eV
    #
    # fns = ['XUV_HOMO_gamma0_2_Y2021M04D28dT12h37m51s.h5', 'XUV_HOMO_gamma0_2_Y2021M04D28dT13h01m38s.h5',
    #        'XUV_HOMO_gamma0_2_Y2021M04D28dT13h25m24s.h5', 'XUV_HOMO_gamma0_2_Y2021M04D28dT13h49m11s.h5']
    # # get and average data from multiple hdf5 files
    # data, ss, cs, p = data_from_hdf(fp, fn)
    # data_avg = avg_array_from_hdfs(fp, fns)
    # datat = np.swapaxes(data_avg, 2, 0)
    # # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # # x_new, y_new, new_data = get2Dslice(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # # data_region, x, y = get_data_region(new_data, xaxis=x_new, yaxis=y_new, xbounds=(-18, 15.5))
    # # plot2D(x=x, y=y, data=data_region)
    # # gauss_data = gaussian_filter(data_region, sigma=2)
    # # plot2D(x=x, y=y, data=gauss_data)
    #
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

    """FS, Gamma0, XUV (April 2021)"""

    # fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/XUV/3D/phi_motor_scan/XUV_FS_gamma0'
    # fn = r'XUV_FS_gamma0_Y2021M04D28dT14h56m26s.h5'
    #
    # slice_val = 18.2  # XUV @ 19H, 22.4eV
    #
    # fns = ['XUV_FS_gamma0_Y2021M04D28dT14h56m26s.h5', 'XUV_FS_gamma0_Y2021M04D28dT15h54m10s.h5',
    #        'XUV_FS_gamma0_Y2021M04D28dT16h51m57s.h5']
    # # get and average data from multiple hdf5 files
    # data, ss, cs, p = data_from_hdf(fp, fn)
    # data_avg = avg_array_from_hdfs(fp, fns)
    # datat = np.swapaxes(data_avg, 2, 0)
    # # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # # x_new, y_new, new_data = get2Dslice(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # # data_region, x, y = get_data_region(new_data, xaxis=x_new, yaxis=y_new, xbounds=(-18, 15.5))
    # # plot2D(x=x, y=y, data=data_region)
    # # gauss_data = gaussian_filter(data_region, sigma=2)
    # # plot2D(x=x, y=y, data=gauss_data)
    #
    # # # export averaged data to hdf5 file
    # avg_data_hdf(fp, 'XUV_FS_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)
    #
    # # load averaged data from hdf5 file
    # data, ss, cs, p = data_from_hdf(fp, 'XUV_FS_averaged.h5')
    #
    # # plot averaged data
    # EF = 18.4
    # slice_val = 18.2
    # val = np.round(slice_val - EF, 2)
    # plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='y', slice_val=slice_val, int_range=0.1,
    #        title=f'FS at gamma XUV {slice_val}eV')
    #
    # # k-correction + hexagons
    # slice_val = 18.0
    # val = np.round(slice_val - 18.2, 1)
    # kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='XUV_FS_averaged.h5', val=val, Eint=0.1, EF=18.2)
    # fig = plot2D(x=kx, y=ky, data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'FS at gamma XUV {val}eV')
    #
    # # add hexagons
    # coords = gen_polygon(6, radius=1.1, rotation=30, translation=(0, 0.25))  # Au at -1.2 eV
    # coords = gen_polygon(6, radius=0.42, rotation=30, translation=(0, -0.5))  # C60 at HOMO
    # new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(coords, radius=0.42,
    #                                                                                      rotation=30,
    #                                                                                      translation=(0, -0.5))
    # plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r], fig=fig, color='black')
    # plot_polygon(coords, fig=fig, color='orange')

    """FS2, Gamma0, Lamp (April 2021)"""

    # fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/FermiEdge_2'
    # fn = r'Lamp_FermiEdge2_Y2021M04D28dT22h33m17s.h5'
    #
    # slice_val = 16.2
    #
    # fns = ['Lamp_FermiEdge2_Y2021M04D28dT22h33m17s.h5', 'Lamp_FermiEdge2_Y2021M04D28dT23h31m01s.h5',
    #        'Lamp_FermiEdge2_Y2021M04D29dT00h28m47s.h5', 'Lamp_FermiEdge2_Y2021M04D29dT01h26m32s.h5',
    #        'Lamp_FermiEdge2_Y2021M04D29dT02h24m17s.h5', 'Lamp_FermiEdge2_Y2021M04D29dT03h22m02s.h5',
    #        'Lamp_FermiEdge2_Y2021M04D29dT04h19m48s.h5', 'Lamp_FermiEdge2_Y2021M04D29dT05h17m33s.h5'
    #        ]
    # # get and average data from multiple hdf5 files
    # data, ss, cs, p = data_from_hdf(fp, fn)
    # data_avg = avg_array_from_hdfs(fp, fns)
    # datat = np.swapaxes(data_avg, 2, 0)
    # # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # # x_new, y_new, new_data = get2Dslice(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # # data_region, x, y = get_data_region(new_data, xaxis=x_new, yaxis=y_new, xbounds=(-18, 15.5))
    # # plot2D(x=x, y=y, data=data_region)
    # # gauss_data = gaussian_filter(data_region, sigma=2)
    # # plot2D(x=x, y=y, data=gauss_data)
    #
    # # # export averaged data to hdf5 file
    # avg_data_hdf(fp, 'FS2_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)
    #
    # # load averaged data from hdf5 file
    # data, ss, cs, p = data_from_hdf(fp, 'FS2_averaged.h5')
    #
    # # plot averaged data
    # EF = 16.8
    # slice_val = 16.6
    # val = np.round(slice_val - EF, 2)
    # plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='y', slice_val=slice_val, int_range=0.1,
    #        title=f'FS2 at gamma lamp {slice_val}eV')

    """HS Cuts Analysis"""

    # # load data
    # path = r'C:\Users\atully\Code\ARPES Code Python\analysis_data'
    # EF = 16.8
    #
    # # XUV
    # # data = Data2D.single_load(month='April', year='2021', light_source='XUV', filepath=path,
    # #                        filename='OMBE_XUV_2D0014_.ibw')  # XUV lots of measurements spot; G -> K -> M
    # # title = 'XUV (hv = 22.65): G -> K -> M (Common measurement spot)'
    # # data = Data2D.single_load(month='April', year='2021', light_source='XUV', filepath=path,
    # #                        filename='OMBE_XUV_2D0015_.ibw')  # XUV new spot; G -> K -> M
    # # title = 'XUV (hv = 22.65): G -> K -> M (New measurement spot)'
    # # data = Data2D.single_load(month='April', year='2021', light_source='XUV', filepath=path,
    # #                        filename='OMBE_XUV_2D0016_.ibw')  # XUV new spot, hv = 20.27 (EF = 16.2); G -> K -> M
    # # title = 'XUV (hv = 20.27): G -> K -> M (New measurement spot)'
    #
    # # Lamp
    # # data = Data2D.single_load(month='April', year='2021', light_source='Lamp', filepath=path,
    # #                        filename='OMBE_Lamp_2D0010_.ibw')  # Lamp long scan HOMO; G -> K -> M
    # # title = 'Lamp: G -> K -> M (New measurement spot)'
    # data = Data2D.single_load(month='April', year='2021', light_source='Lamp', filepath=path,
    #                        filename='OMBE_Lamp_2D0011_.ibw')  # Lamp long scan HOMO; M -> M
    # title = 'Lamp: M -> M (New measurement spot)'
    #
    # # k correction
    # x, y, d = kcorrect2D_general(data=data.data, xaxis=data.xaxis, yaxis=data.yaxis, theta0=-1)
    # xlabel = 'kx [A-1]'
    # ylabel = 'E_B'
    # plot2D(x, y - EF, d, colorscale='RdBu', xlabel=xlabel, ylabel=ylabel, title=title)
    #
    # # take gaussian mask and differentiate data
    # gauss_data = gaussian_filter(d, sigma=5)  # change sigma for different plots
    # diff_data = np.diff(gauss_data, n=2, axis=0)  # this procedure removes 2 indices from whichever axis you're
    #
    # # differentiating over, in this case yaxis.
    # new_yaxis = np.linspace(y[0], y[-1], diff_data.shape[0])  # fix yaxis
    #
    # # adjust region of data and plot
    # dd, dx, dy = get_data_region(d, xaxis=x, yaxis=y, xbounds=(-0.9, 0.8))
    # # xbounds = (-18, 15), ybounds = (-4.75, 0.2)
    # gdd, gdx, gdy = get_data_region(gauss_data, xaxis=x, yaxis=y, xbounds=(-0.9, 0.8))
    # c60, c60x, c60y = get_data_region(diff_data, xaxis=x, yaxis=new_yaxis, xbounds=(-0.9, 0.8))
    #
    # # Final Plots
    # plot2D(dx, dy - EF, dd, title=f'K-Corrected Data {title}', xlabel='kx [A-1]', ylabel='E_B')
    # plot2D(gdx, gdy - EF, gdd, title=f'Gaussian Mask (sigma = 5) {title}', xlabel='kx [A-1]', ylabel='E_B')
    # fig = plot2D(c60x, c60y - EF, c60, title=f'2nd Derivative of (Masked) Data {title}', xlabel='kx [A-1]', ylabel='E_B')
    # fig.update_layout(coloraxis=dict(cmin=-1, cmax=1))
    # fig.show()  # clipped data, everything > 1 = 1 etc.

    """Phi Motor Scan Analysis"""

    # XUV
    fp = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\April_2021\LT\XUV\3D\phi_motor_scan\XUV_FS_gamma0'

    # load averaged data from hdf5 file
    data, ss, cs, p = data_from_hdf(fp, 'XUV_FS_averaged.h5')  # XUV

    EF = 18.4

    # plot averaged data
    slice_val = 18.1
    int_range = 2
    fig = plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='y', slice_val=slice_val, int_range=0.5,
           title=f'FS at gamma XUV {slice_val}eV, int_range:{int_range}', xlabel='')
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
           title=f'FS at gamma XUV {slice_val}eV, int_range:{int_range}', xlabel='')
    transpose_figure(fig).show()
    fig = transpose_figure(fig)

    # k-correct slice in phi FIXME: Not currently working right :(
    # val = slice_val
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

    # Lamp
    fp = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\April_2021\LT\Lamp\3D\phi_motor_scan\FermiEdge_gamma0'
    data, ss, cs, p = data_from_hdf(fp, 'FS_gamma0.h5')  # Lamp
    EF = 16.8
    slice_val = 16.5
    val = np.round(slice_val - EF, 1)
    kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='FS_gamma0.h5', val=val, Eint=0.5, EF=EF, slice_dim='y',
                                      phi_m=6.5, theta_m=-1, phi0=-13, theta0=-1.2, alpha0=0)
    fig = plot2D(x=kx, y=ky, data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'{val}eV (Lamp)')
    # fig.add_vline(x=0)
    # fig.add_hline(y=0)
    # fig.show()

    coords = gen_polygon(6, radius=1.1, rotation=30, translation=(0, 0.25))  # Au at -1.2 eV
    coords = gen_polygon(6, radius=0.42, rotation=30, translation=(0, -0.72))  # C60 at HOMO
    new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(coords, radius=0.42,
                                                                                         rotation=30,
                                                                                         translation=(0, -0.72))
    plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r], fig=fig, color='black')
    plot_polygon(coords, fig=fig, color='orange')


    """Combination Plots!!!"""

    # # load averaged data from hdf5 file
    # fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/XUV/3D/phi_motor_scan/XUV_FS_gamma0'
    # data, ss, cs, p = data_from_hdf(fp, 'XUV_FS_averaged.h5')
    #
    # fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/FermiEdge_gamma0'
    # data, ss, cs, p = data_from_hdf(fp, 'FS_gamma0.h5')
    #
    # fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/XUV/3D/phi_motor_scan/XUV_HOMO_gamma0_2'
    # data, ss, cs, p = data_from_hdf(fp, 'XUV_HOMO_averaged.h5')
    #
    # fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/phi_motor_scan_overnight'
    # data, ss, cs, p = data_from_hdf(fp, 'HOMO_overnight_averaged.h5')
    #
    # # plot averaged data
    # EF = 18.4
    # slice_val = 17.7
    # val = np.round(slice_val - EF, 2)
    # plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='y', slice_val=slice_val, int_range=0.5, xlabel='Theta',
    #        ylabel='Phi', title=f'{val}eV (XUV)')
    #
    # EF = 16.8
    # slice_val = 16.5
    # val = np.round(slice_val - EF, 2)
    # plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='y', slice_val=slice_val, int_range=0.5, xlabel='Theta',
    #        ylabel='Phi', title=f'{val}eV (Lamp)')
    #
    # # from Typing import List
    # # heatmap_info is list of lists: fp, fn, slice_val, EF
    # def getheatmaps(heatmap_info: list, int_range=0.5, color='Plasma'):
    #     heatmaps = []
    #     for h in heatmap_info:
    #         data, ss, cs, p = data_from_hdf(h[0], h[1])
    #         fig = plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='y', slice_val=h[3], int_range=int_range,
    #                      show=False, color=color)
    #         hm = fig.data[0]
    #         heatmaps.append(hm)
    #     return heatmaps
    #
    #
    # def multiplot_3D(heatmap_info: list, rows, cols, int_range=0.5, xlabel=None, ylabel=None, title=None, color='Plasma'):
    #     fig = make_subplots(rows=rows, cols=cols, shared_xaxes=True, shared_yaxes=True)
    #     coords = list(product(range(1, rows + 1), range(1, cols + 1)))
    #     heatmaps = getheatmaps(heatmap_info=heatmap_info, int_range=int_range, color=color)
    #     for coords, heatmap in zip(coords, heatmaps):
    #         fig.add_trace(heatmap, row=coords[0], col=coords[1])
    #
    # hm_info = ['/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/XUV/3D/phi_motor_scan/XUV_HOMO_gamma0_2',
    #            'XUV_HOMO_averaged.h5', 16.9, 18.2], \
    #           ['/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/phi_motor_scan_overnight']\
    #     , []
    #
    #
    #
    #
    # # def multiplot_2D()


