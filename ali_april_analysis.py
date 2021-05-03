import numpy as np
from ali_HDF5_loader import data_from_hdf, avg_array_from_hdfs, avg_data_hdf
from ali_plotting_functions import plot3D, plot2D
from ali_analysis_functions import get_data_region
from ali_k_correction import get2Dslice
from scipy.ndimage import gaussian_filter
from ali_k_correction import kcorrect_phimotor, get2Dslice
from ali_polygons import gen_polygon, gen_tiled_hexagons, plot_polygons, plot_polygon

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

    fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/phi_motor_scan_overnight'
    fn = r'HOMO_FS_overnight_Gamma0_Y2021M04D27dT21h57m10s.h5'

    # slice_val = 15.2  # Lamp: 14.4eV to 16eV
    #
    # fns = ['HOMO_FS_overnight_Gamma0_Y2021M04D27dT21h57m10s.h5', 'HOMO_FS_overnight_Gamma0_Y2021M04D27dT22h54m55s.h5',
    #        'HOMO_FS_overnight_Gamma0_Y2021M04D27dT23h52m41s.h5', 'HOMO_FS_overnight_Gamma0_Y2021M04D28dT00h50m28s.h5',
    #        'HOMO_FS_overnight_Gamma0_Y2021M04D28dT01h48m14s.h5', 'HOMO_FS_overnight_Gamma0_Y2021M04D28dT02h46m00s.h5',
    #        'HOMO_FS_overnight_Gamma0_Y2021M04D28dT03h43m46s.h5']
    #
    # # get and average data from multiple hdf5 files
    # data, ss, cs, p = data_from_hdf(fp, fn)
    # data_avg = avg_array_from_hdfs(fp, fns)
    # datat = np.swapaxes(data_avg, 2, 0)

    # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # x_new, y_new, new_data = get2Dslice(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # data_region, x, y = get_data_region(new_data, xaxis=x_new, yaxis=y_new, xbounds=(-18, 15.5))
    # plot2D(x=x, y=y, data=data_region)
    # gauss_data = gaussian_filter(data_region, sigma=2)
    # plot2D(x=x, y=y, data=gauss_data)

    # # export averaged data to hdf5 file
    # avg_data_hdf(fp, 'HOMO_overnight_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)

    # load averaged data from hdf5 file
    data, ss, cs, p = data_from_hdf(fp, 'HOMO_overnight_averaged.h5')

    # plot averaged data
    slice_val = 15.6
    plot3D(x=ss, y=cs, z=p, data=data, slice_dim='y', slice_val=slice_val, int_range=0.1,
           title=f'HOMO overnight lamp {slice_val}eV')

    # k-correction + hexagons
    val = np.round(slice_val - 16.8, 1)
    kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='HOMO_overnight_averaged.h5', val=val, Eint=0.1)
    fig = plot2D(x=kx, y=ky, data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'HOMO at gamma lamp {val}eV')

    # add hexagons
    coords = gen_polygon(6, radius=1.1, rotation=30, translation=(0, 0.2))  # Au at -1.2 eV
    coords = gen_polygon(6, radius=0.42, rotation=30, translation=(0, -0.5))  # C60 at HOMO
    new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(coords, radius=0.42,
                                                                                         rotation=30,
                                                                                         translation=(0, -0.5))
    plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r], fig=fig, color='black')
    plot_polygon(coords, fig=fig, color='orange')

    """Fermi Surface, Gamma, Lamp (April 2021)"""

    fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/FermiEdge_gamma0'
    fn = r'FS_Gamma0_Y2021M04D28dT06h21m31s.h5'

    # slice_val = 16.8  # Lamp: 15.9eV to 17.5eV
    #
    # fns = ['FS_Gamma0_Y2021M04D28dT06h21m31s.h5', 'FS_Gamma0_Y2021M04D28dT07h19m16s.h5']
    #
    # # get and average data from multiple hdf5 files
    # data, ss, cs, p = data_from_hdf(fp, fn)
    # data_avg = avg_array_from_hdfs(fp, fns)
    # datat = np.swapaxes(data, 2, 0)
    #
    # # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # # x_new, y_new, new_data = get2Dslice(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # # data_region, x, y = get_data_region(new_data, xaxis=x_new, yaxis=y_new, xbounds=(-18, 15.5))
    # # plot2D(x=x, y=y, data=data_region)
    # # gauss_data = gaussian_filter(data_region, sigma=2)
    # # plot2D(x=x, y=y, data=gauss_data)
    #
    # # export averaged data to hdf5 file
    # avg_data_hdf(fp, 'FS_gamma0', datat, p=p, slice_scale=ss, channel_scale=cs)
    #
    # # load averaged data from hdf5 file
    data, ss, cs, p = data_from_hdf(fp, 'FS_gamma0.h5')

    slice_val = 16.6

    # plot averaged data
    # plot3D(x=ss, y=cs, z=p, data=data, slice_dim='y', slice_val=slice_val, int_range=0.1,
    #        title=f'FS at gamma lamp {slice_val}eV')

    # k-correction + hexagons
    val = np.round(slice_val - 16.8, 1)
    kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='FS_gamma0.h5', val=val, Eint=0.1)
    fig = plot2D(x=kx, y=ky, data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'FS at gamma lamp {val}eV')

    # add hexagons
    coords = gen_polygon(6, radius=0.42, rotation=30, translation=(0, -0.5))
    new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(coords, radius=0.42,
                                                                                         rotation=30,
                                                                                         translation=(0, -0.5))
    plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r], fig=fig, color='black')

    """HOMO, Gamma0, XUV (April 2021)"""

    fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/XUV/3D/phi_motor_scan/XUV_HOMO_gamma0_2'
    fn = r'XUV_HOMO_gamma0_2_Y2021M04D28dT12h37m51s.h5'

    slice_val = 17.1  # XUV @ 19H, 22.4eV

    fns = ['XUV_HOMO_gamma0_2_Y2021M04D28dT12h37m51s.h5', 'XUV_HOMO_gamma0_2_Y2021M04D28dT13h01m38s.h5',
           'XUV_HOMO_gamma0_2_Y2021M04D28dT13h25m24s.h5', 'XUV_HOMO_gamma0_2_Y2021M04D28dT13h49m11s.h5']
    # get and average data from multiple hdf5 files
    data, ss, cs, p = data_from_hdf(fp, fn)
    data_avg = avg_array_from_hdfs(fp, fns)
    datat = np.swapaxes(data_avg, 2, 0)
    # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # x_new, y_new, new_data = get2Dslice(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # data_region, x, y = get_data_region(new_data, xaxis=x_new, yaxis=y_new, xbounds=(-18, 15.5))
    # plot2D(x=x, y=y, data=data_region)
    # gauss_data = gaussian_filter(data_region, sigma=2)
    # plot2D(x=x, y=y, data=gauss_data)

    # # export averaged data to hdf5 file
    avg_data_hdf(fp, 'XUV_HOMO_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)

    # load averaged data from hdf5 file
    data, ss, cs, p = data_from_hdf(fp, 'XUV_HOMO_averaged.h5')

    # plot averaged data
    plot3D(x=ss, y=cs, z=p, data=data, slice_dim='y', slice_val=slice_val, int_range=0.1,
           title=f'HOMO at gamma XUV {slice_val}eV')

    # k-correction + hexagons
    slice_val = 16.6
    EF = 18.3
    val = np.round(slice_val - EF, 1)
    kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='XUV_HOMO_averaged.h5', val=val, Eint=0.1, EF=EF)
    fig = plot2D(x=kx, y=ky, data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'HOMO at gamma XUV {val}eV')

    # add hexagons
    coords = gen_polygon(6, radius=1.1, rotation=30, translation=(0, 0.25))  # Au at -1.2 eV
    coords = gen_polygon(6, radius=0.42, rotation=30, translation=(0, -0.5))  # C60 at HOMO
    new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(coords, radius=0.42,
                                                                                         rotation=30,
                                                                                         translation=(0, -0.5))
    plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r], fig=fig, color='black')
    plot_polygon(coords, fig=fig, color='orange')

    """FS, Gamma0, XUV (April 2021)"""

    fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/XUV/3D/phi_motor_scan/XUV_FS_gamma0'
    fn = r'XUV_FS_gamma0_Y2021M04D28dT14h56m26s.h5'

    slice_val = 18.2  # XUV @ 19H, 22.4eV

    fns = ['XUV_FS_gamma0_Y2021M04D28dT14h56m26s.h5', 'XUV_FS_gamma0_Y2021M04D28dT15h54m10s.h5',
           'XUV_FS_gamma0_Y2021M04D28dT16h51m57s.h5']
    # get and average data from multiple hdf5 files
    data, ss, cs, p = data_from_hdf(fp, fn)
    data_avg = avg_array_from_hdfs(fp, fns)
    datat = np.swapaxes(data_avg, 2, 0)
    # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # x_new, y_new, new_data = get2Dslice(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # data_region, x, y = get_data_region(new_data, xaxis=x_new, yaxis=y_new, xbounds=(-18, 15.5))
    # plot2D(x=x, y=y, data=data_region)
    # gauss_data = gaussian_filter(data_region, sigma=2)
    # plot2D(x=x, y=y, data=gauss_data)

    # # export averaged data to hdf5 file
    avg_data_hdf(fp, 'XUV_FS_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)

    # load averaged data from hdf5 file
    data, ss, cs, p = data_from_hdf(fp, 'XUV_FS_averaged.h5')

    # plot averaged data
    EF = 18.4
    slice_val = 18.2
    val = np.round(slice_val - EF, 2)
    plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='y', slice_val=slice_val, int_range=0.1,
           title=f'FS at gamma XUV {slice_val}eV')

    # k-correction + hexagons
    slice_val = 18.0
    val = np.round(slice_val - 18.2, 1)
    kx, ky, kdata = kcorrect_phimotor(fp=fp, fn='XUV_FS_averaged.h5', val=val, Eint=0.1, EF=18.2)
    fig = plot2D(x=kx, y=ky, data=kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title=f'FS at gamma XUV {val}eV')

    # add hexagons
    coords = gen_polygon(6, radius=1.1, rotation=30, translation=(0, 0.25))  # Au at -1.2 eV
    coords = gen_polygon(6, radius=0.42, rotation=30, translation=(0, -0.5))  # C60 at HOMO
    new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(coords, radius=0.42,
                                                                                         rotation=30,
                                                                                         translation=(0, -0.5))
    plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r], fig=fig, color='black')
    plot_polygon(coords, fig=fig, color='orange')

    """FS2, Gamma0, Lamp (April 2021)"""

    fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/FermiEdge_2'
    fn = r'Lamp_FermiEdge2_Y2021M04D28dT22h33m17s.h5'

    slice_val = 16.2

    fns = ['Lamp_FermiEdge2_Y2021M04D28dT22h33m17s.h5', 'Lamp_FermiEdge2_Y2021M04D28dT23h31m01s.h5',
           'Lamp_FermiEdge2_Y2021M04D29dT00h28m47s.h5', 'Lamp_FermiEdge2_Y2021M04D29dT01h26m32s.h5',
           'Lamp_FermiEdge2_Y2021M04D29dT02h24m17s.h5', 'Lamp_FermiEdge2_Y2021M04D29dT03h22m02s.h5',
           'Lamp_FermiEdge2_Y2021M04D29dT04h19m48s.h5', 'Lamp_FermiEdge2_Y2021M04D29dT05h17m33s.h5'
           ]
    # get and average data from multiple hdf5 files
    data, ss, cs, p = data_from_hdf(fp, fn)
    data_avg = avg_array_from_hdfs(fp, fns)
    datat = np.swapaxes(data_avg, 2, 0)
    # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # x_new, y_new, new_data = get2Dslice(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=slice_val, int_range=0.1)
    # data_region, x, y = get_data_region(new_data, xaxis=x_new, yaxis=y_new, xbounds=(-18, 15.5))
    # plot2D(x=x, y=y, data=data_region)
    # gauss_data = gaussian_filter(data_region, sigma=2)
    # plot2D(x=x, y=y, data=gauss_data)

    # # export averaged data to hdf5 file
    avg_data_hdf(fp, 'FS2_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)

    # load averaged data from hdf5 file
    data, ss, cs, p = data_from_hdf(fp, 'FS2_averaged.h5')

    # plot averaged data
    EF = 16.8
    slice_val = 16.6
    val = np.round(slice_val - EF, 2)
    plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='y', slice_val=slice_val, int_range=0.1,
           title=f'FS2 at gamma lamp {slice_val}eV')
    

    """Combination Plots!!!"""

    # load averaged data from hdf5 file
    fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/XUV/3D/phi_motor_scan/XUV_FS_gamma0'
    data, ss, cs, p = data_from_hdf(fp, 'XUV_FS_averaged.h5')

    fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/FermiEdge_gamma0'
    data, ss, cs, p = data_from_hdf(fp, 'FS_gamma0.h5')

    fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/XUV/3D/phi_motor_scan/XUV_HOMO_gamma0_2'
    data, ss, cs, p = data_from_hdf(fp, 'XUV_HOMO_averaged.h5')

    fp = r'/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/phi_motor_scan_overnight'
    data, ss, cs, p = data_from_hdf(fp, 'HOMO_overnight_averaged.h5')

    # plot averaged data
    EF = 18.4
    slice_val = 17.7
    val = np.round(slice_val - EF, 2)
    plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='y', slice_val=slice_val, int_range=0.5, xlabel='Theta',
           ylabel='Phi', title=f'{val}eV (XUV)')

    EF = 16.8
    slice_val = 16.5
    val = np.round(slice_val - EF, 2)
    plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='y', slice_val=slice_val, int_range=0.5, xlabel='Theta',
           ylabel='Phi', title=f'{val}eV (Lamp)')

    # from Typing import List
    # heatmap_info is list of lists: fp, fn, slice_val, EF
    def getheatmaps(heatmap_info: list, int_range=0.5, color='Plasma'):
        heatmaps = []
        for h in heatmap_info:
            data, ss, cs, p = data_from_hdf(h[0], h[1])
            fig = plot3D(x=ss, y=cs, z=np.flip(p), data=data, slice_dim='y', slice_val=h[3], int_range=int_range,
                         show=False, color=color)
            hm = fig.data[0]
            heatmaps.append(hm)
        return heatmaps


    def multiplot_3D(heatmap_info: list, rows, cols, int_range=0.5, xlabel=None, ylabel=None, title=None, color='Plasma'):
        fig = make_subplots(rows=rows, cols=cols, shared_xaxes=True, shared_yaxes=True)
        coords = list(product(range(1, rows + 1), range(1, cols + 1)))
        heatmaps = getheatmaps(heatmap_info=heatmap_info, int_range=int_range, color=color)
        for coords, heatmap in zip(coords, heatmaps):
            fig.add_trace(heatmap, row=coords[0], col=coords[1])

    hm_info = ['/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/XUV/3D/phi_motor_scan/XUV_HOMO_gamma0_2',
               'XUV_HOMO_averaged.h5', 16.9, 18.2], \
              ['/Users/alexandratully/Desktop/ARPES Data/April_2021/LT/Lamp/3D/phi_motor_scan/phi_motor_scan_overnight']\
        , []




    # def multiplot_2D()


