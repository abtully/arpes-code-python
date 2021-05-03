import h5py
import os
import numpy as np


def array_from_hdf(fp_complete: str, dataname: str) -> np.ndarray:
    with h5py.File(fp_complete, 'r') as hdf:  # makes sure hdf file is closed at the end
        dataset = hdf.get(dataname)
        data = dataset[:]
    return data


# def array_from_hdf(fp_complete: str, dataname) -> np.ndarray:
#     with h5py.File(fp_complete, 'r') as hdf:  # makes sure hdf file is closed at the end
#         dataset = hdf.get(d for d in dataname)
#         datas = dataset[:]
#     return datas


def avg_array_from_hdfs(fp: str, fns: list) -> np.ndarray:
    datas = []
    for f in fns:
        filename = os.path.join(fp, f)
        data = array_from_hdf(filename, 'data')
        datas.append(data)
    return np.mean(datas, axis=0)


def data_from_hdf(fp: str, dataname: str):  # returns 4 ndarrays
    p = array_from_hdf(os.path.join(fp, dataname), 'p')
    ss = array_from_hdf(os.path.join(fp, dataname), 'slice_scale')
    cs = array_from_hdf(os.path.join(fp, dataname), 'channel_scale')
    data = array_from_hdf(os.path.join(fp, dataname), 'data')
    return data, ss, cs, p


def avg_data_hdf(fp: str, fn: str, data_avg: np.ndarray, p: np.ndarray, slice_scale: np.ndarray,
                 channel_scale: np.ndarray):
    filepath = os.path.join(fp, f'{fn}.h5')
    with h5py.File(filepath, 'w') as hdf:
        hdf['data'] = data_avg
        hdf['p'] = p
        hdf['slice_scale'] = slice_scale
        hdf['channel_scale'] = channel_scale


if __name__ == '__main__':
    from ali_plotting_functions import plot3D, plot2D
    from ali_analysis_functions import get_data_region
    from ali_k_correction import get2Dslice
    from scipy.ndimage import gaussian_filter
    from ali_k_correction import kcorrect_phimotor, get2Dslice
    from ali_polygons import gen_polygon, gen_tiled_hexagons, plot_polygons, plot_polygon

    """JANUARY 2021"""

    # #
    # # """HOMO, Theta=15"""
    # #
    # # fp = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\January_2021\LT\Lamp\3D\phi_motor_scan\HOMO15_Y2021M01D27'
    # # fn = r'HOMO15_Y2021M01D27dT21h18m41s.h5'
    # # fns = ['HOMO15_Y2021M01D27dT21h18m41s.h5', 'HOMO15_Y2021M01D28dT00h08m04s.h5', 'HOMO15_Y2021M01D28dT02h57m28s.h5',
    # #        'HOMO15_Y2021M01D28dT05h46m53s.h5', 'HOMO15_Y2021M01D28dT08h36m19s.h5']
    # #
    # # # get and average data from multiple hdf5 files
    # # data, ss, cs, p = data_from_hdf(fp, fn)
    # # data_avg = avg_array_from_hdfs(fp, fns)
    # # datat = np.swapaxes(data_avg, 2, 0)
    # # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=15.7)
    # #
    # # # export averaged data to hdf5 file
    # # avg_data_hdf(fp, 'HOMO15_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)
    # #
    # # # load averaged data from hdf5 file
    # # data, ss, cs, p = data_from_hdf(fp, 'HOMO15_averaged.h5')
    # # # plot averaged data
    # # plot3D(x=ss, y=cs, z=p, data=data, slice_dim='y', slice_val=15.7)
    # #
    # # """HOMO, Theta=0"""
    # #
    # # fp = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\January_2021\LT\Lamp\3D\phi_motor_scan\HOMO0_Y2021M01D27'
    # # fn = r'HOMO0_Y2021M01D27dT18h17m03s.h5'
    # #
    # # # load data
    # # data, ss, cs, p = data_from_hdf(fp, fn)
    # # # fix data
    # # datat = np.swapaxes(data, 2, 0)
    # # test = np.flip(datat, 0)
    # # # plot data
    # # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=15.7)
    # # plot3D(x=ss, y=cs, z=p, data=test, slice_dim='y', slice_val=15.7)
    #
    # """FS, Theta=?"""
    # from ali_plotting_functions import plot3D
    #
    # fp = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\January_2021\LT\Lamp\3D\phi_motor_scan\FS_Y2021M01D25'
    # # fn = r'FS_Y2021M01D25dT22h39m08s.h5'
    # # fns = ['FS_Y2021M01D25dT22h39m08s.h5', 'FS_Y2021M01D26dT01h28m27s.h5', 'FS_Y2021M01D26dT04h17m48s.h5']
    # #
    # # # get and average data from multiple hdf5 files
    # # data, ss, cs, p = data_from_hdf(fp, fn)
    # # data_avg = avg_array_from_hdfs(fp, fns)
    # # datat = np.swapaxes(data_avg, 2, 0)
    # # plot3D(x=ss, y=cs, z=p, data=datat, slice_dim='y', slice_val=15.7)
    # #
    # # # export averaged data to hdf5 file
    # # avg_data_hdf(fp, 'FS_averaged', datat, p=p, slice_scale=ss, channel_scale=cs)
    #
    # # load averaged data from hdf5 file
    # data, ss, cs, p = data_from_hdf(fp, 'FS_averaged.h5')
    # # plot averaged data
    # fig = plot3D(x=ss, y=cs, z=p, data=data, slice_dim='y', slice_val=15.7)
    #
    # """K-corrected FS, Theta=?"""
    # # from ali_k_correction import fix_EkatEF, kcorrect3D
    # #
    # # # correct EF
    # # fix_EkatEF(data, 16.8)
    # # # generate and plot k corrected data
    # # kx, ky, kdata = kcorrect3D(self=d, val=15)  # mine
    # # plot2D(kx, ky, kdata, xlabel='kx [A-1]', ylabel='ky [A-1]', title='Constant Energy Slice')  # mine
    #
    # """Add Hexagons"""
    # from ali_polygons import gen_polygon, gen_tiled_hexagons, plot_polygons
    # coords = gen_polygon(6, radius=0.42, rotation=30, translation=(-0.07, 0.1))
    # new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r = gen_tiled_hexagons(coords, radius=0.42,
    #                                                                                      rotation=30,
    #                                                                                      translation=(-0.07, 0.1))
    # plot_polygons([coords, new_coords, coords_bl, coords_tr, coords_l, coords_tl, coords_r], fig=fig)
    #
    # # # datat = np.moveaxis(data, 2, 0)
    # # datat = np.swapaxes(data_avg, 2, 0)
    #
    # # hdf.keys()  # shows groups and datasets
    # # hdf.attrs.keys()  # shows attributes
    # # hdf['channel_scale'].attrs.keys()
    # # cs = hdf.get('channel_scale')  # nice way of getting things from my hdf file
    # # cs  # cs is a dataset, if close hdf file this dataset will close as well (can't do anything with it)
    # # cs_data = cs[:]  # numpy array, this will stay open after I close hdf
    # # cs_data  # energy scale
    # # cs_data.shape  # (1064,)
    # # hdf.close()
