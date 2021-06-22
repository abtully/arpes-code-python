"""
HDF5 file loader (required for phi_motor ARPES data)
@author: Alexandra Tully
@date: April 2021
"""

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
    path = r'C:\Users\atully\Code\ARPES Code Python\analysis_data\January_2021\LT\Lamp\3D\phi_motor_scan' \
           r'\HOMO15_Y2021M01D27\HOMO15_Y2021M01D28dT08h36m19s.h5 '
    hdf = h5py.File(path, 'r')
    hdf.keys()  # shows groups and datasets
    hdf.attrs.keys()  # shows attributes
    hdf['channel_scale'].attrs.keys()
    cs = hdf.get('channel_scale')  # nice way of getting things from my hdf file
    cs  # cs is a dataset, if close hdf file this dataset will close as well (can't do anything with it)
    cs_data = cs[:]  # numpy array, this will stay open after I close hdf
    cs_data  # energy scale
    cs_data.shape  # (1064,)
    hdf.close()
