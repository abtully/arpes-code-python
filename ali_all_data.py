"""
Code to Reshape Multi-Scan FS Data to Work with Loading and Plotting Functions
@author: Alexandra Tully
@date: November 2020

Note:
This rounds down. So try to load the data, then divide array size by (45*1000*1064); seems to
usually be 2.689 -- this code creates a blank array of the right size, and then only uses the amount of
data that fills that array (= 2 full rounds of data).
"""

import numpy as np
import os

# 140K, Fermi_Surface
# path = r'C:\Users\atully\Code\ARPES Code Python\data\October_2020\140K\3D\FS_Fermi_Surface_140K0001'
# binfile = 'Spectrum_SK_Lamp_FS3.bin'

# 9K, Fermi_Edge
# path = r'C:\Users\atully\Code\ARPES Code Python\data\October_2020\9K\3D\FS_Fermi_Edge_9K0001'
# binfile = 'Spectrum_SK_Lamp_FS4.bin'

# RT, Fermi_Edge (1, 2, 3)
path = r'C:\Users\atully\Code\ARPES Code Python\data\October_2020\RT\3D\FS_Fermi_Edge0001'
binfile = 'Spectrum_SK_Lamp_FS0.bin'

layers, columns, rows = 45, 1000, 1064
total = layers*columns*rows*2
data = np.fromfile(path + '/' + binfile, dtype='float32')[:total]
data = np.mean(data.reshape((2, layers, columns, rows)), axis=0)
data.tofile(os.path.join(path, 'reshaped.bin'))

