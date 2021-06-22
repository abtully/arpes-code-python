"""
Testing Imports of Ali's Data with Berend's Code
"""

import lmfit
import sys
import PIL
import scipy
import h5py
import numpy as np
import igor
import matplotlib.pyplot as plt
import matplotlib as mpl
# mpl.use('tkagg')

Aapath = 'C:/Users/atully/Code/ARPES Code Python'
sys.path.append(Aapath)

from zwartsenberg_Data2DObject import AaData2D, AaData3D, UBCGold

datapath = 'C:/Users/atully/Code/ARPES Code Python/data/'

#load 3D data
f06 = AaData3D(datapath + 'October_2020/RT/3D/FS_HOMO0001', datatype = 'APE', zaxis='theta_y')

#slice through x, y, z
fig, ax = f06.show(mode = 'CE', val = 15, Eint = 0.02)  # use f0.xaxis to
fig.show()
f06.show(mode = 'row', val = 15)[0].show()  # same as above but with default = no intergration
f06.show(mode = 'lay', val = -10)  # takes a constant layer (z axis value) phi = -10
f06.show(mode = 'col', val = 0)  # takes a constant y value

#paramters of every object are in obj.p
print(f06.p)

#load gold
gold = UBCGold(datapath + 'October_2020/RT/HS_GM0001_001.ibw', datatype = 'igor')
#gold won't correct without this paramter
gold.p['slitcorrected'] = True

#you can look at the gold correction data like you can with other objects
gold.show()[0].show()

#to find region, plot across gold.data
print(gold.shape)
fig, ax = plt.subplots(1)
ax.plot(gold.data[:,500])  # --> 800 - 1000
fig.show()
plt.plot(gold.data[800,:]) # --> 120 - 840

# region = (x_start, x_end, y_start, y_end)
gold.Analyze(T = 150, region = (800, 1000, 120, 840))
# gold.Analyze(T = 150, region = (16.5, 17.5, -15, 15))  # select the region where you want to fit the data

#show the generated fit
gold.showgold()

#use the gold file to correct objects
f06.Correct(goldobj = gold)  # applies gold correction


f06.data
f06.data.shape
plt.imshow(f06.data[500,:,:])
print(f06.data.shape)
plt.plot(f06.data[:,500,60])
data_bin = f06.data.reshape((-1, 500, 2, 121)).mean(axis=2)  # way of binning data
plt.imshow(data_bin[:,:,60])
f06.xaxis  # energies
f06.yaxis  # y angle that you see on screen while taking a spectrum
f06.zaxis  # z different slices that you're taking sequentially