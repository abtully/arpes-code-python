#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 10:43:23 2019

@author: berend
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

## Example
#python 3.x
#h5py
#igor --> pip install igor
#numpy
#matplotlib
#pillow
#scipy
#lmfit --> pip install lmfit


#append path to Aa files to system path
#this is the path the Aa files are located, change as necessary
# Aapath = 'C:/Users/atully/Desktop'
Aapath = 'C:/Users/atully/Code/ARPES Code Python'
sys.path.append(Aapath)


#import objects
from zwartsenberg_Data2DObject import AaData2D, AaData3D, UBCGold


#variable to keep data path
# datapath = 'C:/Users/atully/Desktop/data/'
datapath = 'C:/Users/atully/Code/ARPES Code Python/data/'

#load 3D data
f06 = AaData3D(datapath + 'FS_01/FS', datatype = 'CLS')

#slice through x, y, z
fig, ax = f06.show(mode = 'CE', val = 43.68, Eint = 0.02)
fig.show()
f06.show(mode = 'row', val = 43.68)[0].show()
f06.show(mode = 'lay', val = -10)
f06.show(mode = 'col', val = 0)

#paramters of every object are in obj.p
print(f06.p)

#load gold
gold = UBCGold(datapath + 'f_0002gold.ibw', datatype = 'igor')
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

#show the generated fit
gold.showgold()


#use the gold file to correct objects
f06.Correct(goldobj = gold)

#now the default plotting behavior changed
f06.show(val = 0.0, Eint = 0.01)


# slitline = (theta_x, theta_y, width)
f06.show(val = 0.0, Eint = 0.01, slitline = (12, -13, 30))


#use it to correct angles (theta_x, theta_y, rot)
f06.setangles(2, -2, -2)



#path through k-space
f06.show(mode = 'k2', points = 200, path = [[-0.8,0],[0.,0.],[0.8,-0.8]])



#similar for 2D datamaps
spec = AaData2D(datapath + 'f_0001SRO_spec.ibw', datatype = 'igor')
spec.show()
spec.Correct(goldobj = gold)

#use it for custom plotting, obj.data has the numpy array that holds the actual data
plt.plot(spec.xaxis[500:], spec.data[500:,300:400].mean(axis = 1))

#save and load
f06.save(path = datapath + 'f06.h5')
f03 = AaData3D(datapath + 'f06.h5')


#use numpy functions for normalizations etc:
f06.data /= f06.data.max()




