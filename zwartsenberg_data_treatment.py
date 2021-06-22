#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 10:57:25 2018

@author: berend
"""

### data treatment functions
from scipy.signal import savgol_filter,convolve2d,fftconvolve
import numpy as np


def Gauss(x, A, x0, sig):
    """Gaussian function, sig is the FWHM"""
    Sig = sig/(2*np.sqrt(2*np.log(2)))
    return A * np.exp(-(x-x0)**2 /(2*Sig*Sig))
    

def Gauss2D(x,y, A, x0,y0, sigx,sigy):
    """Gaussian function, sig is the FWHM"""
    Sigx = sigx/(2*np.sqrt(2*np.log(2)))
    Sigy = sigy/(2*np.sqrt(2*np.log(2)))
    return A * np.exp(-(x-x0)**2 /(2*Sigx*Sigx)-(y-y0)**2 /(2*Sigy*Sigy))
    
        


def smooth(data,x,y, dx = None, dy = None, axis = 'x', kernel = 'savgol', **kwargs):
    """Smooth data with respect to an axis
    Args:
        data: 2D data input
        x: xaxis (row scaling)
        y: yaxis (column scaling)
        dx: xaxis filter window
        dy: yaxis filter window
        axis: axis to be smoothed along"""
        
        
    if kernel == 'savgol':
        
        porder = kwargs['porder'] if 'porder' in kwargs else 2
        
                
        if axis == 'x':
            if dx is not None:
                filt_size = 2*int(abs(0.5*dx/(x[1]-x[0])))+1
                
            else:
                filt_size = 2*int(0.05*x.shape[0])+1

            smth = np.zeros_like(data)
            for i in range(smth.shape[1]):
                smth[:,i] = savgol_filter(data[:,i],filt_size,porder,**kwargs)

        elif axis == 'y':
            if dy is not None:
                filt_size = 2*int(abs(0.5*dy/(y[1]-y[0])))+1
                
            else:
                filt_size = 2*int(0.05*y.shape[0])+1
            smth = np.zeros_like(data)
            for i in range(smth.shape[1]):
                smth[i] = savgol_filter(data[i],filt_size,porder,**kwargs)
                
        return smth
        
    if kernel == 'gauss':
        
        if axis == 'x':
        
            filt_axis = np.arange(-4*dx, +4*dx + x[1]-x[0], x[1]-x[0])
            sig = dx
            rs = (-1,1)
            
        elif axis == 'y':

            filt_axis = np.arange(-4*dy, +4*dy + y[1]-y[0], y[1]-y[0])
            sig = dy
            rs = (1,-1)
            
            
        filt = Gauss(filt_axis, 1.0, 0.0, sig)/np.sum(Gauss(filt_axis, 1.0, 0.0, sig))
        filt = filt.reshape(rs)
        print(filt.shape)
        
        return fftconvolve(pad_data(data,filt.shape),filt, mode = 'valid')
        
    if kernel == 'gauss2d':
        
        filt_xaxis = np.arange(-4*dx, +4*dx + x[1]-x[0], x[1]-x[0])
        filt_yaxis = np.arange(-4*dy, +4*dy + y[1]-y[0], y[1]-y[0])
        FX,FY = np.meshgrid(filt_xaxis,filt_yaxis)
        filt = Gauss2D(FX,FY, 1.0, 0.0,0.0, dx,dy)/np.sum(Gauss2D(FX,FY, 1.0, 0.0,0.0, dx,dy))
        print(filt.shape)

        return fftconvolve(pad_data(data,filt.shape),filt, mode = 'valid')
        
    
def second_d(data,x,y, mode = 'x'):
    """Make second derivative of data with scaling x and y
    
    Args:
        data: 2D image data to take the 2nd derivative on
        x: xaxis of 2D image (rows)
        y: yaxis of 2D image (columns)
        mode: x,y or xy, indicating x 2nd derivative, y 2nd derivative or laplacian"""
        
        
        
    filt_dic = {'x' : np.array([1,-2,1]).reshape((-1,1)),
                'y' : np.array([1,-2,1]).reshape((1,-1)),
                'xy' : np.array([[0,1,0],[1,-4,1],[0,1,0]])}


    return convolve2d(data, filt_dic[mode],mode = 'same',boundary = 'symm')
   
    
def pad_data(data, filtshape, mode = 'edge'):
    """Pad data for convolution using filtshape
    Args:
        data: 2d data array to be padded
        filtshape: tuple describing the filter to be applied to data
    Returns
        padded array (same)"""
        
    odd = [(filt+1)%2 for filt in filtshape]
                           
    pad_width = tuple([(int((filtlen-1)/2),int((filtlen-1)/2)+o) for filtlen,o in zip(filtshape,odd)])

    return np.pad(data,pad_width,mode)
    

    
        
        
    