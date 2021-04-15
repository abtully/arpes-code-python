# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 14:25:55 2017

@author: Berend
"""

import numpy as np
import matplotlib.pyplot as plt
# from arpes.kcorrection import kfromangles, kzfromhv
from kcorrection_new_new import kfromangles, kzfromhv


################## 1D plotters ###############################################


def plot1D(Aa1D, color='black', ax=None, plotrange=None, size=3.0, asp=3/4, plotkwargs={}):
    """Plot a 1D dataset


    Args:
        color (str): line color
        ax (axis): axis object to plot in
        plotrange (list): set the plotrange
        size (float): size of the generated figure
        asp (float): aspect ratio, when passed as 'true', calculate based on kx/kz range
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """    
    
    if plotrange is not None:
        xstart,xend,ystart,yend = plotrange
    else:
        ystart = np.min(Aa1D.data)
        yend = 1.2*np.max(Aa1D.data)
        xstart = Aa1D.xaxis[0]
        xend = Aa1D.xaxis[-1]


    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))
    
    ax.plot(Aa1D.xaxis, Aa1D.data, **plotkwargs, color = color)
    ax.set_ylabel('Counts [arb. un.]', fontsize = 12)
    ax.set_xlabel('Energy [eV]', fontsize = 12)  
    ax.axis([xstart, xend, ystart, yend])

    return ax

    
    
    
    
    
    
    
    
    
################## 2D plotters ###############################################
    


        
def plot2D(Aa2D, vmin=None, vmax=None, cmap='Spectral_r', ax=None, plotrange=None, size=3.0, asp=3/4, plotkwargs={}):
    """plot a 2D dataset

    

    Args:
        vmin (float): start of the colormap
        vmax (float): end of the colormap
        cmap (str): matplotlib colormap
        ax (axis): axis object to plot in
        plotrange (list): set the plotrange
        size (float): size of the generated figure
        asp (float): aspect ratio
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """

    
    if vmin is None:
        vmin = Aa2D.data.min()
    if vmax is None:
        vmax = Aa2D.data.max()
    
    if plotrange is not None:
        xstart,xend,ystart,yend = plotrange
    else:
        ystart = Aa2D.ystart
        yend = Aa2D.ystart+(Aa2D.data.shape[1]-1)*Aa2D.ydelta
        xstart = Aa2D.xstart
        xend = Aa2D.xstart+(Aa2D.data.shape[0]-1)*Aa2D.xdelta
    ratio = abs((yend-ystart)/(xend-xstart)*asp)
    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))
    
    ax.imshow(np.flipud(Aa2D.data), cmap=cmap, interpolation = 'none', extent = [ystart,yend,xstart,xend], aspect=ratio, vmin=vmin, vmax=vmax, **plotkwargs)
    ax.set_xlabel('Theta [deg]', fontsize = 12)
    ax.set_ylabel('Energy [eV]', fontsize = 12)    
    return ax
        
def kplot2D(Aa2D, kmode='edep', vmin=None, vmax=None, cmap='Spectral_r', ax=None, 
            plotrange=None, size=3.0, asp=3/4, plotkwargs={}, rasterized=True):
    """plot a 2D dataset with k-correction applied

    

    Args:
        kmode (string): 'edep' (dependent on energy) or 'simple' (taken at 0 energy). 
                    Always use 'edep' in normal operation
        vmin (float): start of the colormap
        vmax (float): end of the colormap
        cmap (str): matplotlib colormap
        ax (axis): axis object to plot in
        plotrange (list): set the plotrange
        size (float): size of the generated figure
        asp (float): aspect ratio
        plotkwargs (dict): additional arguments to be passed to the plot function
        rasterized (bool): plot the mesh rasterized, recommended for saving pdfs
    Returns:
        ax (axis): axis object that was plotted in

    """


    if Aa2D.p['type'] == 'rawdata':
        raise RuntimeError('k plots only work on corrected data')

    if vmin is None:
        vmin = Aa2D.data.min()
    if vmax is None:
        vmax = Aa2D.data.max()

        
    if(kmode == 'simple'):
        kx,ky = kfromangles(Aa2D.yaxis, 0., Aa2D.p['EkatEF'], theta_m = Aa2D.p['theta_m'], phi_m = Aa2D.p['phi_m'], theta0 = Aa2D.p['theta0'], phi0 = Aa2D.p['phi0'], alpha0 = Aa2D.p['alpha0']) ## phi is set to zero, because manipulator offset is already captured by phi_m
        dkx = np.zeros_like(kx)
        dky = np.zeros_like(ky)
        dkx[1:] = kx[1:]-kx[:-1]
        dky[1:] = ky[1:]-ky[:-1]
        dk = np.sqrt(dkx*dkx + dky*dky)
        kax = np.cumsum(dk)
        argzero = np.argmin(np.sqrt(kx**2 + ky**2))
        kax -= kax[argzero]
        kmesh,Emesh = np.meshgrid(kax,Aa2D.xaxis)
    elif(kmode == 'edep'):
        thmesh,Emesh = np.meshgrid(Aa2D.yaxis,Aa2D.xaxis+Aa2D.p['EkatEF'])
        kx,ky = kfromangles(thmesh, 0., Emesh, theta_m = Aa2D.p['theta_m'], phi_m = Aa2D.p['phi_m'], theta0 = Aa2D.p['theta0'], phi0 = Aa2D.p['phi0'], alpha0 = Aa2D.p['alpha0'])
        dkx = np.zeros_like(kx).astype('float64') #this is to prevent a bug in np that raises an error taking the sqrt of 0
        dky = np.zeros_like(ky).astype('float64')
        dkx[:,1:] = kx[:,1:]-kx[:,:-1]
        dky[:,1:] = ky[:,1:]-ky[:,:-1]        

        dk = np.sqrt(dkx*dkx + dky*dky)
        kmesh = np.cumsum(dk, axis = 1)
        
        argzero = np.argmin(np.sqrt(kx**2 + ky**2), axis = 1)
        
        for i in range(kmesh.shape[0]):
            kmesh[i] -= kmesh[i,argzero[i]]
        
        Emesh -= Aa2D.p['EkatEF']
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))
    
    ax.pcolormesh(kmesh, Emesh, Aa2D.data, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=rasterized, **plotkwargs)
    ax.set_xlabel('kx [A-1]', fontsize = 12)
    ax.set_ylabel('Energy [eV]', fontsize = 12)    
    if plotrange is None:
        
        ax.axis([np.min(kmesh),np.max(kmesh),np.min(Aa2D.xaxis),np.max(Aa2D.xaxis)])
    else:
        ax.axis(plotrange)


    return ax





def plotEDC2D(Aa2D,path, pathtype='k', NEDC=10, avg=1, offset=0.3, Estart=None, Eend=None, color='black',
              ax=None, plotrange=None, size=3.0, asp=3/4, plotkwargs={}):
    """plot a range of EDCs along a path

    Args:
        path (list): defines the kpath, can be [kx0,ky0,kx1,ky1] with the path
                     going from k0 to k1, or a list containing arrays [kx,ky] with the full path
                     Pass k if pathtype 'k', or angle if pathtype 'a'
        pathtype (str): 'k' or 'a' type of the path (see AaData2D.getEDCs for info)
        NEDC (int): number of EDCs
        avg (int): number of EDCs to average
        offset (float): offset between curves in the EDC stack
        Estart (float): energy range lower limit
        Eend (float): energy range upper limit
        color (str): color of k-waypoint indicator lines
        ax (axis): axis object to plot in
        plotrange (list): set the plotrange        
        size (float): size of the generated figure
        asp (float): aspect ratio
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """
    
    if Estart is None:
        Estart = Aa2D.xaxis[0]
    if Eend is None:
        Eend = Aa2D.xaxis[-1]
                                    

    EDCs, Escale, path = Aa2D.getEDCs(path=path, pathtype=pathtype, NEDC=NEDC, avg=avg, Estart=Estart, Eend=Eend)
    
    avg_offset = np.mean(np.array([np.mean(EDC[0:10]) for EDC in EDCs]))
        
    doffset = offset*avg_offset

    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))

    


    for i in range(NEDC):
        ax.plot(Escale,(EDCs[i]+i*doffset), *plotkwargs)
       
            
    ax.set_xlabel('Energy [eV]', fontsize = 12)
    ax.set_ylabel('Photocurrent [arb. un.]', fontsize = 12) 
    ax.set_xlim(Estart,Eend)      
    
    return ax
    
    

    
def plotMDC2D(Aa2D):
    """Plot a stack of MDCs. Not implemented"""
    raise NotImplementedError('2D MDCs not yet implemented')
  



################## 3D plotters ###############################################





def plotCE3D(Aa3D, num=None, val=None, cmap='Spectral_r', vmin=None, vmax=None, Eint_n=1, Eint=None,
             ax=None, plotrange=None, size=3.0, asp=3/4, plotkwargs={}):
    """plot a constant energy map at num or val

    

    Args:
        num (int): the xaxis index at which energy the slice is taken
        val (float): energy at which the slice is taken, overrides num
        cmap (str): matplotlib colormap
        vmin (float): start of the colormap
        vmax (float): end of the colormap
        Eint_n (int): number of slices to average over
        Eint (float): energy range to average over, overrides Eint_n
        ax (axis): axis object to plot in
        plotrange (list): set the plotrange
        size (float): size of the generated figure
        asp (float): aspect ratio, when passed as 'true', calculate based on kx/ky range
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """
    
    
    if val is not None:
        num = int(round((val- Aa3D.xstart)/Aa3D.xdelta))
    elif num is None:
        raise RuntimeError('Provide either num or val')

    if((num < 0) or (num >= Aa3D.data.shape[0])):
        raise ValueError('Val or num outside of plot range')   
    
        

    if vmin is None:
        vmin = Aa3D.data[num].min()
    if vmax is None:
        vmax = Aa3D.data[num].max()



    if Eint is not None:
        Eint_n = int(Eint/Aa3D.xdelta)

    
    
    pyaxis = np.linspace(Aa3D.yaxis[0]-0.5*Aa3D.ydelta, Aa3D.yaxis[-1]+0.5*Aa3D.ydelta, (Aa3D.yaxis.shape[0]+1))
    pzaxis = np.linspace(Aa3D.zaxis[0]-0.5*Aa3D.zdelta, Aa3D.zaxis[-1]+0.5*Aa3D.zdelta, (Aa3D.zaxis.shape[0]+1))
    
    if asp == 'true':
        asp = (pzaxis[-1]-pzaxis[0])/(pyaxis[-1] - pyaxis[0])
    elif asp == 'ptrue':
        asp = float(Aa3D.yaxis.size)/(Aa3D.zaxis.size)

    
    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))
 
    meshy, meshz = np.meshgrid(pyaxis, pzaxis)
    ax.pcolormesh(meshy, meshz, np.mean(Aa3D.data[num-Eint_n:num+Eint_n], axis = 0).T, 
                  cmap=cmap, vmin=vmin, vmax=vmax,**plotkwargs)
    ax.set_xlabel('Theta [deg]', fontsize = 12)
    ax.set_ylabel('Phi [deg]', fontsize = 12)                


    if plotrange is None:
        ax.axis([pyaxis[0],pyaxis[-1], pzaxis[0],pzaxis[-1]])
    else:
        ax.axis(plotrange)


    return ax
    



        
    
def plotcol3D(Aa3D,num=None, val=None, cmap='Spectral_r', vmin=None, vmax=None,
             ax=None, plotrange=None, size=3.0, asp=3/4, plotkwargs={}):
    """plot an energy (xaxis) vs zaxis slice (constant yaxis)



    Args:
        num (int): the yaxis index at which energy the slice is taken
        val (float): yaxis value at which the slice is taken, overrides num
        cmap (str): matplotlib colormap
        vmin (float): start of the colormap
        vmax (float): end of the colormap
        ax (axis): axis object to plot in
        plotrange (list): set the plotrange
        size (float): size of the generated figure
        asp (float): aspect ratio
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """
    if val is not None:
        num = int(round((val- Aa3D.ystart)/Aa3D.ydelta))
    elif num is None:
        raise RuntimeError('Provide either num or val')

    if((num < 0) or (num >= Aa3D.data.shape[0])):
        raise ValueError('Val or num outside of plot range')   
    
        

    if vmin is None:
        vmin = Aa3D.data[num].min()
    if vmax is None:
        vmax = Aa3D.data[num].max()
    
    pxaxis = np.linspace(Aa3D.xaxis[0]-0.5*Aa3D.xdelta, Aa3D.xaxis[-1]+0.5*Aa3D.xdelta, (Aa3D.xaxis.shape[0]+1))
    pzaxis = np.linspace(Aa3D.zaxis[0]-0.5*Aa3D.zdelta, Aa3D.zaxis[-1]+0.5*Aa3D.zdelta, (Aa3D.zaxis.shape[0]+1))
    
    if asp == 'true':
        asp = abs((pxaxis[-1]-pxaxis[0])/pzaxis[-1]-pzaxis[0])*(40/1.2)*(1376/1024)
    elif asp == 'ptrue':
        asp = float(Aa3D.xaxis.size)/(Aa3D.zaxis.size)

    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))

    ax.set_xlabel('zaxis', fontsize = 12)
    ax.set_ylabel('Energy [eV]', fontsize = 12)         
    meshz, meshx = np.meshgrid(pzaxis, pxaxis)
    ax.pcolormesh(meshz, meshx, Aa3D.data[:,num,:], cmap=cmap, vmin=vmin, vmax=vmax, **plotkwargs)
    
    if plotrange is None:
        ax.axis([pzaxis[0],pzaxis[-1], pxaxis[0],pxaxis[-1]])
    else:
        ax.axis(plotrange)
        
    return ax
  
    
    
def plotlay3D(Aa3D,num=None, val=None, cmap='Spectral_r', vmin=None, vmax=None,
             ax=None, plotrange=None, size=3.0, asp=3/4, plotkwargs={}):
    """plot an energy (xaxis) vs yaxis slice (constant zaxis)



    Args:
        num (int): the zaxis index at which energy the slice is taken
        val (float): zaxis value at which the slice is taken, overrides num
        cmap (str): matplotlib colormap
        vmin (float): start of the colormap
        vmax (float): end of the colormap
        ax (axis): axis object to plot in
        plotrange (list): set the plotrange
        size (float): size of the generated figure
        asp (float): aspect ratio
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """
              
    
    if val is not None:
        num = int(round((val- Aa3D.zstart)/Aa3D.zdelta))
    elif num is None:
        raise RuntimeError('Provide either num or val')

    if((num < 0) or (num >= Aa3D.data.shape[0])):
        raise ValueError('Val or num outside of plot range')   
    
        

    if vmin is None:
        vmin = Aa3D.data[num].min()
    if vmax is None:
        vmax = Aa3D.data[num].max()
    
    pxaxis = np.linspace(Aa3D.xaxis[0]-0.5*Aa3D.xdelta, Aa3D.xaxis[-1]+0.5*Aa3D.xdelta, (Aa3D.xaxis.shape[0]+1))
    pyaxis = np.linspace(Aa3D.yaxis[0]-0.5*Aa3D.ydelta, Aa3D.yaxis[-1]+0.5*Aa3D.ydelta, (Aa3D.yaxis.shape[0]+1))
    
    
    if asp == 'true':
        asp = abs((pxaxis[-1]-pxaxis[0])/pyaxis[-1]-pyaxis[0])*(40/1.2)*(1376/1024)
    elif asp == 'ptrue':
        asp = float(Aa3D.xaxis.size)/(Aa3D.yaxis.size)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))

    meshy, meshx = np.meshgrid(pyaxis, pxaxis)
    ax.pcolormesh(meshy, meshx, Aa3D.data[:,:,num], cmap=cmap, vmin=vmin, vmax=vmax, **plotkwargs)
    ax.set_xlabel('Theta [deg]', fontsize = 12)
    ax.set_ylabel('Energy [eV]', fontsize = 12)    

    if plotrange is None:
        ax.axis([pyaxis[0],pyaxis[-1], pxaxis[0],pxaxis[-1]])
    else:
        ax.axis(plotrange)
    
    return ax

        
        
def plotkpath3D(Aa3D, path, points=100, pathtype='kcor', Estart=None, Eend=None,
                cmap='Spectral_r', vmin=None, vmax=None,
              ax=None, size=3.0, asp=3/4, plotkwargs={}):
    """plot a slice along a kpath vs energy

    Args:
        path (list): defines the kpath, can be [kx0,ky0,kx1,ky1] with the path
                     going from k0 to k1, or a list containing arrays [kx,ky] with the full path
        points (int): number of points along the path, only necessary if path is not given as arrays
        pathtype (str): currently only 'kcor' supported
        cmap (str): matplotlib colormap
        vmin (float): start of the colormap
        vmax (float): end of the colormap
        ax (axis): axis object to plot in
        size (float): size of the generated figure
        asp (float): aspect ratio
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """
    
    if Estart is None:
        Estart = Aa3D.xaxis[0]
    if Eend is None:
        Eend = Aa3D.xaxis[-1]
    
    if len(path) == 4:
        kx = np.linspace(path[0],path[2],points)
        ky = np.linspace(path[1],path[3],points)
    elif len(path) == 2:
        kx = path[0]
        ky = path[1]
    else:
        raise ValueError('Invalid path')
        
 
    Estartn,Eendn = int((Estart-Aa3D.xstart)/Aa3D.xdelta),int((Eend-Aa3D.xstart)/Aa3D.xdelta)     
    Eb = np.array(Aa3D.xaxis[Estartn:Eendn])            

    if pathtype == 'k':
        raise NotImplementedError('k path is not implemented')
    elif pathtype == 'kcor':
        #make meshgrid:
        Eb_m,kx_m = np.meshgrid(Eb,kx)
        Eb_m,ky_m = np.meshgrid(Eb,ky)
        #get data:
        data = Aa3D.atk(Eb_m,kx_m,ky_m, makegrid = False).T
    else:
        raise ValueError('pathtype unknown')
        
    #make k_abs
    if (kx[0] < 0) or (ky[0] < 0):
        k0 = -np.sqrt(kx[0]**2 + ky[0]**2)
    else:
        k0 = np.sqrt(kx[0]**2 + ky[0]**2)
    kabs = np.sqrt((kx-kx[0])**2 + (ky-ky[0])**2) + k0

        
    if vmin is None:
        vmin = data.min()
    if vmax is None:
        vmax = data.max()


    #plotting: (change for p-colormesh?)
    ystart = kabs[0]
    yend = kabs[-1]
    xstart = Eb[0]
    xend = Eb[-1]
    ratio = (yend-ystart)/(xend-xstart)*asp

    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))

    ax.imshow(np.flipud(data), cmap=cmap, vmin=vmin, vmax=vmax, interpolation='none', extent=[ystart,yend,xstart,xend], aspect=ratio, **plotkwargs)      
    ax.set_xlabel('kx [A-1]', fontsize = 12)
    ax.set_ylabel('Eb [eV]', fontsize = 12)   
    return ax
    
def plotkpath2_3D(Aa3D, path, points=200, pathtype='kcor', Estart=None, Eend=None, color='grey',
                line=False, klabels=None, cmap='Spectral_r', vmin=None, vmax=None,
              ax=None, size=3.0, asp=3/4, plotkwargs={}):
    """plot a slice along a multipoint kpath vs energy

    Args:
        path (list): defines the kpath, is a list of [kx,ky] waypoints, minimum length 2
                    as [[kx0,ky0],[kx1,ky1],[kx2,ky2],...]
        points (int): number of points along every path segment
        Estart (float): energy range lower limit
        Eend (float): energy range upper limit
        color (str): color of k-waypoint indicator lines
        line (bool): mark k-waypoints with lines if true
        klabels (list): list of labels to mark the waypoint lines with
        cmap (str): matplotlib colormap
        vmin (float): start of the colormap
        vmax (float): end of the colormap
        ax (axis): axis object to plot in
        size (float): size of the generated figure
        asp (float): aspect ratio
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """
    
    if Estart is None:
        Estart = Aa3D.xaxis[0]
    if Eend is None:
        Eend = Aa3D.xaxis[-1]

    kx = np.zeros(points*(len(path)-1)+1)
    ky = np.array(kx)
    kabs = np.array(kx) #for plotting later
    kticks = np.zeros(len(path))
    if klabels is None:
        klabels = ['']*len(path) ##number of points going from point to point        

    for i in range(len(path)-1):
        kxsub = np.linspace(path[i][0],path[i+1][0],points+1)
        kysub = np.linspace(path[i][1],path[i+1][1],points+1)
        kx[i*points:(i+1)*points+1] = kxsub
        ky[i*points:(i+1)*points+1] = kysub
        kabs[i*points:(i+1)*points+1] = np.sqrt((kxsub-kxsub[0])**2+(kysub-kysub[0])**2) + kabs[i*points]
             
    kticks = kabs[points*np.arange(len(path))]
 
    Estartn,Eendn = int((Estart-Aa3D.xstart)/Aa3D.xdelta),int((Eend-Aa3D.xstart)/Aa3D.xdelta)     
    Eb = np.array(Aa3D.xaxis[Estartn:Eendn])            

    if pathtype == 'k':
        raise NotImplementedError('k path is not implemented')
    elif pathtype == 'kcor':
        #make meshgrid:
        Eb_m,kx_m = np.meshgrid(Eb,kx)
        Eb_m,ky_m = np.meshgrid(Eb,ky)
        #get data:
        data = Aa3D.atk(Eb_m,kx_m,ky_m, makegrid = False).T
    else:
        raise ValueError('pathtype unknown')
        
    Eb_m,kabs_m = np.meshgrid(Eb,kabs)

    if vmin is None:
        vmin = data.min()
    if vmax is None:
        vmax = data.max()

    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))

    ax.pcolormesh(kabs_m,Eb_m,data.T, cmap=cmap, vmin=vmin, vmax=vmax, **plotkwargs)     
    ax.set_xlabel('k', fontsize = 12)
    ax.set_ylabel('Eb [eV]', fontsize = 12)   
    ax.set_xticks(kticks)
    ax.set_xticklabels(klabels, size = 15)
    if line:
        if len(kticks) > 2:
            for i in range(len(kticks)-2):
                ax.plot(np.array([kticks[i+1]]*2),np.array([-10e10,10e10]), ls = 'dashed', c = color)
    ax.axis([kabs[0],kabs[-1],Estart,Eend])
    return fig,ax        
    
def plotEDC3D(Aa3D, path, pathtype='k', NEDC=20, avg=1, avgtype='l', offset=0.3, 
              Estart=None, Eend=None, color='black',
              ax=None, plotrange=None, size=3.0, asp=3/4, plotkwargs={}):
    """plot a range of EDCs along a path

    Args:
        path (list): defines the kpath, can be [kx0,ky0,kx1,ky1] with the path
                     going from k0 to k1, or a list containing arrays [kx,ky] with the full path
        pathtype (str): type of the path (see AaData3D.getEDCs for info)
        NEDC (int): number of EDCs
        avg (int): number of EDCs to average
        avgtype (str): average type l (line), s (square) or c (circle)
        offset (float): offset between curves in the EDC stack
        Estart (float): energy range lower limit
        Eend (float): energy range upper limit
        color (str): color of k-waypoint indicator lines
        ax (axis): axis object to plot in
        plotrange (list): set the plotrange        
        size (float): size of the generated figure
        asp (float): aspect ratio
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """
    
    
    if Estart is None:
        Estart = Aa3D.xaxis[0]
    if Eend is None:
        Eend = Aa3D.xaxis[-1]
    
    
    if not ((len(path) == 4) or (len(path) == 2)):
        raise ValueError('Invalid path')
    
    if not ((avgtype == 'l') or (pathtype == 'sq') or (pathtype == 'c')):
        raise ValueError('Invalid avg type')
    
    
    EDCs,Escale,thpath,phipath = Aa3D.getEDCs(path,pathtype,NEDC,avg,avgtype,Estart,Eend)
    offsetavg = 0.0
    for i in range(NEDC):
        offsetavg += np.sum(EDCs[i][:10])
    doffset = offset*offsetavg/(10*NEDC)   
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))
    

    Estartn,Eendn = int((Estart-Aa3D.xstart)/Aa3D.xdelta),int((Eend-Aa3D.xstart)/Aa3D.xdelta)

    for i in range(NEDC):
        ax.plot(Aa3D.xaxis[Estartn:Eendn],(EDCs[i]+i*doffset), **plotkwargs)
    plotmin,plotmax = ax.axis()[2],ax.axis()[3]

    ax.set_xlabel('Energy [eV]', fontsize = 12)
    ax.set_ylabel('Photocurrent [arb. un.]', fontsize = 12)            
    ax.axis([Estart,Eend,plotmin,plotmax])
    
    return fig,ax
    
def get_kmeshCE(Aa3D, atE, extend = True):
    """Get kmesh at energy from angles in Aa3D
    Args:
        Aa3D: Aa3D object
        atE: energy to get the kmesh at
        extend: if True, return the dimensions+1 as kmesh, for use with pcolormesh
        
    Returns:
        kmeshx, kmeshy, a tuple of arrays with the plot mesh
        
        """
    
    if extend:
        pyaxis = np.linspace(Aa3D.yaxis[0]-0.5*Aa3D.ydelta, Aa3D.yaxis[-1]+0.5*Aa3D.ydelta, (Aa3D.yaxis.shape[0]+1))
        pzaxis = np.linspace(Aa3D.zaxis[0]-0.5*Aa3D.zdelta, Aa3D.zaxis[-1]+0.5*Aa3D.zdelta, (Aa3D.zaxis.shape[0]+1))
    else:
        pyaxis = Aa3D.yaxis
        pzaxis = Aa3D.zaxis
    ameshth, ameshphi = np.meshgrid(pyaxis, pzaxis)

    kmeshx, kmeshy = kfromangles(ameshth,ameshphi,(Aa3D.p['EkatEF']+atE),
                    **{k:Aa3D.p[k] for k in ['theta_m','phi_m','theta0','phi0','alpha0']})
    return kmeshx, kmeshy

    

 
def plotCEk3D(Aa3D, num=None, val=None, cmap='Spectral_r', vmin=None, vmax=None, 
              Eint_n=1, Eint=None, pltbzn=True, rasterized=True, slitline=None,square=None,
             ax=None, plotrange=None, size=3.0, asp='true', plotkwargs={}):
    """k-corrected constant energy plot function


    Args:
        num (int): the xaxis index at which energy the slice is taken
        val (float): energy at which the slice is taken, overrides num
        cmap (str): matplotlib colormap
        vmin (float): start of the colormap
        vmax (float): end of the colormap
        Eint_n (int): number of slices to average over
        Eint (float): energy range to average over, overrides Eint_n
        pltbzn (bool): plot a square first Brillouin Zone
        rasterized (bool): plot the mesh rasterized, recommended for saving pdfs
        slitline (tuple): (th_x, th_y, width) plot a curve where the slit measures given angles
        square (str): if 'crop', crop the image so axis are proportional, if not None: extend the axes until proportional
        ax (axis): axis object to plot in
        plotrange (list): set the plotrange
        size (float): size of the generated figure
        asp (float): aspect ratio, when passed as 'true', calculate based on kx/ky range
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """
    
    if(not (Aa3D.p['type'] == 'corrected')):
        raise RuntimeError('k plots only work on corrected data')
                  
    if val is not None:
        num = int(round((val- Aa3D.xstart)/Aa3D.xdelta))
    elif num is None:
        raise RuntimeError('Provide either num or val')

    if((num < 0) or (num >= Aa3D.data.shape[0])):
        raise ValueError('Val or num outside of plot range')   

    if vmin is None:
        vmin = Aa3D.data[num].min()
    if vmax is None:
        vmax = Aa3D.data[num].max()

    if Eint is not None:
        Eint_n = int(Eint/Aa3D.xdelta)

    print(f'num = {num}, {Aa3D.xaxis[num]}')
    
    kmeshx, kmeshy = get_kmeshCE(Aa3D, Aa3D.xaxis[num], extend = True)

    kxstart,kxend,kystart,kyend = np.min(kmeshx),np.max(kmeshx),np.min(kmeshy),np.max(kmeshy)


    if asp == 'true':
        asp = abs((kxend-kxstart)/(kystart-kyend))
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))
    
    ax.pcolormesh(kmeshx, kmeshy, np.mean(Aa3D.data[num-Eint_n:num+Eint_n], axis = 0).T, 
                  cmap=cmap, vmin=vmin, vmax=vmax,rasterized=rasterized, **plotkwargs)
    ax.set_xlabel('kx [A-1]', fontsize = 12)
    ax.set_ylabel('ky [A-1]', fontsize = 12)    
    if pltbzn: #plot Brillouin zone
        if(Aa3D.p['bz'] == 'square'):
            bsize = np.pi/Aa3D.p['alatt']
            bzkx = np.array([bsize, -bsize, -bsize, bsize, bsize])
            bzky = np.array([bsize, bsize, -bsize, -bsize, bsize])
            gridkx = np.array([-5*bsize, -5*bsize,-4*bsize, -4*bsize,-3*bsize, -3*bsize,-2*bsize, -2*bsize,-1*bsize, -1*bsize,0, 0,1*bsize, 1*bsize,2*bsize, 2*bsize,3*bsize, 3*bsize,4*bsize, 4*bsize,5*bsize, 5*bsize]) #make 10x10
            gridky = np.array([-50*bsize, 50*bsize, 50*bsize, -50*bsize,-50*bsize, 50*bsize, 50*bsize, -50*bsize,-50*bsize, 50*bsize, 50*bsize, -50*bsize,-50*bsize, 50*bsize, 50*bsize, -50*bsize,-50*bsize, 50*bsize, 50*bsize, -50*bsize,-50*bsize, 50*bsize])
            ax.plot(gridkx, gridky, color = 'black')
            ax.plot(gridky, gridkx, color = 'black')
            ax.plot(bzkx, bzky, color = 'red')
    if slitline is not None: #slitline is (th,phi,dth)
        th,phi,dth = slitline
        thstart, thend = th-dth/2, th+dth/2
        slitth = np.linspace(thstart,thend,50)
        slitphi = np.linspace(phi,phi,50)
        slitkx, slitky = kfromangles(slitth,slitphi,(Aa3D.p['EkatEF']+Aa3D.xstart+Aa3D.xdelta*num), theta_m = Aa3D.p['theta_m'], phi_m = Aa3D.p['phi_m'], theta0 = Aa3D.p['theta0'], phi0 = Aa3D.p['phi0'], alpha0 = Aa3D.p['alpha0'])
        ax.plot(slitkx,slitky,color = 'red')
    if plotrange is None:
        if pltbzn:
            plotrange = [min(-bsize*1.1,kxstart),max(bsize*1.1,kxend),min(kystart,-bsize*1.1),max(kyend,bsize*1.1)]
        else:
            plotrange = [kxstart,kxend,kystart,kyend]
    if square is not None:
        ax_win = ax.get_window_extent()
        ax_asp = ax_win.width/ax_win.height
        plot_width,plot_height = plotrange[1] - plotrange[0],plotrange[3] - plotrange[2]
        x_center,y_center = (plotrange[1] + plotrange[0])/2,(plotrange[3] + plotrange[2])/2
        plot_asp = plot_width/plot_height
        if square == 'crop':
            if plot_asp <= ax_asp: #if plot is thinner than axis, cut height:
                plot_height = plot_width/ax_asp
            else:
                plot_width = ax_asp * plot_height #else cut width
        else:#if not cropping but extending: 
            if plot_asp <= ax_asp: #if plot is thinner than axis, extend width:
                plot_width = ax_asp * plot_height 
            else: #else extend length
                plot_height = plot_width/ax_asp
        plotrange = [x_center - plot_width/2, x_center + plot_width/2, y_center - plot_height/2, y_center + plot_height/2]
                

    ax.axis(plotrange)

    
    return ax

def plotkz3D(Aa3D, num=None, val=None, cmap='Spectral_r', vmin=None, vmax=None, 
              Eint_n=1, Eint=None, pltbzn=True, rasterized=True,
             ax=None, plotrange=None, size=3.0, asp='true', plotkwargs={}):
    """k-corrected constant energy plot function with zaxis set up as hv (kz dependence)


    Args:
        num (int): the xaxis index at which energy the slice is taken
        val (float): energy at which the slice is taken, overrides num
        cmap (str): matplotlib colormap
        vmin (float): start of the colormap
        vmax (float): end of the colormap
        Eint_n (int): number of slices to average over
        Eint (float): energy range to average over, overrides Eint_n
        pltbzn (bool): plot a square first Brillouin Zone
        rasterized (bool): plot the mesh rasterized, recommended for saving pdfs
        ax (axis): axis object to plot in
        plotrange (list): set the plotrange
        size (float): size of the generated figure
        asp (float): aspect ratio, when passed as 'true', calculate based on kx/kz range
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """
    
    
    if(not (Aa3D.p['type'] == 'corrected')):
        raise RuntimeError('k plots only work on corrected data')
                  
    if val is not None:
        num = int(round((val- Aa3D.xstart)/Aa3D.xdelta))
    elif num is None:
        raise RuntimeError('Provide either num or val')

    if((num < 0) or (num >= Aa3D.data.shape[0])):
        raise ValueError('Val or num outside of plot range')   

    if vmin is None:
        vmin = Aa3D.data[num].min()
    if vmax is None:
        vmax = Aa3D.data[num].max()

    if Eint is not None:
        Eint_n = int(Eint/Aa3D.xdelta)
    

    if not 'V0' in Aa3D.p:
        Aa3D.p['V0'] = 10.
        print('Using V0 = 10.0 because no V0 was set')
    if not 'WF' in Aa3D.p:
        Aa3D.p['WF'] = 4.3
        print('Using WF = 4.3 because no WF was set')
    
    pyaxis = np.linspace(Aa3D.yaxis[0]-0.5*Aa3D.ydelta, Aa3D.yaxis[-1]+0.5*Aa3D.ydelta, (Aa3D.yaxis.shape[0]+1))
    pzaxis = np.linspace(Aa3D.zaxis[0]-0.5*Aa3D.zdelta, Aa3D.zaxis[-1]+0.5*Aa3D.zdelta, (Aa3D.zaxis.shape[0]+1))
    ameshth, ameshhv = np.meshgrid(pyaxis, pzaxis)
    kmeshx, kmeshy, kmeshz = kzfromhv(ameshhv,ameshth,Aa3D.p['phi'],E = Aa3D.xaxis[num], Aaobj = Aa3D)
    kxstart,kxend,kzstart,kzend = np.min(kmeshx),np.max(kmeshx),np.min(kmeshz),np.max(kmeshz)

    if asp == 'true':
        asp = abs((kxend-kxstart)/(kzstart-kzend))
        
    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))
    
    ax.pcolormesh(kmeshx, kmeshz, np.mean(Aa3D.data[num-Eint_n:num+Eint_n], axis = 0).T, 
                  cmap=cmap, vmin=vmin, vmax=vmax,rasterized=rasterized, **plotkwargs)
    ax.set_xlabel('kx [A-1]', fontsize = 12)
    ax.set_ylabel('kz [A-1]', fontsize = 12)    
    if pltbzn: #plot Brillouin zone
        if(Aa3D.p['bz'] == 's' or Aa3D.p['bz'] == 'square'):
            bsize_z = np.pi/Aa3D.p['clatt']
            bsize = np.pi/Aa3D.p['alatt']
            ax.vlines([bsize*i for i in range(-20,20)], -2*bsize_z, 30*bsize_z)
            ax.hlines([bsize_z*i for i in range(0,30)], -20*bsize, 20*bsize)

    if plotrange is None:
        if pltbzn:
            ax.axis([min(-bsize*1.1,kxstart),max(bsize*1.1,kxend),kzstart,max(kzend,bsize_z*1.1)])
        else:
            ax.axis([kxstart,kxend,kzstart,kzend])
    else:
        ax.axis(plotrange)
    
    return ax            

        
        
        
        
################## Hybrid class plotters #######################################

################## Spin Pol ##################################################
        
def plotcpspol(CPSAobj, integrated=False, zeroline=True, color='green', ax=None, 
               plotrange=None, size=3.0, asp=3/4., plotkwargs={}):
    """Plot the circularly polarized spin arpes


    Args:
        integrated (bool): plot the integrated (sum) of the signals in the background
        zeroline (bool): plot a line marking the Fermi energy
        color (str): color of the markers
        ax (axis): axis object to plot in
        plotrange (list): set the plotrange
        size (float): size of the generated figure
        asp (float): aspect ratio, when passed as 'true', calculate based on kx/kz range
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """

    
    CP_up = CPSAobj.data['CP_up'].data
    CP_dn = CPSAobj.data['CP_dn'].data
    CM_up = CPSAobj.data['CM_up'].data
    CM_dn = CPSAobj.data['CM_dn'].data

    xaxis = CPSAobj.data['CP_up'].xaxis
    
    N = np.sqrt(CP_dn*CM_up)-np.sqrt(CP_up*CM_dn)
    D = np.sqrt(CP_up*CM_dn)+np.sqrt(CP_dn*CM_up)    
    spin_int = (CM_up + CM_dn + CP_up + CP_dn)    

    SPA = N/D

    if plotrange is None:
        plotrange = [xaxis[0],xaxis[-1], -0.1, 0.1]


    dSPA = SPA*np.sqrt(0.25*(CP_up+CP_dn+CM_up+CM_dn)*(N**-2+D**-2))
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))

    ax.errorbar(xaxis, SPA, yerr = dSPA, marker = 'o', color = color, capsize = 4, **plotkwargs)
    ax.axis(plotrange)
    
    if zeroline:
        ax.hlines(0.0, xaxis[0]-1., xaxis[-1]+1.)
    
    if integrated == True:
        ax_t = ax.twinx()
#        ax_t.fill_between(xaxis, 0, np.sum(CPSAobj.data, axis = 1), color = 'grey')
        ax_t.fill_between(xaxis, 0, spin_int/np.max(spin_int), color = 'grey', alpha = 0.3)
        ax_t.axis([xaxis[0],xaxis[-1],0,1.2])
        ax_t.set_ylabel('Spin integrated', fontsize = 12)   
        ax.set_zorder(ax_t.get_zorder() + 1)
        ax.patch.set_visible(False)

    
    ax.set_xlabel('Energy [eV]', fontsize = 12)
    ax.set_ylabel('Polarization', fontsize = 12)   

    
    return ax


def plotcoan(CPSAobj, ax=None, plotrange=None, size=3.0, asp=3/4., plotkwargs={}):
    """Plot the colinear and anti parallel components


    Args:
        ax (axis): axis object to plot in
        plotrange (list): set the plotrange
        size (float): size of the generated figure
        asp (float): aspect ratio, when passed as 'true', calculate based on kx/kz range
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """    
    
    
    CP_up = CPSAobj.data['CP_up'].data
    CP_dn = CPSAobj.data['CP_dn'].data
    CM_up = CPSAobj.data['CM_up'].data
    CM_dn = CPSAobj.data['CM_dn'].data

    xaxis = CPSAobj.data['CP_up'].xaxis

    co = np.sqrt(CM_up*CP_dn)
    an = np.sqrt(CM_dn*CP_up)
    
    if plotrange is None:
        plotrange = [xaxis[0],xaxis[-1], min(np.min(co), np.min(an)), max(np.max(co), np.max(an))] ##std plotrange
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))
    

    ax.plot(xaxis, co, marker = 'o', color = 'red', label = 'co', **plotkwargs)
    ax.plot(xaxis, an, marker = 'o', color = 'blue', label = 'an', **plotkwargs)
    ax.legend(loc = 'upper right', shadow = False, fontsize = 12)
    
    ax.axis(plotrange)
    
    ax.set_xlabel('Energy [eV]', fontsize = 12)
    ax.set_ylabel('Intensity [arb. un.]', fontsize = 12)       
    
    return ax
    
    

def plotspinpol(SpinPolobj, normtype='sens', integrated=False, zeroline=True, ax=None,
                Estart=None, Eend=None, Estart_n=None, Eend_n=None, color='green',
               plotrange=None, size=3.0, asp=3/4., plotkwargs={}):
    """Plot linear spin polarization (spinUp-spinDn)


    Args:
        normtype (str): type of normalization 'sens' takes sensitivity, 'int' integrates the total signal    
        integrated (bool): plot the integrated (sum) of the signals in the background
        zeroline (bool): plot a line marking the Fermi energy
        ax (axis): axis object to plot in
        Estart (float): integration energy range lower limit
        Eend (float): integration energy range upper limit        
        Estart_n (int): integration energy range lower limit bin number
        Eend_n (int): integration energy range upper limit bin number
        color (str): color of the markers        
        plotrange (list): set the plotrange
        size (float): size of the generated figure
        asp (float): aspect ratio, when passed as 'true', calculate based on kx/kz range
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """    
    

    up = SpinPolobj.data['up'].data
    dn = SpinPolobj.data['dn'].data
    
    xaxis = SpinPolobj.data['up'].xaxis

 
    if normtype == 'sens':
        up_fac = SpinPolobj.p['sens_up']
        dn_fac = SpinPolobj.p['sens_dn']
    elif normtype == 'int':
        if Estart is None: 
            Estart = xaxis[0]
        if Eend is None: 
            Eend = xaxis[10]
        nstart = max(0,np.argmin(np.abs(xaxis-Estart)))
        nend = min(xaxis.size,np.argmin(np.abs(xaxis-Eend)))
        up_fac = np.sum(up[nstart:nend])/(nend-nstart)
        dn_fac = np.sum(dn[nstart:nend])/(nend-nstart)
    elif normtype == 'int_n':
        if Estart_n is None: 
            Estart_n = 0
        if Eend_n is None: 
            Eend_n = 10
        up_fac = np.sum(up[nstart:nend])/(nend-nstart)
        dn_fac = np.sum(dn[nstart:nend])/(nend-nstart)        

    up = up_fac*SpinPolobj.data['up'].data
    dn = dn_fac*SpinPolobj.data['dn'].data
        
    spin_int = up + dn
    
    
    dup = np.sqrt(up)
    ddn = np.sqrt(dn)
    
    dpol = np.sqrt((2*dn*dup)**2 + (2*up*ddn)**2)/spin_int**2
    
    pol = (up-dn)/spin_int
    
    if plotrange is None:
        plotrange = [xaxis[0],xaxis[-1], -0.1, 0.1] ##std plotrange
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))
    

    ax.errorbar(xaxis, pol, yerr = dpol, marker = 'o', color = 'green', capsize = 4, **plotkwargs)
    ax.axis(plotrange)        
        
    if zeroline:
        ax.hlines(0.0, xaxis[0]-1., xaxis[-1]+1.)
    
    if integrated == True:
        ax_t = ax.twinx()
        ax_t.fill_between(xaxis, 0, spin_int/np.max(spin_int), color = 'grey', alpha = 0.3)
        ax_t.axis([xaxis[0],xaxis[-1],0,1.2])
        ax_t.set_ylabel('Spin integrated', fontsize = 12)   
        ax.set_zorder(ax_t.get_zorder() + 1)
        ax.patch.set_visible(False)        
    
        

    ax.set_xlabel('Energy [eV]', fontsize = 12)
    ax.set_ylabel('Polarization [arb. un.]', fontsize = 12)       
    
    return ax    


def plotspin(SpinPolobj, normtype='sens', integrated=False, zeroline=True, ax=None,
                Estart=None, Eend=None, Estart_n=None, Eend_n=None, color='green',
               size=3.0, asp=3/4., plotkwargs={}):
    """Plot linear spin polarization components


    Args:
        normtype (str): type of normalization 'sens' takes sensitivity, 'int' integrates the total signal    
        integrated (bool): plot the integrated (sum) of the signals in the background
        zeroline (bool): plot a line marking the Fermi energy
        ax (axis): axis object to plot in
        Estart (float): integration energy range lower limit
        Eend (float): integration energy range upper limit        
        Estart_n (int): integration energy range lower limit bin number
        Eend_n (int): integration energy range upper limit bin number
        color (str): color of the markers        
        plotrange (list): set the plotrange
        size (float): size of the generated figure
        asp (float): aspect ratio, when passed as 'true', calculate based on kx/kz range
        plotkwargs (dict): additional arguments to be passed to the plot function
    Returns:
        ax (axis): axis object that was plotted in

    """    
    
    up = SpinPolobj.data['up'].data
    dn = SpinPolobj.data['dn'].data

    xaxis = SpinPolobj.data['up'].xaxis

    if normtype == 'sens':
        up_fac = SpinPolobj.p['sens_up']
        dn_fac = SpinPolobj.p['sens_dn']
    elif normtype == 'int':
        if Estart is None: 
            Estart = xaxis[0]
        if Eend is None: 
            Eend = xaxis[10]
        nstart = max(0,np.argmin(np.abs(xaxis-Estart)))
        nend = min(xaxis.size,np.argmin(np.abs(xaxis-Eend)))
        up_fac = np.sum(up[nstart:nend])/(nend-nstart)
        dn_fac = np.sum(dn[nstart:nend])/(nend-nstart)
    elif normtype == 'int_n':
        if Estart_n is None: 
            Estart_n = 0
        if Eend_n is None: 
            Eend_n = 10
        up_fac = np.sum(up[nstart:nend])/(nend-nstart)
        dn_fac = np.sum(dn[nstart:nend])/(nend-nstart)        
 

    
    if ax is None:
        fig, ax = plt.subplots(figsize=(size,size*asp))

    
    

    ax.plot(xaxis, up_fac*up, marker = 'o', color = 'blue', label = 'up')
    ax.plot(xaxis, dn_fac*dn, marker = 'o', color = 'red', label = 'dn')
    ax.legend(loc = 'upper right', shadow = False, fontsize = 12)
    
    ax.set_xlabel('Energy [eV]', fontsize = 12)
    ax.set_ylabel('Intensity [arb. un.]', fontsize = 12)       
    
    return ax
