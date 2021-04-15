# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:48:41 2016

@author: berend
"""

import numpy as np
from scipy.interpolate import RectBivariateSpline, interp2d, RegularGridInterpolator, interpn




def kfromangles(theta,phi,Ek,theta_m = 0, phi_m = 0, alpha0 = 0, theta0 = 0, phi0 = 0, Aaobj = None):
    """return (kx, ky) from angles (theta, phi), alpha0,theta0,phi0 are offset angles for the sample,
    theta_m, phi_m are offset angles for the manipulator, all angles in degrees"""

    # sin = np.sin
    # cos = np.cos

    if Aaobj is not None:
        theta_m = Aaobj.p['theta_m']
        phi_m = Aaobj.p['phi_m']
        theta0 = Aaobj.p['theta0']
        phi0 = Aaobj.p['phi0']
        alpha0 = Aaobj.p['alpha0']

    
    theta = theta + theta_m ## add offset angles for manipulator to theta and phi
    phi = phi + phi_m ##26 for manipulator 0
    

    
    #translate to radians:
    phi /= 180/np.pi  
    theta /= 180/np.pi  
    phi0 /= 180/np.pi  
    theta0 /= 180/np.pi  
    alpha0 /= 180/np.pi

    
    #calculate k_m (hypothetical k for if the sample was perfectly aligned) (k_m is a unit vector)
    kxm = np.sin(theta)
    kym = np.sin(phi)*np.cos(theta)
    kzm = np.sqrt(1-kxm*kxm-kym*kym)


    #now calculate the actual kx, ky, kz
    # Ek = Ek.astype('float64')  # this is to prevent a bug in np that raises an error taking the
    # sqrt of 0
    kx = 0.5123*np.sqrt(Ek)*(kxm*np.cos(alpha0)*np.cos(theta0) + kym*np.sin(alpha0)*np.cos(theta0) - kzm*np.sin(theta0))
    ky = 0.5123*np.sqrt(Ek)*(kxm*(-np.sin(alpha0)*np.cos(phi0) + np.sin(phi0)*np.sin(theta0)*np.cos(alpha0)) +
                             kym*(np.sin(alpha0)*np.sin(phi0)*np.sin(theta0) + np.cos(alpha0)*np.cos(phi0)) +
                             kzm*np.sin(phi0)*np.cos(theta0))
    
    #return values:
    return (kx, ky)
   

def anglesfromk(kx,ky,Ek,theta_m = 0, phi_m = 0, 
                alpha0 = 0, theta0 = 0, phi0 = 0, Aaobj = None):
    """return (theta, phi) from k (kx, ky), alpha0,theta0,phi0 are offset angles for the sample,
    theta_m, phi_m are offset angles for the manipulator, all angles in degrees"""
    sin = np.sin
    cos = np.cos
    
    if Aaobj is not None:
        theta_m = Aaobj.p['theta_m']
        phi_m = Aaobj.p['phi_m']
        theta0 = Aaobj.p['theta0']
        phi0 = Aaobj.p['phi0']
        alpha0 = Aaobj.p['alpha0']
    
    #translate to radians:
    phi0 /= 180/np.pi  
    theta0 /= 180/np.pi  
    alpha0 /= 180/np.pi  
    
    #turn k into a unit vector
    kx = kx/(0.5123*np.sqrt(Ek))
    ky = ky/(0.5123*np.sqrt(Ek))
    
    #calculate kz
    kz = np.sqrt(1-kx*kx-ky*ky)

    #calculate km (the hypothetical k for if the sample was perfectly aligned)
    kxm = kx*cos(alpha0)*cos(theta0) + ky*(-sin(alpha0)*cos(phi0) + sin(phi0)*sin(theta0)*cos(alpha0)) + kz*(sin(alpha0)*sin(phi0) + sin(theta0)*cos(alpha0)*cos(phi0))
    kym = kx*sin(alpha0)*cos(theta0) + ky*(sin(alpha0)*sin(phi0)*sin(theta0) + cos(alpha0)*cos(phi0)) + kz*(sin(alpha0)*sin(theta0)*cos(phi0) - sin(phi0)*cos(alpha0))        
    
    #now calculate angles from kxm, kym:
    theta = np.arcsin(kxm)
    phi = np.arcsin(kym/np.cos(theta))

    
    return (theta*180/np.pi - theta_m),(phi*180/np.pi - phi_m)
    
def kzfromhv(hv,theta,phi, E = 0.0, V0 = None, WF = None, theta_m = 0, phi_m = 0, 
             alpha0 = 0, theta0 = 0, phi0 = 0, Aaobj = None):
    """return k (kx, ky, kz) from (theta, hv). E is defined as negative energies are below EF (filled states)
    alpha0,theta0,phi0 are offset angles for the sample,
    theta_m, phi_m are offset angles for the manipulator, all angles in degrees"""
    sin = np.sin
    cos = np.cos
    
    if Aaobj is not None:
        theta_m = Aaobj.p['theta_m']
        phi_m = Aaobj.p['phi_m']
        theta0 = Aaobj.p['theta0']
        phi0 = Aaobj.p['phi0']
        alpha0 = Aaobj.p['alpha0']
        WF = Aaobj.p['WF']
        V0 = Aaobj.p['V0']
    else:
        if (V0 is None) or (WF is None):
            raise RuntimeError('Need V0 and WF to run, or valid Aaobj')

    
    theta = theta + theta_m ## add offset angles for manipulator to theta and phi
    phi = phi + phi_m ##26 for manipulator 0

    
    #translate to radians:
    phi /= 180/np.pi  
    theta /= 180/np.pi  
    phi0 /= 180/np.pi  
    theta0 /= 180/np.pi  
    alpha0 /= 180/np.pi      

    
    #calculate k_m (hypothetical k for if the sample was perfectly aligned) (k_m is a unit vector)
    kxm = np.sin(theta)
    kym = np.sin(phi)*np.cos(theta)
    kzm = np.sqrt(1-kxm*kxm-kym*kym)
    
    Ek = hv-WF+E
    
    #now calculate the actual kx, ky, kz
    kx = 0.5123*np.sqrt(Ek)*(kxm*cos(alpha0)*cos(theta0) + kym*sin(alpha0)*cos(theta0) - kzm*sin(theta0))
    ky = 0.5123*np.sqrt(Ek)*(kxm*(-sin(alpha0)*cos(phi0) + sin(phi0)*sin(theta0)*cos(alpha0)) + kym*(sin(alpha0)*sin(phi0)*sin(theta0) + cos(alpha0)*cos(phi0)) + kzm*sin(phi0)*cos(theta0))
    kz = np.sqrt(0.5123**2*(Ek+V0)-kx*kx-ky*ky)
    
    
    #return values:
    return (kx,ky,kz)
    
    #note: see details for this in lab book Nov 2018
def anglesfromkz(kx,kz,E, V0 = None, WF = None, theta0 = 0, Aaobj = None):
    """Return angles/hv from  (kx,kz,E) grid
    Assumes that phi0 and ky are zero, otherwise transformation is impossible
    Args:
        kx: np array of kx points
        kz: np array of kz points
        E: np array of energies (E - EF)
        WF: work function
        V0: inner potential""" 
    
    if Aaobj is not None:
        theta0 = Aaobj.p['theta0']
        WF = Aaobj.p['WF']
        V0 = Aaobj.p['V0']
    else:
        if (V0 is None) or (WF is None):
            raise RuntimeError('Need V0 and WF to run, or valid Aaobj')
    
    #translate to radians:
    theta0 /= 180/np.pi  
    
    theta = np.arctan2(kx, np.sqrt(kz**2 - 0.5123**2*V0))
    hv = (kx*kx + kz*kz)/0.5123**2 - V0 + WF - E

    #add offset to theta0:
    theta += theta0
    
    return (theta*180/np.pi, hv)

            


    
    
def findkrange(Eaxis, thetaaxis, phiaxis,theta_m = 0, phi_m = 0, alpha0 = 0, theta0 = 0, phi0 = 0):
    """Return the limits to the kx, ky range from Eaxis, thetaaxis, phiaxis and offset angles"""
    #make subsampling of each axis:
    try:
        Eax = np.linspace(Eaxis[0], Eaxis[-1], 10)
    except:
        Eax = np.linspace(Eaxis, Eaxis, 10)
    thax = np.linspace(thetaaxis[0], thetaaxis[-1], 10)
    phiax = np.linspace(phiaxis[0], phiaxis[-1], 10)
    inrange = np.meshgrid(thax,phiax,Eax,indexing='ij')
    inrange = (inrange[0],inrange[1],inrange[2]) #turn into tuple
    outrange = kfromangles(*inrange, theta_m = theta_m, phi_m = phi_m, alpha0 = alpha0, theta0 = theta0, phi0 = phi0)

    kxmin = outrange[0].min()
    kxmax = outrange[0].max()    
    kymin = outrange[1].min()
    kymax = outrange[1].max()  

    return kxmin,kxmax,kymin,kymax
    
def findsampling(thetaaxis, phiaxis, sampletype, size = None, alpha0 = 0):
    """Find a suitable number of points and return as (kxpoints, kypoints)"""
    thetapoints = thetaaxis.shape[0]
    phipoints = phiaxis.shape[0]
    if(sampletype == 'fixed'):#fixed type
        return (512,512)
    elif(sampletype == 'square'):
        points = max(thetapoints,phipoints)
        points = min(512,points) #cut off at 512
        return (points, points)
    elif(sampletype == 'prop'):#proportional
        kxpoints = abs(np.cos(alpha0)*thetapoints - np.sin(alpha0)*phipoints)
        kypoints = abs(np.sin(alpha0)*thetapoints + np.cos(alpha0)*phipoints)
        return (kxpoints, kypoints)
    elif(sampletype == 'custom'):
        if not (size == None):
            if (len(size) == 2):
                return (size[0], size[1])
            else:
                print('Provide size as a tuple as (kxn, kyn)')
        else:
            print('Please provide size')

def convertslicetok(data, Ek, thetaaxis, phiaxis, sampletype = 'prop', size = None, theta_m = 0, phi_m = 0, alpha0 = 0, theta0 = 0, phi0 = 0):
    """Convert a constant energy slice to k"""
    kxmin, kxmax, kymin, kymax = findkrange(Ek, thetaaxis, phiaxis,theta_m = theta_m, phi_m = phi_m, alpha0 = alpha0, theta0 = theta0, phi0 = phi0)
    kxpoints, kypoints = findsampling(thetaaxis, phiaxis, sampletype, size = size, alpha0 = alpha0)
    
    datak = np.zeros((kxpoints, kypoints))
    kxaxis = np.linspace(kxmin, kxmax, kxpoints)
    kyaxis = np.linspace(kymin, kymax, kypoints)

    
    


    #now interpolate    
    interp = interp2d(phiaxis, thetaaxis, data, bounds_error=False, fill_value = 0.)    
    
    for i in range(datak.shape[0]):
        for j in range(datak.shape[1]):
            angles = anglesfromk(kxaxis[i], kyaxis[j], Ek, theta_m = theta_m, phi_m = phi_m, alpha0 = alpha0, theta0 = theta0, phi0 = phi0)
            datak[i,j] = interp(*angles)
    
    datak = interp(kxaxis, kyaxis)
    

    return datak, kxaxis, kyaxis
    
def atk(kx, ky, Eb, thetaaxis, phiaxis, data, Ek, theta_m = 0, phi_m = 0, alpha0 = 0, theta0 = 0, phi0 = 0):
    """"""
    # flatten kx, ky (for when they are meshgrids), then have a vectorized
    # function available that converts these vectors to angle space, and do a linear interpolation
    # via RectBivariateSpline (grid turned off) and after that use a vectorized
    # numpy function to make points zero if they fall outside the input angle range



    if(not (type(kx) == np.ndarray)):
        kx = np.array([kx])
    if(not (type(ky) == np.ndarray)):
        ky = np.array([ky])
    

    th,phi = anglesfromk(kx,ky,Ek+Eb,theta_m = theta_m, phi_m = phi_m, alpha0 = alpha0, theta0 = theta0, phi0 = phi0)

    thmin = min(thetaaxis[0],thetaaxis[-1])
    thmax = max(thetaaxis[0],thetaaxis[-1])
    phimin = min(phiaxis[0],phiaxis[-1])
    phimax = max(phiaxis[0],phiaxis[-1])



    #interpolator:
    interp = RectBivariateSpline(thetaaxis, phiaxis, data)
    
        
    dataout = interp(kx, ky, grid = False)
    dataout[np.where((th < thmin) | (th > thmax) | (phi < phimin) | (phi > phimax))] = 0

    return dataout
    
def atk2(kx, ky, Eb, thetaaxis, phiaxis, Eaxis, data, Ek, theta_m = 0, phi_m = 0, alpha0 = 0, theta0 = 0, phi0 = 0):
    """"""
    # flatten kx, ky (for when they are meshgrids), then have a vectorized
    # function available that converts these vectors to angle space, and do a linear interpolation
    # via RectBivariateSpline (grid turned off) and after that use a vectorized
    # numpy function to make points zero if they fall outside the input angle range



    if(not (type(kx) == np.ndarray)):
        kx = np.array([kx])
    if(not (type(ky) == np.ndarray)):
        ky = np.array([ky])
    if(not (type(Eb) == np.ndarray)):
        Eb = np.array([Eb])
    

    th,phi = anglesfromk(kx,ky,Ek+Eb,theta_m = theta_m, phi_m = phi_m, alpha0 = alpha0, theta0 = theta0, phi0 = phi0)

    thmin = min(thetaaxis[0],thetaaxis[-1])
    thmax = max(thetaaxis[0],thetaaxis[-1])
    phimin = min(phiaxis[0],phiaxis[-1])
    phimax = max(phiaxis[0],phiaxis[-1])
    Emin = min(Eaxis[0], Eaxis[-1])
    Emax = max(Eaxis[0], Eaxis[-1])



    #interpolator:
    interp = RegularGridInterpolator((Eaxis,thetaaxis, phiaxis), data)
    
        
    dataout = interp(Eb, kx, ky)
    dataout[np.where((th < thmin) | (th > thmax) | (phi < phimin) | (phi > phimax) |(Eb < Emin) | (Eb > Emax))] = 0

    return dataout

    
def atk3(dataobj, Eb, kx, ky, makegrid = True):
    
    if makegrid:
        if(not (type(kx) == np.ndarray)):
            kx = np.array([kx])
        if(not (type(ky) == np.ndarray)):
            ky = np.array([ky])
        if(not (type(Eb) == np.ndarray)):
            Eb = np.array([Eb])
        gridEb,gridkx,gridky = np.meshgrid(Eb,kx,ky)
    else: #if providing own grid
        gridEb = Eb
        gridkx = kx
        gridky = ky
        
    gridth,gridphi = anglesfromk(gridkx,gridky,dataobj.EkatEF+gridEb,theta_m = dataobj.theta_m, phi_m = dataobj.phi_m, alpha0 = dataobj.alpha0, theta0 = dataobj.theta0, phi0 = dataobj.phi0)

    intgrid = np.array([gridEb,gridth,gridphi]).transpose((2,1,3,0))
    
    return interpn((dataobj.xaxis,dataobj.yaxis,dataobj.zaxis),dataobj.data, intgrid, bounds_error = False, fill_value = 0.)
    

