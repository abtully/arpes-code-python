#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 19:27:58 2017

@author: berend
"""

import numpy as np
from lmfit import Parameters, minimize
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve, argrelextrema
import scipy.integrate as inte

#### Curve fitting

#fitfuncs

##lorentzian
def Lorentz(x, A, x0, g):
    """Lorentzian function, g is the FWHM"""
    g /= 2
    return A*(g*g)/((x-x0)*(x-x0)+g*g)
    
    

def Gauss(x, A, x0, sig):
    """Gaussian function, sig is the FWHM"""
    Sig = sig/(2*np.sqrt(2*np.log(2)))
    return A * np.exp(-(x-x0)**2 /(2*Sig*Sig))

    
def Voigt(x, A, x0, g, sig):
    """Voigt profile, sig and g are gaussian and lorentzian fwhm respectively"""
    dx = x[1]-x[0]
    sig_n = int(5*sig/dx)
    xgauss = np.linspace(-sig_n*dx,sig_n*dx,2*sig_n + 1) #minimal 1 pt, same dx scale. Lenght is N, need to extend x by (N-1)/2 on each side
    x_n = np.linspace(x[0]-dx*(xgauss.shape[0]-1)/2, x[-1]+dx*(xgauss.shape[0]-1)/2, x.shape[0]+xgauss.shape[0]-1)
    vgt = np.convolve(Gauss(xgauss,1.0, 0.0, sig), Lorentz(x_n,1.0,x0,g), mode = 'valid')
    return A*vgt/np.max(vgt)

def MPeak(x, *p, ptype = 'voigt', N_po = -1, N_pe = 1, sbg = False):
    # order for p is: 
    # (a0,a1,a2,..., p0A, p0x0, p0g, p0sig, p1A, p1x0, ...)
    if type(x) == np.ndarray:
        y = np.zeros(x.shape)
    else:
        y = 0.0
    for i in range(N_po+1): #background polynomial
        y += p[i]*x**i

    for i in range(N_pe):
        if ptype == 'voigt':
            y += Voigt(x, *p[N_po+1 + 4*i: N_po+1 + 4*(i+1)])
        elif ptype == 'gauss':
            y += Gauss(x, *(p[N_po+1 + 4*i: N_po+1 + 4*i + 2] + p[N_po+1 + 4*i + 3: N_po+1 + 4*i + 4]))   
        elif ptype == 'lorentz':
            y += Lorentz(x, *p[N_po+1 + 4*i: N_po+1 + 4*i + 3])
            
    if sbg:
        sbg_fac = p[-1]
        y_nobg = MPeak(x, *p[N_po+1:], ptype = ptype, N_po = -1, N_pe = N_pe, sbg = False)
        y += get_SBG(y_nobg, x, sbg_fac)            
            
    return y

    
    
def MPeakff(p, x, data, ptype = 'voigt', N_po = -1, N_pe = 1, sbg = False):
    
    pp = ()
    for i in range(N_po + 1):
        pp += (p['a%s' % i].value,)

    for i in range(N_pe):
        pp += (p['A%s' % i].value,)
        pp += (p['x0%s' % i].value,)
        pp += (p['g%s' % i].value,)
        pp += (p['sig%s' % i].value,)

    if sbg:
        pp += (p['sbg_fac'].value,)


    return MPeak(x, *pp, ptype = ptype, N_po = N_po, N_pe = N_pe, sbg = sbg) - data
    

def get_SBG(spec, x, sbg_fac):
    """Return shirley background for SBG, multiply by factor a"""
    return sbg_fac*inte.cumtrapz(spec[::-1], x = x, initial = 0)[::-1]

    



def MPeakGuess(x, data, N_po = 1, N_pe = 1, w = None, ptype = 'voigt', sbg = False):
    #N is the number of peaks to look for, w is the minimum width of peaks to look for

    if w is None:
        w = 0.05*abs(x[-1]-x[0]) #assume 5% of window as smoothing
        
    dx = x[1]-x[0]
    w_n = int(5*w/dx)
    xgauss = np.linspace(-w_n*dx,w_n*dx,2*w_n + 1) #minimal 1 pt, same dx scale.
    gauss = Gauss(xgauss, 1.0, 0.0, w)
    data_s = fftconvolve(data, gauss, 'same')/np.sum(gauss) #smoothed data
    
    data_s_mask = data_s[int(w/dx):-int(w/dx)] #exclude points near boundary
    x_mask = x[int(w/dx):-int(w/dx)]
    #get linear BG:
    lbg = np.linspace(data_s_mask[0],data_s_mask[-1],data_s_mask.size)
    
    x0 = np.average(lbg)
    x1 = (lbg[1]-lbg[0])/(x_mask[1]-x_mask[0])
    data_s_mask -= lbg
    
    p = ()
    for i in range(N_pe):
        peak_i = argrelextrema(data_s_mask, np.greater)[0] # find peak maxima
        peak_A = data_s_mask[peak_i] # find values
        peak_i = peak_i[np.argsort(peak_A)[::-1]][0:N_pe] # sort by highest max
        peak_A = peak_A[np.argsort(peak_A)[::-1]][0:N_pe] # same
        #argrelextrema does not give boundary points
        
        #get width at half max:
        halfmax_i = np.where(np.diff(np.signbit(data_s_mask-peak_A[0]/2)))[0]
        halfmax_x = x_mask[halfmax_i]
        #find closes zerocrossing negative:
        xdist = halfmax_x-x_mask[peak_i[0]]
        xnegs = xdist[np.argwhere(xdist < 0)]
        xneg = xnegs.min()
        xposs = xdist[np.argwhere(xdist > 0)]
        xpos = xposs.max()
        
        w = xpos - xneg
        
        if ptype == 'voigt':
            p += (peak_A[0], x_mask[peak_i[0]], w/np.sqrt(2), w/np.sqrt(2)) #share widht between gauss and lorentz
        elif ptype == 'gauss':
            p += (peak_A[0], x_mask[peak_i[0]], 0.0, w) #weight all gauss
        elif ptype == 'lorentz':
            p += (peak_A[0], x_mask[peak_i[0]], w, 0.0) #weight all lorentz
                  
        data_s_mask -=  MPeak(x_mask, *p[-4:], ptype = ptype, N_po = -1, N_pe = 1)
    
    if sbg:
        p += (1.,)
            
    ppo = ()
    for i in range(N_po+1):
        if i == 0:
            ppo += (x0,)
        elif i == 1:
            ppo += (x1,)
        else:
            ppo += (0.0,) #add 0 as first app for higher powers

    return ppo + p

def MPeakFit(x, data, N_pe = None, N_po = -1, ptype = 'voigt', w = None, pguess = None, 
             constraints = {}, sbg = False, maxfev = 400, plot = False, ax = None):
    """Multi peak fit, N_pe is number of peaks, N_po is polynomial background order
    w is approximate width, ptype is the type of peak"""
        
    if pguess is None:
        pguess =  MPeakGuess(x, data, N_po = N_po, N_pe = N_pe, w = w, ptype = ptype, sbg = sbg)
    

    
    p = Parameters()

    for i in range(N_po+1):
        #for constraints: set a default, and if a set is given in "constraints" dict, overwrite whichever are given
        #that means defaults stay in place unless you specifically overwrite them
        sub_constraints = {'vary' : True}
        if ('a%s' % i) in constraints:
            sub_constraints = dict(sub_constraints, **constraints[('a%s' % i)])
        p.add(('a%s' % i), value=pguess[i], **sub_constraints)
        
    for i in range(N_pe):
        #amplitude:
        sub_constraints = {'vary' : True, 'min' : 0.}
        if ('A%s' % i) in constraints:
            sub_constraints = dict(sub_constraints, **constraints[('A%s' % i)])
        p.add(('A%s' % i), value=pguess[4*i + N_po + 1 + 0], **sub_constraints)
        
        #x0:
        sub_constraints = {'vary' : True}
        if ('x0%s' % i) in constraints:
            sub_constraints = dict(sub_constraints, **constraints[('x0%s' % i)])
        p.add(('x0%s' % i), value=pguess[4*i + N_po + 1 + 1], **sub_constraints)
        
        
        #lorentz width:
        sub_constraints = {'vary' : ((ptype == 'voigt') or (ptype == 'lorentz')), 'min' : 0.0} # only vary for lorentz and voigt
        if ('g%s' % i) in constraints:
            sub_constraints = dict(sub_constraints, **constraints[('g%s' % i)])
        p.add(('g%s' % i), value=pguess[4*i + N_po + 1 + 2], **sub_constraints) 

        #gauss width:
        sub_constraints = {'vary' : ((ptype == 'voigt') or (ptype == 'gauss')), 'min' : 0.0}  #only vary for gauss and voigt
        if ('sig%s' % i) in constraints:
            sub_constraints = dict(sub_constraints, **constraints[('sig%s' % i)])        
        p.add(('sig%s' % i), value=pguess[4*i + N_po + 1 + 3], **sub_constraints) #only vary for gauss and voigt
        
        
        if sbg:
            sub_constraints = {'value' : 3.0, 'min' : 0., 'vary' : True}
            if 'sbg_fac' in constraints:
                sub_constraints = dict(sub_constraints, **constraints['sbg_fac'])        
            
            p.add('sbg_fac', **sub_constraints)
              
    kws = {'ptype' : ptype, 'N_po' : N_po, 'N_pe' : N_pe, 'sbg' : sbg}

    result = minimize(MPeakff, p, args=(np.array(x, dtype = 'float64'), np.array(data, dtype = 'float64')), kws = kws, maxfev = 2*maxfev)

    pp = ()
    for i in range(N_po+1):
        pp += (result.params['a%s' % i].value,)

    for i in range(N_pe):
        pp += (result.params['A%s' % i].value,)
        pp += (result.params['x0%s' % i].value,)
        pp += (result.params['g%s' % i].value,)
        pp += (result.params['sig%s' % i].value,)
    if sbg:
        pp += (result.params['sbg_fac'].value,)
        
    if plot and (result.nfev <= maxfev):
        MPeakPlot(x, data, pp, dict(N_pe = N_pe, N_po = N_po, ptype = ptype, sbg = sbg), ax = ax, indiv_pks = True)
        
    if result.nfev > maxfev:
        return None
    else:
        return pp, dict(N_pe = N_pe, N_po = N_po, ptype = ptype, sbg = sbg)
        
        
        
    
def MPeakPlot(x, data, pp, ppkwarg, ax = None, indiv_pks = True):
    """Plot the output of Mpeakfit"""
    
    ppkwarg_std = dict(sbg = False)
    ppkwarg = dict(ppkwarg_std, **ppkwarg)
    
    if ax is None:
        fig,ax = plt.subplots(figsize = (4,4))
        

    ax.plot(x, data, color = 'black')
    if indiv_pks:
        ax.plot(x, MPeak(x, *pp, **ppkwarg), color = 'red')
        for i in range(ppkwarg['N_pe']):
            pk = pp[ppkwarg['N_po']+1+4*(i+0):ppkwarg['N_po']+1+4*(i+1)]
            ax.plot(x, MPeak(x, *pk, ptype = ppkwarg['ptype'], N_pe = 1, N_po = -1)+np.min(data), color = 'green')
        if ppkwarg['sbg']:            
            ax.plot(x, get_SBG(MPeak(x, *pp[ppkwarg['N_po']+1:], ptype = ppkwarg['ptype'], N_po = -1, N_pe = ppkwarg['N_pe'], sbg = False), x, pp[-1]), color = 'orange')
    else:
        ax.plot(x, MPeak(x, *pp, **ppkwarg), color = 'red')

            
    #pass individual peaks without background in green
    ax.axis([x[0],x[-1],np.min(data),np.max(data) + 0.2*(np.max(data)-np.min(data))])
    
    
    
def UnpackMPeak(fit_result_list, add_widths = True):
    """Unpacks results for a list of MPeak results
    Args:
        result_list: list of the form [(N_pe, N_po, ptype, pp), ...]
    Returns:
        a list of np arrays as (c0,c1,...,A0,x0,fwhm0,... or if add_widths is false:
        a list of np arrays as (c0,c1,...,A0,x0,g0,sig0,... """
        
    returnlist = []
    N_po = fit_result_list[0][1]
    N_pe = fit_result_list[0][0]     
    for i in range(N_po+1):
        returnlist.append(np.array([res[3][i] for res in fit_result_list]))
                 
    for i in range(N_pe): 
        returnlist.append(np.array([res[3][N_po+1+4*i+0] for res in fit_result_list]))
        returnlist.append(np.array([res[3][N_po+1+4*i+1] for res in fit_result_list]))
        if add_widths:
            lorwidth = np.array([res[3][N_po+1+4*i+2] for res in fit_result_list])
            gausswidth = np.array([res[3][N_po+1+4*i+3] for res in fit_result_list])
            returnlist.append(0.5346*lorwidth + np.sqrt(0.2166*lorwidth**2 + gausswidth**2))
        else:
            returnlist.append(np.array([res[3][N_po+1+4*i+2] for res in fit_result_list]))
            returnlist.append(np.array([res[3][N_po+1+4*i+3] for res in fit_result_list]))
               
    return returnlist

            

def FitBands():
    """
    """
    pass


def Poly(x, *p, N_po = 0):
    if type(x) == np.ndarray:
        y = np.zeros(x.shape)
    else:
        y = 0.0
    for i in range(N_po+1):
        y += p[i]*x**i
    return y
    
def Polyff(p, x, data, N_po = 0):
    
    pp = ()
    for i in range(N_po + 1):
        pp += (p['a%s' % i].value,)

    return Poly(x, *pp, N_po = N_po) - data


def PolyFit(x, data, N_po = 2):
    pguess = PolyGuess(x,data,N_po)
    
    p = Parameters()

    for i in range(N_po+1):
        p.add(('a%s' % i), value=pguess[i], vary=True)    
    kws = {'N_po' : N_po}

    result = minimize(Polyff, p, args=(np.array(x, dtype = 'float64'), np.array(data, dtype = 'float64')), kws = kws)

    pp = ()
    for i in range(N_po+1):
        pp += (result.params['a%s' % i].value,)

        

    return N_po, pp
    
def PolyGuess(x,data,N_po):        
    return (0,)*(N_po+1)

    
class MDCfitter():
    """Class used for MDC fits of 2D data
    to do:
        update fitting so that it has a "results" attribute
        update it so that it has a "fit_one" method
        update it so that it has a "show_one" method
        update methods to work with new style data saving
        update fit so that you can choose between "makeguess, giveguess, giveseed, useseed" 
        update AaData2D so that you have a method that returns a kaxis for y, so no need to transform the data
    """
    def __init__(self, obj2D, fittype = 'a', fitrange = None, rebin = (1,1), ptype = 'lorentz', N_po = 1, N_pe = 1):
        """Class initializer
        
        Args:
            obj2D: the AaData2D object to be fit
            fittype: string, either 'a' or 'k' selecting raw data or k-corrected data
            fitrange: 4 element list or tuple indicating the fitrange as (xmin, xmax, ymin, ymax)
            ptype: the peak type, 'lorentz', 'gauss' or 'voigt'
            N_po: order of the polynomial, -1 for no polynomial background
            N_pe: number of peaks
        """
        self.obj2D = obj2D
        self.fittype = fittype
        if fitrange is None:
            if self.fittype == 'a':
                self.xmin = self.obj2D.xaxis[0]
                self.xmax = self.obj2D.xaxis[-1]
                self.ymin = self.obj2D.yaxis[0]
                self.ymax = self.obj2D.yaxis[-1]
            else:
                print('Not implemented yet')
        else:
            self.xmin, self.xmax, self.ymin, self.ymax = fitrange
        
        self.fitobj = self.obj2D.crop(self.xmin, self.xmax, self.ymin, self.ymax).rebin(rebin)
        
        self.setpeaks(N_po = N_po, N_pe = N_pe, ptype = ptype)

    def setpeaks(self, N_po = None, N_pe = None, ptype = None):
        """setup the polynomial order, number of peaks and peaktype for the fit
        
        Args:
            N_po: polynomial background order, -1 for no background
            N_pe: number of peaks in the fit
            ptype: peaktype to be fit, 'gaussian' 'voigt' or 'lorentz'
        Returns:
            None
        """
        if N_po is not None:
            self.N_po = N_po
        if N_pe is not None:
            self.N_pe = N_pe
        if ptype is not None:
            self.ptype = ptype
            
        #fits is an array length of number of MDCs, width is [a0, a1, a2, ..., an, A0, x0, sl0, sg0, A1, x1, sl1, sg1, ...] 
        #an are polynomial coefs, A0, x0, sl0, sg0 are amplitude, location, lorentzian and gaussian width of the paeks
        self.fits = np.zeros((self.fitobj.xaxis.size, self.N_po + 1 + 4*self.N_pe))
        self.amplitudes = np.zeros((self.fitobj.xaxis.size, self.N_pe))
        self.positions = np.zeros((self.fitobj.xaxis.size, self.N_pe))
        self.lorwidths = np.zeros((self.fitobj.xaxis.size, self.N_pe))
        self.gausswidths = np.zeros((self.fitobj.xaxis.size, self.N_pe))
        self.fwhms = np.zeros((self.fitobj.xaxis.size, self.N_pe))
        
        
    def fitone(self, n, pguess = None, w = None, showfit = True):
        """Fit one MDC, which can be used as a seed for fitting others
        
        Args:
            n: the number of the MDC to fit
            w: approximate width of the peaks for filtering and peak detection (if None, w is set to 10% of the window)
            pguess: guess parameters (if None, they are generated automatically with w
            showfit: boolean, if True it shows the result of the fit
        Returns:
            None
        """
        self.fits[n] = MPeakFit(self.fitobj.yaxis, self.fitobj.data[n], N_pe = self.N_pe, N_po = self.N_po, ptype = self.ptype, w = w, pguess = pguess)[3]
        self.show_n(n)
        
        
    def unpackresults(self):
        
        for i in range(self.N_pe):                  
            self.amplitudes[:,i] = self.fits[:,self.N_po + 1 + 4*i + 0]
            self.positions[:,i] = self.fits[:,self.N_po + 1 + 4*i + 1]
            self.lorwidths[:,i] = self.fits[:,self.N_po + 1 + 4*i + 2]
            self.gausswidths[:,i] = self.fits[:,self.N_po + 1 + 4*i + 3]
            self.fwhms = 0.5346*self.lorwidths + np.sqrt(0.2166*self.lorwidths**2 + self.gausswidths**2)
        
        
        
    def fit(self, showfit = True, w = None, mode = 'guess', pguess = None, seed = 0):
        """Fit the MDCs as defined in the class
        
        Args:
            showfit: show results at the end
            w: approximate width of the peaks for guessing, if None it is set to 10% of the window
            mode: the mode to use, modes are 'guess' (runs the guesser for every MDC with w, or if pguess is not None it uses pguess for every MDC)
                                            'seed' (runs using fitresults from mdc number in seed, then uses results as a guess for the next fit, runs to the postive and negative end)
                                            'seedup' (see seed, only runs from n+1 to the end)
                                            'seeddown' (see seed, only runs from n-1 to zero)
            pguess: tuple of initial parameters
            seed: MDC number to use as a seed
        """
        if self.fittype == 'a':
            if mode == 'guess':
                for i in range(self.fitobj.data.shape[0]):
                    self.fits[i] = MPeakFit(self.fitobj.yaxis, self.fitobj.data[i], N_pe = self.N_pe, N_po = self.N_po, ptype = self.ptype, w = w, pguess = pguess)[3]
            if (mode == 'seedup') or (mode == 'seed'):
                for i in range(seed+1, self.fitobj.data.shape[0]):
                    self.fits[i] = MPeakFit(self.fitobj.yaxis, self.fitobj.data[i], N_pe = self.N_pe, N_po = self.N_po, ptype = self.ptype, pguess = tuple(self.fits[i-1]))[3]
            if (mode == 'seeddown') or (mode == 'seed'):
                for i in range(seed-1, -1,-1):
                    self.fits[i] = MPeakFit(self.fitobj.yaxis, self.fitobj.data[i], N_pe = self.N_pe, N_po = self.N_po, ptype = self.ptype, pguess = tuple(self.fits[i+1]))[3]

        self.unpackresults()
        self.showfit()
            
    def getresult(self,n, plot = True):
        """Get the result of the nth fit, optionally show a plot of said fit"""
        print("\nMDC fit at Energy %s " % self.fitobj.xaxis[n])
        print("")
        print("Peak type: %s" % self.ptype)
        print("Number of peaks: %s" % self.N_pe)
        print("Background polynomial order: %s\n\n" % self.N_po)
        for i in range(self.N_pe):
            print("Peak %s Amplitude: %s" % (i,self.fits[n,self.N_po + 1 + (4*i) + 0]))
            print("Peak %s Location: %s" % (i,self.fits[n,self.N_po + 1 + (4*i) + 1]))
            print("Peak %s Lorentzian width: %s" % (i,self.fits[n,self.N_po + 1 + (4*i) + 2]))
            print("Peak %s Gaussian width: %s" % (i,self.fits[n,self.N_po + 1 + (4*i) + 3]))

        if plot:
            self.show_n(n)
            
            
            
    def show_n(self,n, peaks = True):
        """Show the nth fit in a plot, overlaid with data, optionally with separate peaks"""
        fig,ax = plt.subplots(figsize = (4,4))
        
        pp_n = tuple(self.fits[n])

        ax.plot(self.fitobj.yaxis, self.fitobj.data[n]-np.min(self.fitobj.data[n]), color = 'black')
        ax.plot(self.fitobj.yaxis, MPeak(self.fitobj.yaxis, *pp_n, ptype = self.ptype, N_pe = self.N_pe, N_po = self.N_po)-np.min(self.fitobj.data[n]), color = 'red')
        

        
        for i in range(self.N_pe):
            pk = pp_n[self.N_po+1+4*(i+0):self.N_po+1+4*(i+1)]
            ax.plot(self.fitobj.yaxis, MPeak(self.fitobj.yaxis, *pk, ptype = self.ptype, N_pe = 1, N_po = -1), color = 'green')
            
        ax.axis([self.fitobj.yaxis[0],self.fitobj.yaxis[-1],0,np.max(self.fitobj.data[n])*1.2])
        #pass individual peaks without background in green
        
        
    def showfit(self, plotrange = None, cstart = None, cend = None, colors = None, nMDCfits = 5, offset = 1.5):
        """show a fit overview, fits overlaid with data plus a sample of the fits"""
        
        fig,ax = plt.subplots(1,2,figsize = (7,4)) #show fits overlayed on spec on left, show MDCs + fits on right
        
        if cstart is None:
            cstart = np.min(self.obj2D.data)
        if cend is None:
            cend = np.max(self.obj2D.data)
        
        if plotrange is None:
            plotrange = [self.obj2D.yaxis[0],self.obj2D.yaxis[-1],self.obj2D.xaxis[0],self.obj2D.xaxis[-1]]
        
        ax[0].pcolormesh(*np.meshgrid(self.obj2D.yaxis,self.obj2D.xaxis), self.obj2D.data, cmap = 'bone_r', vmin = cstart, vmax = cend)
        if colors is None:
            colors = (int(np.ceil(self.N_pe/4.0))*['red','black','blue','green'])[:self.N_pe]
        for pos,fwhm,c in zip(self.positions.T,self.fwhms.T,colors):
            ax[0].errorbar(pos,self.fitobj.xaxis, xerr = fwhm/2.0, color = c, marker = 's', linestyle = 'none', alpha = 0.3)
        ax[0].axis(plotrange)
        
        fitEs = np.linspace(self.fitobj.xaxis[0],self.fitobj.xaxis[-1], nMDCfits)
        
        delta = 0.0
        for E in fitEs:
            n = np.absolute(self.fitobj.xaxis-E).argmin()
            ax[1].plot(self.fitobj.yaxis, delta + self.fitobj.data[n], color = 'black')
            ax[1].plot(self.fitobj.yaxis, delta + MPeak(self.fitobj.yaxis, *tuple(self.fits[n]), ptype = self.ptype, N_pe = self.N_pe, N_po = self.N_po), color = 'red')
            delta += offset*(np.max(self.fitobj.data[n])-np.min(self.fitobj.data[n]))
            
        ax[1].axis([self.fitobj.yaxis[0],self.fitobj.yaxis[-1],np.min(self.fitobj.data),np.max(delta + self.fitobj.data)])
        