# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 09:48:21 2016

@author: berend
"""
import numpy as np
from lmfit import minimize, Parameters#, Parameter
from scipy.interpolate import interp1d


##NOTE ABOUT MINIMIZE:
#There is a bug that if you fit a dataset that is float32 it will not fit

def FD(x, EF, T): #24.8 us per loop 
    """Fermi Dirac distribution function EF Fermi energy in eV, T temp in K"""
    return 1./(1+np.exp((x-EF)/(T*8.61733E-5)))
    
def FD2(x, *params):
    """FD distr, with EF Fermi energy, T temperature, SIG resolution taken into
    account as effective temperature DY1 slope below EF, DY2 slope above EF, C 
    offset"""
    (EF, T, sig, A, dy1, dy2, c) = params
    Teff = np.sqrt(T*T + 8*8.41646E6*sig*sig)
    
    slope1 = 1. + dy1*(x-EF)
    slope2 = dy2*(x-EF)
    
    return slope1*A*FD(x,EF,Teff)+slope2+c
    

def FD2ff(p, x, data):#fitfunc
    return FD2(x,p['EF'],p['T'],p['sig'],p['A'],p['dy1'],p['dy2'],p['c']) - data

    
def FD2p(x, data, T): #generate parameters for FD2
    params = FD2Guess(x,data,T)
    p = Parameters()
    p.add_many(('EF',   params[0],  True,   x[0],   x[-1],  None),
               ('T',    params[1],  False,  None,   None,   None),
               ('sig',  params[2],  True,   0.0,    None,   None),
               ('A',    params[3],  True,   0.0,    None,   None),
               ('dy1',  params[4],  True,   None,   None,   None),
               ('dy2',  params[5],  False,  None,   None,   None),
               ('c',    params[6],  True,   None,   None,   None))
    return p
    
def FD2fit(x,data,T):
    p = FD2p(x,data,T)
    result = minimize(FD2ff, p, args=(np.array(x, dtype = 'float64'), np.array(data, dtype = 'float64')))
    EF = result.params['EF'].value
    T = result.params['T'].value
    sig = result.params['sig'].value
    if((EF < x[0]) or (EF > x[-1]) or (sig > 0.5*(x[-1]-x[0]))):
        raise RuntimeError('Fit quality inadequate')
    A = result.params['A'].value
    dy1 = result.params['dy1'].value
    dy2 = result.params['dy2'].value
    c = result.params['c'].value
    return (EF,T,sig,A,dy1,dy2,c)

def FD2Guess(x, y, T):
    """Make a guess for FD2 parameters from data in Y and X with temp T"""
    c = y[-1]
    dy1 = 0
    dy2 = 0
    A = y.max()
    sig = 0.025 #fixed guess, should yield right value
    EF = x[np.argmin(np.absolute(y-y.max()*0.5))]
    ## incorp error handling here so faulty ones don't get fit. Don't know what error could be here? EF cannot be outside the region
    #only sig can be wrong but that is taken care off statically and through Find region
    return (EF, T, sig, A, dy1, dy2, c)
  


def FPoly(x, *params):
    Np = params[0] #Np is the poly order
    cn = params[1:] #the rest of the params (Np+1) are the coefs
    
    y = np.array(x)
    
    for ii in range(Np + 1):
        y += np.power(x,ii)*cn[ii]
    
    return y


def FPolyp(x, data, Np): #generate parameters for FD2
    params = FPolypGuess(x, data, Np)
    p = Parameters()
    p.add('Np', value=Np, vary=False)

    for ii in range(Np + 1):
        p.add('c%s' % ii, value = params[ii+1], vary=True)

    return p


def FPolyff(p, x, data):#fitfunc
    Np = p['Np'].value
    params = [Np]
    for ii in range(Np+1):
        params.append(p['c%s' % ii].value)

    params = tuple(params) #convert back to tuple
    
    return FPoly(x, *params) - data

def FPolypGuess(x, data, Np):
    params = [Np]
    
    #guess c0 to be the average, the rest 0
    
    for ii in range(Np+1):
        if(ii == 0):
            params.append(np.average(data))
        else:
            params.append(0)
    
    return tuple(params)
    

def FPolyfit(x,data,Np):
    p = FPolyp(x,data,Np)
    result = minimize(FPolyff, p, args=(np.array(x,dtype = 'float64'), np.array(data, dtype = 'float64')))

    Np = result.params['Np'].value
    result_p = [Np]
    for ii in range(Np+1):
        result_p.append(result.params['c%s' % ii].value)

    result_p = tuple(result_p) #convert back to tuple

    return result_p
  
## Error handling done
def AnalyzeGold(x, data, T, Np = 2, region = None):
    """Find a suitable region to fit, then use this region to fit gold
    then fit 3rd order poly to the spectrum, return average Ek and polynomial"""
    if region == None:
        try:
            region = FindRegion(x, data)
        except:
            print('Unable to find region, please provide suitable fit region')
    elif region == 'full':
        region = (0,data.shape[0],0,data.shape[1])
            
    #make results dictionary
    fititems = region[3]-region[2]
    amin_n = region[2]
    Emin_n = region[0]
    Emax_n = region[1]
    fitresults = []
    
    print('Found region (%s,%s,%s,%s)' % region)

    
    #do the fit of the gold spectra (use try statements here as well?)
    for ii in range(fititems):
        try:
            #save both the column number and the results if fitting succeeds
            fitresults.append((ii+amin_n, *FD2fit(x[Emin_n:Emax_n], data[Emin_n:Emax_n,ii+amin_n], T)))
        except:
            print('Failed fit at %s' % (ii+amin_n))
            pass
        
    if(len(fitresults) < 0.8*fititems):
        raise RuntimeError('Review region, too many failed fits')
    
    #convert to numpy array
    fitresults = np.array(fitresults, dtype = 'float')
    
    fitdict = {}

    fitdict['col'] = np.array(fitresults[:,0])
    fitdict['EF'] = np.array(fitresults[:,1])
    fitdict['T'] = np.array(fitresults[:,2])
    fitdict['sig'] = np.array(fitresults[:,3])
    fitdict['A'] = np.array(fitresults[:,4])
    fitdict['dy1'] = np.array(fitresults[:,5])
    fitdict['dy2'] = np.array(fitresults[:,6])
    fitdict['c'] = np.array(fitresults[:,7])

    #then fit a poly to col vs EF
    
    EFparam = FPolyfit(fitdict['col'], fitdict['EF'], Np)
    
    fitEdge = FPoly(np.arange(data.shape[1], dtype = 'float'), *EFparam) #make the curve to return

    return fitdict, EFparam, fitEdge
    
def AnalyzeGold1D(x, data, T, region = None):

    if region is None:
        region = (0,x.shape[0]) #middle fifth
    
    
    
    fitres = FDexpfit(x[region[0]:region[1]], data[region[0]:region[1]], T)
    
    print('Resolution is %s mV' % (1000*fitres[2]))
    
    return fitres    
    
    

def GoldCorrect(data, xaxis, Edge):
    """" t is transformation matrix"""
    out = np.array(data)
    xaxisout = np.array(xaxis)

    #dterine EF
    EF = np.average(Edge)
    
    #correct xaxis:
    xaxisout = xaxis - EF
    
    
    for i in range(out.shape[1]):
        interp = interp1d((xaxis-Edge[i]), data[:,i],bounds_error = False, fill_value = 0.0)
        out[:,i] = interp(xaxisout) #  actually interp returns an array with x_new and y_new as input arays interp(array_x, array_y)
            
    return out, xaxisout   
    
def Gold1DCorrect(data, xaxis, EF):
    out = np.array(data)
    xaxisout = np.array(xaxis)    
    #correct xaxis:
    xaxisout = xaxis - EF 
    # transformation for 1D is trivial
            
    return out, xaxisout       
        

def FindRegion(x, data):
    """Find a suitable region to fit"""
    #raise errors so other functions can make use of it
    #start in the middle:
    sample = np.sum(data[:,512-10:512+10]/20, axis = 1) #integrate a bit
    EF_n = np.argmin(np.absolute(sample-sample.max()*0.3))
    EF = x[EF_n]
    Emax = x[-1]
    deltaE = Emax - EF
    Erange = min(deltaE, (0.1*(x[-1]-x[0]))) # the fitrange should be the EF+/- 0.1 total range, unless the range is smaller than that
    dx = x[1]-x[0]
    Emax_n = EF_n + int(Erange/dx)
    Emin_n = EF_n - int(Erange/dx)
    if((Emax_n-Emin_n) < 100):
        raise RuntimeError('Can\'t find energy range')
    #then find a cut at the standard 0.025 resolution below the fermi level
#    a_ns = np.argsort(np.absolute(data[EF_n-int(0.025/dx)]-0.1*data.max()))#get size in angles
#    amin_n = a_ns[0]
#    i = 1
#    amax_n = a_ns[i]
#    while abs(amax_n-amin_n) < 200: #make sure amin and amax are at least 200 apart
#        amax_n = a_ns[i]
#        i += 1
#    if(amin_n > amax_n):
#        temp = amin_n
#        amin_n = amax_n
#        amax_n = temp
    amin_n = 0
    amax_n = data.shape[1]

    if((amax_n - amin_n) < 300): #if the range found was less than 200 the data was probably too noisy to make a good estimate
        raise RuntimeError('Unable to find suitable angle range')
    
    return (Emin_n, Emax_n, amin_n, amax_n)

    
def FitRegion(x, data, T):
    """Fit series of gold specs, X being the x axis of DATA"""
    paramlist = []
    anglelist = []
    for i in range(data.shape[1]):
        print(i)
        try:
            fitparams = FD2fit(x, data[:,i],T) #figure out T, also figure out proper error handling
            paramlist.append(fitparams)
            anglelist.append(i)
        except RuntimeError as a:
            print('For slice i = %s' % i)
            print(a)
    return np.array(paramlist),np.array(anglelist)

    
    
    


def Gauss(x, A, x0, sig):
    """Gaussian function, sig is the FWHM"""
    Sig = sig/(2*np.sqrt(2*np.log(2)))
    return A * np.exp(-(x-x0)**2 /(2*Sig*Sig))

    
    
def FDexp(x, EF, T, sig, A = 1, dy1 = 0, dy2 = 0, c = 0):
    """Fermi Dirac fit function using a convolution"""
    dx = x[1]-x[0]
    sig_n = int(5*sig/dx)
    xgauss = np.linspace(-sig_n*dx,sig_n*dx,2*sig_n + 1) #minimal 1 pt, same dx scale. Lenght is N, need to extend x by (N-1)/2 on each side
    x_n = np.linspace(x[0]-dx*(xgauss.shape[0]-1)/2, x[-1]+dx*(xgauss.shape[0]-1)/2, x.shape[0]+xgauss.shape[0]-1)
    return np.convolve(Gauss(xgauss,1.0, 0.0, sig), FD2(x_n,EF,T,0,A,dy1,dy2,c), mode = 'valid')/np.sum(Gauss(xgauss,1.0, 0.0, sig))



def FDexpff(p, x, data):#fitfunc
    return FDexp(x,p['EF'],p['T'],p['sig'],p['A'],p['dy1'],p['dy2'],p['c']) - data

    
def FDexpp(x, data, T, guess = None): #generate parameters for FD2
    if not guess:
        params = FD2Guess(x,data,T) #use FD2 guess function, should work fine
    else:
        params = guess
    p = Parameters()
    p.add_many(('EF',   params[0],  True,   x[0],   x[-1],  None),
               ('T',    params[1],  False,  None,   None,   None),
               ('sig',  params[2],  True,   0.0,    None,   None),
               ('A',    params[3],  True,   0.0,    None,   None),
               ('dy1',  params[4],  True,   None,   None,   None),
               ('dy2',  params[5],  False,  None,   None,   None),
               ('c',    params[6],  True,   None,   None,   None))
    return p
    
def FDexpfit(x,data,T):
    p = FDexpp(x,data,T)
    result = minimize(FDexpff, p, args=(np.array(x, dtype = 'float64'), np.array(data, dtype = 'float64')))
    EF = result.params['EF'].value
    T = result.params['T'].value
    sig = result.params['sig'].value
    if((EF < x[0]) or (EF > x[-1]) or (sig > 0.5*(x[-1]-x[0]))):
        raise RuntimeError('Fit quality inadequate')
    A = result.params['A'].value
    dy1 = result.params['dy1'].value
    dy2 = result.params['dy2'].value
    c = result.params['c'].value
    return (EF,T,sig,A,dy1,dy2,c)
