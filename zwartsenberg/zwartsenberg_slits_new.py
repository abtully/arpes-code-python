# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 10:08:24 2016

@author: berend
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.signal import argrelextrema


#Slits work well for 12 eV pass energy, WAM. Not so well for laser data or 
#LAD data, have not tried HeII data yet (guessing the same problem)
#Approaches: 
#1. incorporate a script that filters out curves that are intersecting for the mask
#2. possibly incorporate a "mask" for ranges of kinetic energies
#3. make different approaches for LAD and laser data

def AnalyzeSlits(slits, N = None):
    """Analyze slits"""
    
    #determine N
    if(N == None):
        N = GuessPeakNum(slits)
    
    thresh = 0.05 #threshold for digitizing slits
    mask_x, mask_peak, mask_antipeak = FindRegions(digitize(slits, thresh), N = N)
    maskparams = FitSlitCurves(mask_antipeak, mask_x, N_p = 1)
    x = np.arange(slits.shape[0])
    N_p = 2 # slit poly order
    curveparams = np.zeros((N, N_p+1))
    slitcurves = np.zeros((N, slits.shape[0]))
    print('Mask found')
    for i in range(N):
        y_lo = Poly(x, *maskparams[i])
        y_hi = Poly(x, *maskparams[i+1])
        y_curve, x_curve = FindSlitCurve(slits, y_lo, y_hi)
        curveparam = FitSlitCurve(y_curve, N_p, x_curve)
        curveparams[i,:] = curveparam
        slitcurves[i] = Poly(x,*curveparam)
        print('Curve %s Fit' % i)
    T,x,y = MakeTransformationMatrix(curveparams, slits.shape)
    return T, x, y, curveparams, slitcurves
        


def digi(val, thresh):
    if(val < thresh):
        return 0
    else: return 1

    
def digitize(A, thresh):
    """return 1 if A[i,j] greater then thresh*A.max(), else return 0"""
    Thresh = thresh * A.max()
    return np.where(np.greater(A, Thresh),1,0)




def FindRegions(data, N = 17, sep = 10, n = 20):
    """Find N regions in data with minimal separation sep (in pixels). Input
    should be digitized (0 or 1 input). return lines that are fitted to the
    middle between the peaks."""
    i = 0
    peaks = [] # list with peak positions in the line
    antipeaks = [] # list with positions between the peaks
    xvalues = [] #list with array positions of where the line was taken
    while True:
        try:
            xvalue,peak, antipeak = FindDigPeaksCluster(data[n*i:n*(i+1)], N = N, sep = sep)
            peaks.append(peak)
            antipeaks.append(antipeak)
            xvalues.append(i*n + xvalue) # append position in the array
        except RuntimeError:
            pass
        i += 1
        if(i*(n+1) >= data.shape[0]):
            break
    return np.array(xvalues),np.array(peaks).transpose(), np.array(antipeaks).transpose()


def FindDigPeaks(data, N, sep):
    """Find locations of peaks in 1D array data, minimal separation sep (pix)
    raise runtime error if not enough peaks found"""
    level = 0.5 * data.max()
    risexlist = []
    fallxlist = []
    peakpos = []
    antipeakpos = []

    for i in range(len(data)-1): # find all peaks
        if((data[i]<level)&(data[i+1]>level)):
            risexlist.append(i)
        if((data[i]>level)&(data[i+1]<level)):
            fallxlist.append(i)   
    #check if they're further apart then separation
    i = 0
    
    while True:    
        if(i >= (len(risexlist)-1)):
            break
        if((risexlist[i+1] - fallxlist[i]) < sep): #if they're too close together
            risexlist.pop(i+1)
            fallxlist.pop(i)    #remove peaks that are too close together
        else:
            i += 1


        
    #check if right number of peaks    
    if(len(risexlist) != N):
        raise RuntimeError('Not the right number of peaks')
    
    for i in range(len(risexlist)):
        peakpos.append((fallxlist[i]+risexlist[i])/2.) # calc mid peakposition
        if(i != 0): antipeakpos.append((fallxlist[i-1]+risexlist[i])/2.) # and middle between the peak
    
    #extend antipeakpos:
    antipeakpos.insert(0,2*peakpos[0]-antipeakpos[0])
    antipeakpos.append(2*peakpos[-1]-antipeakpos[-1])    
    
    return peakpos, antipeakpos

def FindDigPeaksCluster(data, N, sep):
    """Find the first line in the array where FindDigPeaks does not error"""
    """ return the positions of the first line that fits without generating an error"""
    i = 0    
    while True:
        try:
            return (i,)+FindDigPeaks(data[i], N, sep)
        except RuntimeError:
            pass
        i += 1
        if(i == data.shape[0]):
            raise RuntimeError('No suitable line found')
            
 

def MaskByLines(data,line1, line2):
    """mask data so it is zero outside line1 and line2. line1 and line2 are
    tuples to go into Poly() line1 is lower than line2"""
    out = np.array(data)
    x = np.arange(data.shape[0])
    y1 = Poly(x, *line1)
    y2 = Poly(x, *line2)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if(not((j > y1[i]) & (j < y2[i]))):
                out[i,j] = 0.
    return out
    


def GetPeak(data, y1, y2, x):
    """"""
    return np.array(data[x,y1:y2])  
            


def GuessPeakNum(data, smoothwidth = 10, atrow = 300, thresh = 0.05):
    """Guess the number of slitlines in data """
    #first smooth the data:
    con_x = np.linspace(-50,50,101)
    con_y = np.exp(-con_x*con_x/(0.5*smoothwidth*smoothwidth))
    con_y /= np.sum(con_y)
    
    peaks = np.sum(data[atrow-20:atrow+20], axis = 0)
    
    smooth_peaks = np.convolve(peaks, con_y, mode = 'same')
    
    max_x = argrelextrema(smooth_peaks, np.greater)[0]
    max_y = smooth_peaks[max_x]
    max_y /= max_y.max()
    return np.argwhere(max_y>thresh).size
    



def GuessSinglePeak(data):
    """Make intital guesses for a single gaussian in data"""
    height = data.max()
    
    widthrange = np.argwhere(data > (0.5*height)).transpose()[0]
    if((len(widthrange) < 2)):
        raise RuntimeError('No peak found')
    elif(np.gradient(widthrange).max()>2): ##separated peaks means gradient higher than 2
        raise RuntimeError('Too many peaks found')
    else:
        width = (widthrange[-1] - widthrange[0] + 1)
        pos = (widthrange[-1] + widthrange[0])/2.
        return (height, pos, width)



def FindSlitCurve(data, y1, y2, n = 5, fit = True, relcutoff = 0.02):
    """track the curve in data that is between y1 and y2 (np arrays), just find
    the max of the curve if fit = False, otherwise fit the peak with a Gaussian
    return a np array with x values and a np array with y values"""
    cutoff = relcutoff*data.max()
    i = 0
    xvalues = []
    peaks = []
    while True:
        if(i >= data.shape[0]):
            break
        try:    
            yval, xval = FindPeakInCluster(data[i:i+n],y1[i:i+n], y2[i:i+n], fit = fit, cutoff = cutoff)
            peaks.append(yval + y1[i + xval])
            xvalues.append(i + xval)
        except RuntimeError:
            pass
        i += n
    return np.array(peaks), np.array(xvalues)



def FindPeakInCluster(data, y1, y2, fit = True, cutoff = -1.): #-1 means no cutoff
    """Find the peak position in a cluster, if fit is true with fitting, if false, no fit"""
    if(cutoff == -1): #cutoff
        cutoff = data.max()
    if fit:
        i = 0
        while True:
            if (i >= data.shape[0]):
                raise RuntimeError('No fittable peak found')
            try:
                guess = GuessSinglePeak(data[i,int(y1[i]):int(y2[i])])
                avg = np.average(data[i])
                if((len(guess) == 3) & ((guess[0]/avg) > 3.) & (guess[0] > cutoff)): ## check if guesses are suitable (there can only be one peak)
                    x = np.arange(int(y2[i]) - int(y1[i]))
                    try:
                        fitparams, covariance = curve_fit(MPGauss, x, data[i,int(y1[i]):int(y2[i])], p0=guess, maxfev = 100)
                        if((fitparams[1] > x[0]) & (fitparams[1] < x[-1])): #check position is between boundaries and weight is at least 3 times the BG
                            return fitparams[1], i
                    except RuntimeError:
                        pass
            except RuntimeError:
                pass
            i += 1
    else:
        return np.argmax(data[int(data.shape[0]/2.)]), int(data.shape[0]/2.)
        
        


       
def Gauss(x, A, x0, sig):
    return A * np.exp(-(x-x0)**2 /sig)
    


def MPGauss(x, *params):
    """Multi peak gaussian for fitting purposes, params/3 is the number of peaks
    params[i+0] = A_i, params[i+1] = x0_i, params[i+2] = sig_i"""
    n_peaks = int(len(params)/3)
    value = np.zeros(len(x)) #copy array
    
    for i in range(n_peaks):
        value += Gauss(x, *(params[3*i:3*i+3]))
        
    return value  

def Poly(x, *params):
    """Polynomial for fitting. returns params[0] + params[1]*x + params[2]*x^2 etc"""
    N = len(params)
    
    y = np.zeros(x.shape)    
    
    for i in range(N):
        y += params[i]*(x**i)
    
    return y
    
    
    
def NPeakGuess(data, N = 0, step = 0.01, ratio = 0.1, absmin = 0.): # check if 0.01 is a consistent value
    """ Find N peaks in data and return their parameters in a list of
    step is the ratio between the maximum of data and the sampling interval.
    For low peaks, lower values might be needed. Data is a np.array
    If N = 0, keep finding peaks until ratio is reached. absmin is the absolute
    minimum value after which the algorithm stops looking
    """

    risexlist = []
    fallxlist = []
    level = data.max()
    peakheights = [] 
    peakpos = []
    width = []

    
    while True:
        risexlist = []
        fallxlist = []
    
        
        for ii in range(len(data)-1):
            if((data[ii]<level)&(data[ii+1]>level)):
                risexlist.append(ii)
            if((data[ii]>level)&(data[ii+1]<level)):
                fallxlist.append(ii)    
                
        #append height for number of found peaks        
        peakheights += (len(risexlist)-len(peakheights))*[level]

        
        #append position of new peaks
        if len(peakpos)<len(risexlist): # if new position found
            for pos in peakpos: #loop over old positions
                for (risex,fallx) in zip(risexlist,fallxlist): #and check where they are in the list
                    if (pos > risex) & (pos < fallx): # if found: remove them
                        risexlist.remove(risex)
                        fallxlist.remove(fallx)
            
            for (risex,fallx) in zip(risexlist,fallxlist): #now add the remaining positions
                peakpos.append( ( risex + fallx ) / 2. )
    
    
        level -= step*data.max()
            
        if(level < absmin): ## stop looking if level is < absmin and error
            return ()
        
        if (N != 0) & (len(peakheights) > N):## return if too many peaks are found
            return ()
        
        if ( (N == 0) & ( level < ratio*data.max() )): #stop looking if level < ratio level
            break # only if unknown number of peaks
 
        if (N != 0) & (len(peakheights) == N):## break if all peaks are found
            break
  

    #find peakwidths
    for (pos,height) in zip(peakpos,peakheights):
        #right side
        ii = round(pos)
        while (data[ii] > height / 2. ):
            ii += 1
        jj = round(pos)
        while (data[jj] > height / 2. ):
            jj -= 1
        if((ii-jj-1) > 0): # check if it isn't a lonely spike
            width.append(ii-jj-1)
        else:
            return () # if it is, it means the data is too noisy to find the right amount of peaks
    #return params as a tuple
    peakparams = tuple(sum([[peakheights[i],peakpos[i],width[i]] for i in range(len(peakheights))],[]))    
    
    return peakparams





def SlitCorrect(data, t):
    """" t is transformation matrix"""
    out = np.array(t)

    x = np.arange(data.shape[0])
    y = np.arange(data.shape[1])
    interp = interp2d(y, x, data, bounds_error=False)    
    
    for i in range(out.shape[0]):
        out[i,:] = interp(t[i],i) #  actually interp returns an array with x_new and y_new as input arays interp(array_x, array_y)
            
    return out





#maybe divide this into smaller subfunctions
#MultiGaussFit should only contain the guess function and the fit function
# for one line, then implement this into a determine slits function
#This should not loop over data, but just take one line and guess and fit it.
def MultiGaussFit(data, N = 17, step = 0.02, ratio = 0.2, absmin = 0.):
    """ Fit N Gaussians to data, or as much as peakfind can find if N = 0
    
    """

    #guess
    guesses = NPeakGuess(data, N = N, step = step, ratio = ratio, absmin = absmin)

    if(guesses == ()):# if something went wrong
        #raise error
        raise RuntimeError('No Guesses')
        return [] #return empty list
        

    x = np.arange(1024)   

    #fit, include error handling
    try:    
        fitparams, covariance = curve_fit(MPGauss, x, data, p0=guesses, maxfev = 1000) 
    except RuntimeError as error:
        print(error)
        raise RuntimeError('Fit not converging')
    return fitparams
 

########
#to do: make a transformation array: on x-axis energy pix, y-axis is angle, 
#the value of each pixel should be the pixel y-number where that angle is in the slit array image
#The array can be found by interpolating between the fitcurves
# better would be to fit the fitcurves as a function of E-pix and then calculate instead of interpolating
# after this we get a transformation array, which is the size of the new image, which contains the
# angle-pixle values from the raw data
#####   

######
#for correcting the data later: you can make a similar array like that from the gold
#correction data, this time for the x-pixel values (energy)
# the interpolation can be done as a 2D interpolation at the end with the (x,y) values
# from the two arrays


####
#after that write 2D data and 3D data objects and make slit and gold objects
#slitcorrection and goldcorrection should also be objects
#then two objects: rawdata (inherits from 2Ddata)
# and angledata (inherits from 2Ddata)
# then you can transform rawdata into an angledata object by correcting with slits and gold
# Then make an export to hdf5 function that writes the same hdf5 files as Giorgio
# so I can actually use this code for corrections and compare
# Then I can also k-correct the data in IDL
def AnalyzeSlitsOld(data, N = 17): ## put everything together
    """Make a slit transformation matrix from DATA, with N slit lines"""        
    xvalues, peaks = FitSlitLines(data, N = N)
    SortPos2(peaks)
    curveparams = FitSlitCurves(peaks,xvalues)
    T,x,y = MakeTransformationMatrix(curveparams, data.shape)
    return T



def AnalyzeSlitsNoFit(data, N = 17):
    """ analyze slits without fitting"""
    xvalues, peaks = FindSlitLines(data, N = N)
    SortPos2(peaks)
    curveparams = FitSlitCurves(peaks, xvalues)
    T,x,y = MakeTransformationMatrix(curveparams, data.shape)
    return T


# make such that this one takes curveparams instead of fitcurves
def MakeTransformationMatrix(curveparams, datashape): ## **give this a better name too
    """Define a 2D image that has pixel values for a matrix of angles.
    CURVEPARAMS contains the results of the fits, the middle slit is at CENTER"""

    dd = 1. #** distance per slit
    dl = 33.7  #should be 1 mm per slit at a distance of 33.7 mm
    
    n_slits = curveparams.shape[0] # number of slits
    dstart = -((n_slits-1)*dd/2.) # angle of first slit
    dend = dstart + dd*(n_slits-1)
    #make xscaling for "fitcurves": (angle dimension)
    slit_x_scale = np.tan(np.linspace(dstart, dend, n_slits)/dl)*180/np.pi

    

    ysize = 1024 # checked from some gold file
    x = np.arange(datashape[0])
    y = np.linspace(-20, 20, ysize) #these angles are emperical and hard coded for the UBC 245 system
    A = np.zeros((datashape[0],ysize))
    
    #take params and convert back to lines using the x-scaling
    fitcurves = np.zeros((n_slits,x.shape[0]))
    for i in range(n_slits):
        fitcurves[i,:] = Poly(x, *curveparams[i])
 

    

    for i in range(len(x)):
        #interp1d returns a function that interpolates a whole array at once
        interp = interp1d(slit_x_scale,fitcurves[:,i],fill_value='extrapolate') 
        A[i,:] = interp(y)
    
    return A, x, y




def FitSlitCurves(curves, x, N_p = 2):
    """Fit CURVES (n*curve) of the slits and return the parameters in an array
    with X as xdata."""
    N = curves.shape[0]# number of slit curves
    #polynomial fit order. 2 seems good
    result = np.zeros((N, N_p+1))
    for i in range(N):
        result[i] = FitSlitCurve(curves[i], N_p, x)
    
    return result

def FitSlitCurve(curve, N_p, x):
    """Fit the CURVE with a N order polynomial with X as xdata, return fit parameters"""
    guesses = (0,)*(N_p+1)
    try:    
        fitparams, covariance = curve_fit(Poly, x, curve, p0=guesses, maxfev = 1000) 
    except RuntimeError as error:
        print(error)
        raise RuntimeError('Fit not converging')    
    return fitparams


def SortPos(temp, xvalues, max_d = 0.3):
    
    ii = 0
    peaks = np.array(temp)
    while True: #loop until all curves are sorted
        lastval = peaks[ii,0] # init vals
        lastxval = xvalues[0]
        maxskips = 20
        skips = 0
        start = 1 #value is 1 when starting, set to 0 when first value is met
        for jj in range(peaks.shape[1]-1):
            #when starting:
            if((start == 1) & (peaks[ii,jj] == -1)): # do nothing, but swap to put -1 at end for jj+1
                if(jj == 0):
                    for kk in range(peaks.shape[1]):
                        swap1s(peaks[ii:,kk])#put -1's of jj+1 at the end
                lastval = peaks[ii,jj+1] # already set lastval to the right value for next time
                lastxval = xvalues[jj+1] # and lastxval
            else:#this means we aren't at the start
                start = 0
                swap(peaks[ii:,jj+1],lastval) # try to find something similar to the last value
                d = abs((peaks[ii,jj+1]-lastval)/(xvalues[jj+1]-lastxval)) # check if it belongs to the same curve
                if(d > max_d):
                    if(not (peaks[-1,jj+1] == -1)): #check if there is an empty element at the end of the series
                        peaks = pad(peaks)
                    peaks[-1,jj+1] = peaks[ii,jj+1] #swap with -1
                    peaks[ii,jj+1] = -1
                    skips += 1 # increaese skips

                    if(skips > maxskips): #if more then 3 skips, replace the rest by -1's
                        for kk in range(jj+2, peaks.shape[1]):
                            peaks[-1,kk] = peaks[ii,kk] 
                            peaks[ii,kk] = -1                
                        break # and exit the for loop

                if(peaks[ii,jj] == -1): #if we still find a -1 it means it was created just before
                    pass # I don't think anything needs to be done here now.
                else:
                    lastval = peaks[ii,jj+1]
                    lastxval = xvalues[jj+1] #update x and y values to compare to
                    skips = 0 # reset skips

        ii += 1    # to get out of the while loop
        #print(ii)
        if(ii>=peaks.shape[0]):
            break
    #count -1's
    curvelength = np.zeros((peaks.shape[0]))
    curves = []
    for ii in range(peaks.shape[0]):
        for jj in range(peaks.shape[1]):
            if(peaks[ii,jj] == -1):
                curvelength[ii] += 1
        if (curvelength[ii] < 100):
            curves.append(ii)
            
    
    
    return peaks, curves




def SortPos1(peaks, xvalues, max_d = 0.3):
    """Sort positions from fitting into lists that belong togehter so they can be fitted"""
    ii = 0
    # try sort-split-combine
    #sort: put elements in increasing order
    # split: split wherever d > max_d
    # combine: try to combine the shorter separate curves

    while True: #loop until all curves are sorted
        lastval = peaks[ii,0]
        lastxval = xvalues[0]
        maxskips = 3
        skips = 0
        ######
        ## if the -1's are at the start, skip over them until the first value is found 
        #######        
        for jj in range(peaks.shape[1]-1):
            # to check if it belongs to the same curve
            if(peaks[ii,jj] == -1): # if get to an "undefined" value at start
                swap1s(peaks[ii:,jj+1])#put -1's of jj+1 at the end
                if(ii == 17):
                     #print(peaks[ii,jj+1])
                     print(abs((peaks[ii,jj+1]-lastval)/(xvalues[jj+1]-lastxval)))
            else:   #if the next one is too big or too small
                #first put peak at jj+1,ii that is closest to jj,ii
                swap(peaks[ii:,jj+1],peaks[ii,jj])
                #######
                #compare peaks to the last known -1 value
                #######
                #######
                #only update this value (lastval) if it wasn't -1
                #######
                if(ii == 17):
                    print(peaks[ii,jj])
            ## add if too many -1's, should switch to next curve
            d = abs((peaks[ii,jj+1]-lastval)/(xvalues[jj+1]-lastxval))
            if(d > max_d): #if this isn't the same curve
                if(not (peaks[-1,jj+1] == -1)): #check if there is an empty element at the end of the series
                    #if there isn't, make space:
                    temp = np.zeros((peaks.shape[0]+1, peaks.shape[1]))
                    temp[:-1,:] = peaks
                    temp[-1,:] = -1
                    peaks = np.array(temp)
                                
                peaks[-1,jj+1] = peaks[ii,jj+1] #and swap
                peaks[ii,jj+1] = -1
                skips += 1 # increaese skips
                    
                if(skips > maxskips): #if more then 3 skips, replace the rest by -1's
                    for kk in range(jj+2, peaks.shape[1]):
                        peaks[-1,kk] = peaks[ii,kk] 
                        peaks[ii,kk] = -1                
                    break # and exit the for loop
                    
                
            else: # if it is the same curve:
                lastval = peaks[ii,jj+1]
                lastxval = xvalues[jj+1] #update x and y values to compare to
                skips = 0 # reset skips

#                if(ii == peaks.shape[1]-1):
#                   pass #
#                else: #swap it with the right in peaks[jj,ii+1:]
#                   pass #make func that swaps
        ii += 1    
        #print(ii)
        if((ii == 25) | (ii>peaks.shape[0])):
            break
        
    return peaks
                
                    
def SortPos2(peaks):
    for ii in range(peaks.shape[1]):
        peaks[:,ii].sort()
    
    

def score(array1, array2):
    out = []
    for ii in range(len(array1),0):
        out.append(0)
        for jj in range(0,len(array1)-ii):
            print('ii = ', ii)
            print('jj = ', jj)
            out[ii] = abs(array1[jj]-array2[ii+jj])    
    for ii in range(1, len(array1)):
        out.append(0)
        for jj in range(0,len(array1)-ii):
            print('ii = ', ii)
            print('jj = ', jj)
            out[ii] = abs(array1[ii+jj]-array2[jj])
    return out
        

def padlow(array):
    """put a row of zeros at the start of the array"""
    temp = np.zeros((array.shape[0]+1, array.shape[1]))
    temp[1:,:] = array
    temp[0,:] = -1
    return temp

def padhigh(array):
    """put a row of zeros at the end of the array"""
    temp = np.zeros((array.shape[0]+1, array.shape[1]))
    temp[:-1,:] = array
    temp[-1,:] = -1
    return temp

def pad(array):
    """put a row of zeros at the end of the array"""
    temp = np.zeros((array.shape[0]+1, array.shape[1]))
    temp[:-1,:] = array
    temp[-1,:] = -1
    return temp


def swap(array, val): 
    """return array like array where array[0] contains the value closest to val"""
    ii_close = 0
    for ii in range(array.shape[0]):
        if (abs(array[ii]-val) < abs(array[ii_close]-val)):
            ii_close = ii
    #swap:
    temp = array[0]
    array[0] = array[ii_close]
    array[ii_close] = temp
    
    
def swap1s(array):
    """put -1's at the end of the array"""
    ii = 0
    jj = 0
    while(jj < array.shape[0]):
        if (array[ii] == -1):
            array[ii:-1] = array[ii+1:]
            array[-1] = -1
        else:
            ii += 1
        jj += 1


def FitSlitLines(data, n = 5, N = 17, step = 0.02, ratio = 0.2, absminratio = 0.01):
    """loop over 2D data, and fit all lines with multipeaks """
    #call GetLinePositions for a set of lines. Once every 5 lines seems OK
    absmin = absminratio*data.max()
    ii = 0
    peaks = [] # list with peak positions in the line
    xvalues = [] #list with array positions of where the line was taken
    while True:
        try:
            xvalue,peak = GetPosFromCluster(data[n*ii:n*(ii+1),:], N = N, step = step, ratio = ratio, absmin = absmin)
            peaks.append(peak)
            xvalues.append(ii*n + xvalue) # append position in the array
        except RuntimeError:
            pass
        ii += 1
        print(n*ii)
        if(ii*n > data.shape[0]):
            break
    return np.array(xvalues),np.array(peaks).transpose()
 

def FindSlitLines(data, n = 5, N = 17, step = 0.02, ratio = 0.2, absminratio = 0.01):
    """loop over 2D data, and fit all lines with multipeaks """
    #call GetLinePositions for a set of lines. Once every 5 lines seems OK
    absmin = absminratio*data.max()
    ii = 0
    peaks = [] # list with peak positions in the line
    xvalues = [] #list with array positions of where the line was taken
    while True:
        try:
            xvalue,peak = GuessPosFromCluster(data[n*ii:n*(ii+1),:], N = N, step = step, ratio = ratio, absmin = absmin)
            peaks.append(peak)
            xvalues.append(ii*n + xvalue) # append position in the array
        except RuntimeError:
            pass
        ii += 1
        print(n*ii)
        if(ii*n > data.shape[0]):
            break
    return np.array(xvalues),np.array(peaks).transpose()
 

def GuessPosFromCluster(data, N = 17, step = 0.02, ratio = 0.2, absmin = 0.):
    """ return the positions of the first line that fits without generating an error"""
    ii = 0    
    while True:
        try:
            return ii,GuessLinePositions(data[ii], N = N, step = step, ratio = ratio, absmin = absmin)
        except RuntimeError:
            pass
        ii += 1
        if(ii == data.shape[0]):
            raise RuntimeError('No Fittable line found')
     
def GetPosFromCluster(data, N = 17, step = 0.02, ratio = 0.2, absmin = 0.):
    """ return the positions of the first line that fits without generating an error"""
    ii = 0    
    while True:
        try:
            return ii,GetLinePositions(data[ii], N = N, step = step, ratio = ratio, absmin = absmin)
        except RuntimeError:
            pass
        ii += 1
        if(ii == data.shape[0]):
            raise RuntimeError('No Fittable line found')
            
            

def GetLinePositions(data, N = 17, step = 0.02, ratio = 0.2, absmin = 0.):
    """Get peak positions from a line"""
    fitresults = MultiGaussFit(data, N = N, step = step, ratio = ratio, absmin = absmin)    
    # extract peak positions and return
    positions = [fitresults[3*ii+1] for ii in range(round(len(fitresults)/3))]
    return positions
    
def GuessLinePositions(data, N = 17, step = 0.02, ratio = 0.2, absmin = 0.):
    """Get peak positions from a line"""
    guesses = NPeakGuess(data, N = N, step = step, ratio = ratio, absmin = absmin)
    if(guesses == ()):# if something went wrong
        #raise error
        raise RuntimeError('No Guesses')
    positions = [guesses[3*ii+1] for ii in range(round(len(guesses)/3))]
    return positions
    

    
    
    
    