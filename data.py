# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 20:12:25 2015

@author: berend
"""


import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
# from arpes.slits import AnalyzeSlits, SlitCorrect
from slits_new import AnalyzeSlits, SlitCorrect
# from arpes.gold import AnalyzeGold, AnalyzeGold1D, GoldCorrect, Gold1DCorrect, FDexpfit, FDexp
from gold_new import AnalyzeGold, AnalyzeGold1D, GoldCorrect, Gold1DCorrect, FDexpfit, FDexp
# from arpes.kcorrection import anglesfromk,kfromangles,anglesfromkz
from kcorrection_new_new import anglesfromk,kfromangles,anglesfromkz
import h5py
from scipy.interpolate import RectBivariateSpline,interpn
# import arpes.data_loaders as dl
import data_loaders as dl
# import arpes.plotters as pl
import plotters as pl
from matplotlib.path import Path
from warnings import warn



class Data1D():
    def __init__(self, shape = (1,), name = 'NoName', note = ''):
        self.data = np.zeros(shape, dtype = 'float32')
        self.datalabel = ''
        self.shape = shape
        self.xaxis = np.arange(shape[0], dtype = 'float32')
        self.xstart = 0
        self.xdelta = 1        
        self.xlabel = ''
        self.name = name
        self.note = note


    def reshape(self, newshape):
        tempdata = np.zeros(newshape, dtype = 'float32')
        xlen = min(newshape[0],self.shape[0])
        tempdata[0:xlen] = self.data[0:xlen]
        self.shape = newshape
        self.data = tempdata
        self.setXScale(self.xstart, self.xdelta)
    
    def setXScale(self,xstart, xdelta):
        self.xstart = xstart
        self.xdelta = xdelta
        self.xaxis = (np.arange(self.data.shape[0], dtype = 'float32')*self.xdelta)+self.xstart
              
    def plot(self, **kwargs):
        fig, ax = plt.subplots(figsize=(6,4))
        ax.plot(self.xaxis,self.data)
        ax.set_ylabel = self.datalabel
        ax.set_xlabel = self.xlabel



class Data2D():
    """A class containing 2D data with x and y axes scaling"""
    def __init__(self, shape = (1,1), name = "NoName", note = ""):
        self.data = np.zeros(shape, dtype = 'float32')
        self.datalabel = ''
        self.shape = shape
        self.xaxis = np.arange(shape[0], dtype = 'float32')
        self.yaxis = np.arange(shape[1], dtype = 'float32')
        self.xstart = 0
        self.xdelta = 1        
        self.ystart = 0
        self.ydelta = 1
        self.xlabel = ''
        self.ylabel = ''
        self.name = name
        self.note = note


    def reshape(self, newshape):
        tempdata = np.zeros(newshape, dtype = 'float32')
        xlen = min(newshape[0],self.shape[0])
        ylen = min(newshape[1],self.shape[1])
        tempdata[0:xlen,0:ylen] = self.data[0:xlen,:ylen]
        self.shape = newshape
        self.data = tempdata
        self.setXScale(self.xstart, self.xdelta)
        self.setYScale(self.ystart, self.ydelta)    
    
    def setXScale(self,xstart, xdelta):
        self.xstart = xstart
        self.xdelta = xdelta
        self.xaxis = (np.arange(self.data.shape[0], dtype = 'float32')*self.xdelta)+self.xstart
        
    def setYScale(self,ystart, ydelta):
        self.ystart = ystart
        self.ydelta = ydelta        
        self.yaxis = (np.arange(self.data.shape[1], dtype = 'float32')*self.ydelta)+self.ystart

                                



class Data3D():
    """A class containing 3D data with x, y and z axes scaling"""
    def __init__(self, shape = (1,1,1), name = "NoName", note = ""):
        self.data = np.zeros(shape, dtype = 'float32')
        self.datalabel = ''
        self.shape = shape
        self.xaxis = np.arange(shape[0], dtype = 'float32')
        self.yaxis = np.arange(shape[1], dtype = 'float32')
        self.zaxis = np.arange(shape[2], dtype = 'float32')
        self.xstart = 0
        self.xdelta = 1        
        self.ystart = 0
        self.ydelta = 1
        self.zstart = 0
        self.zdelta = 1
        self.xlabel = ''
        self.ylabel = ''
        self.zlabel = ''
        self.name = name
        self.note = note


    def reshape(self, newshape):
        tempdata = np.zeros(newshape, dtype = 'float32')
        #fill in tempdata with self.data
        xlen = min(newshape[0],self.shape[0])
        ylen = min(newshape[1],self.shape[1])
        zlen = min(newshape[2],self.shape[2])
        tempdata[0:xlen,0:ylen,0:zlen] = self.data[0:xlen,:ylen,0:zlen]
        self.data = tempdata
        self.shape = newshape
        self.setXScale(self.xstart, self.xdelta)
        self.setYScale(self.ystart, self.ydelta)
        self.setZScale(self.zstart, self.zdelta)
    
    
    def setXScale(self,xstart, xdelta):
        self.xstart = xstart
        self.xdelta = xdelta
        self.xaxis = (np.arange(self.data.shape[0], dtype = 'float32')*self.xdelta)+self.xstart
        
    def setYScale(self,ystart, ydelta):
        self.ystart = ystart
        self.ydelta = ydelta        
        self.yaxis = (np.arange(self.data.shape[1], dtype = 'float32')*self.ydelta)+self.ystart
    
    def setZScale(self,zstart, zdelta):
        self.zstart = zstart
        self.zdelta = zdelta        
        self.zaxis = (np.arange(self.data.shape[2], dtype = 'float32')*self.zdelta)+self.zstart
        
     

#class for 1D ARPES data (EDCS)
class AaData1D(Data1D):
    def __init__(self, path, objtype = 'AaData1D', datatype = None, name = '', note = '', pin = dl.std_p()):
        super(AaData1D, self).__init__(name = name, note = note)
        
        if datatype is None:
            datatype = 'load'


            
        if datatype == 'load':
            self.p = dl.loadH5(self, path, objtype = objtype) #loading already saved data
        elif datatype == 'empty':
            self.p = {}
        elif datatype == 'APE':
            self.p = dl.loadData1DAPE(self,path = path)
        else:
            raise RuntimeError('invalid datatype')




    def crop(self, xmin = None, xmax = None):
        if xmin is None:
            xmin = self.xstart
        if xmax is None:
            xmax = self.xstart + self.xaxis.size * self.xdelta

            
        #get min_n:
        xmin_n = int((xmin-self.xstart)/self.xdelta)
        xmax_n = int((xmax-self.xstart)/self.xdelta)

        
        return self.cropN(xmin_n,xmax_n)
    
    
    def cropN(self, xmin_n = None, xmax_n = None):
        if xmin_n is None:
            xmin_n = 0
        if xmax_n is None:
            xmax_n = self.xaxis.size
  
        xmin_n,xmax_n = (min(xmin_n,xmax_n),max(xmin_n,xmax_n))
            
        newobj = AaData1D('', datatype = 'empty', shape = (xmax_n-xmin_n), name = self.name, note = self.note)
        newobj.data = self.data[xmin_n:xmax_n]
        newobj.p = self.p.copy()

        newobj.setXScale(self.xstart + self.xdelta*xmin_n, self.xdelta)        
        return newobj    

        
    def rebin(self, binning = (1,)):
        
        shape = self.data.shape  
        #reduce to nearest integer multiple
        redshape = (int(shape[0]//binning[0])*binning[0],)

        self.reshape(redshape)

        newshape = (int(redshape[0]//binning[0]),)
        compression_pairs = [(d, int(c//d)) for d,c in zip(newshape,redshape)]


        flattened = [l for p in compression_pairs for l in p]

        newdata = self.data.reshape(flattened)
        for i in range(len(newshape)):
            newdata = newdata.mean(-1*(i+1))
        
        newobj = AaData1D('', datatype = 'empty', name = self.name, note = self.note)
        newobj.reshape(newshape)
        newobj.p = self.p.copy()
        newobj.setXScale(self.xstart,binning[0]*self.xdelta)
        
        newobj.data = newdata
        
        return newobj
                                
                
    def setangles(self, theta0 = None, phi0 = None, alpha0 = None, theta_m = None, phi_m = None):
        if not (theta0 == None):
            self.p['theta0'] = theta0
        if not (phi0 == None):
            self.p['phi0'] = phi0
        if not (alpha0 == None):
            self.p['alpha0'] = alpha0
        if not (theta_m == None):
            self.p['theta_m'] = theta_m
        if not (phi_m == None):
            self.p['phi_m'] = phi_m
        print('Set angles:')
        self.printangles()
        
    def getangles(self):
        return (self.p['theta0'],self.p['phi0'],self.p['alpha0'], self.p['theta_m'],self.p['phi_m'])
        
    def printangles(self):
        print('theta offset: %s' % self.p['theta0'])
        print('phi offset: %s' % self.p['phi0'])
        print('rotation: %s' % self.p['alpha0'])
        print('theta manipulator: %s' % self.p['theta_m'])
        print('phi manipulator: %s' % self.p['phi_m'])        
        
            
    # from arpes.plotters import plot1D
    from plotters import plot1D

    def show(self, mode='row', **kwargs):
        """General plotfunction
        Args:
            mode: string indicating which plotfunction to run. Run as 'help' for available modes
            kwargs: keyword args necessary to run the chosen plot function
        
        """
        
        warn("show is deprecated, use the direct functions prepended with plot")


        if mode is None:
            if(self.p['type'] == 'rawdata') or (self.p['type'] == 'angledata'):
                mode = 'col'
            else:
                mode = 'k'


        plotfuncdict = {
                        'row' : self.plot1D,
                        'plot' : self.plot1D,
                        'EDC' : self.plot1D,
                        }
                        
        if mode == 'help':
            for k,v in plotfuncdict.items():
                print('Available mode \'{}\', can be called directly as obj.{}, callsign and docstring: '.format(k,v.__name__))
                help(v)
                print('\n\n\n')
            return
        
                        
        if(not mode in plotfuncdict):
            raise ValueError('Unkown plotfunction mode')
            

        ax = plotfuncdict[mode](**kwargs) #run the plot function

        return ax


    
    def QuickCorrect(self, EF):
        """Quick method of correcting so k-plotting is available, only use offset, no calibration
        Args:
            EF: energy to which to set the Fermi level"""

        self.setXScale(self.xaxis[0] - EF, self.xaxis[1] - self.xaxis[0])
        self.p['type'] = 'corrected'
        self.p['EkatEF'] = EF
            
        


    def Correct(self, goldobj = None):

            
        if goldobj:
            #gold correction        
            if(not goldobj.p['analyzed']):
                raise RuntimeError('Analyze gold first')
            goldcorrected, newxaxis = Gold1DCorrect(self.data, self.xaxis, goldobj.p['EkatEF'])
            self.setXScale(newxaxis[0], newxaxis[1]-newxaxis[0])
    
            #copy and finish up
            np.copyto(self.data,goldcorrected)   
            self.p['type'] = 'corrected'
            self.p['EkatEF'] = goldobj.p['EkatEF']

    def CastGold1D(self):
        pass
                    

            
    def save(self, path = '', objtype = 'AaData1D'):
        if path == '':
            if not 'savepath' in self.p:
                raise RuntimeError('No previously known save path')
            else:
                path = self.p['savepath']
        else:
            if((path.find('.h5') == -1) and (path.find('.hdf5') == -1)):
                path += '.h5'
            self.p['savepath'] = path

                
        with h5py.File(path,'w') as hf:
            #header:
            hf.clear()
            hf.create_dataset('h/ver', data =  np.array('1.0'.encode('utf-8')))
            hf.create_dataset('h/obj', data =  np.array(objtype.encode('utf-8')))
            #data:    
            hf.create_dataset('d/data', data = self.data)
            hf.create_dataset('d/datalabel', data = np.array(self.datalabel.encode('utf-8')))
            hf.create_dataset('d/xstart', data = np.array(self.xstart))
            hf.create_dataset('d/xdelta', data = np.array(self.xdelta))
            hf.create_dataset('d/xlabel', data = np.array(self.xlabel.encode('utf-8')))
            hf.create_dataset('d/name', data = np.array(self.name.encode('utf-8')))
            hf.create_dataset('d/note', data = np.array(self.note.encode('utf-8')))                
            #parameters:            
            for k,v in self.p.items():
                if (type(v) == str) or (type(v) == np.str_):
                    hf.create_dataset('p/str/' + k, data = np.array(v.encode('utf-8')))
                elif v is None:                   
                    hf.create_dataset('p/none/' + k, data =  np.array(''.encode('utf-8'))) #save None's as an empty string
                elif type(v) is bool:                   
                    hf.create_dataset('p/bool/' + k, data =  np.array([v])) #save bools separately
                else:
                    hf.create_dataset('p/num/' + k, data =  np.array([v]))
            
            

                              

        
        

#class for loading 2D data
# has attributes like history, properties, measurement parameters
# loads binary, loads raw (?)
class AaData2D(Data2D):
    def __init__(self, path, objtype = 'AaData2D', datatype = None, name = '', note = '', mode = 'avg', pin = dl.std_p()):
        super(AaData2D, self).__init__(name = name, note = note)
        
        if datatype is None:
            datatype = 'load'

            
        if datatype == 'load':
            self.p = dl.loadH5_2D(self,path, objtype = objtype) #loading already saved data
        elif datatype == 'empty':
            self.p = {}
        elif datatype == 'UBC':
            self.p = dl.loadData2DUBC(self,path = path, mode = mode, pin = pin)            
        elif datatype == 'SLS':
            self.p = dl.loadData2DSLS(self,path = path, pin = pin)
        elif datatype == 'igor':
            self.p = dl.loadData2DIgor(self,path = path, mode = mode, pin = pin)
        elif datatype == 'ALS':
            self.p = dl.loadData2DALS(self,path = path, mode = mode, pin = pin)
        elif datatype == 'APE':
            self.p = dl.loadData2DAPE(self,path = path, pin = pin)
        elif datatype == 'CLS':
            self.p = dl.loadData2DCLS(self,path = path, pin = pin)
        else:
            raise RuntimeError('invalid datatype')


    def crop(self, xmin = None, xmax = None, ymin = None, ymax = None):
        if xmin is None:
            xmin = self.xstart
        if xmax is None:
            xmax = self.xstart + self.xaxis.size * self.xdelta
        if ymin is None:
            ymin = self.ystart
        if ymax is None:
            ymax = self.ystart + self.yaxis.size * self.ydelta

            
        #get min_n:
        xmin_n = int((xmin-self.xstart)/self.xdelta)
        xmax_n = int((xmax-self.xstart)/self.xdelta)
        ymin_n = int((ymin-self.ystart)/self.ydelta)
        ymax_n = int((ymax-self.ystart)/self.ydelta)

        
        return self.cropN(xmin_n,xmax_n,ymin_n,ymax_n)
    
    
    def cropN(self, xmin_n = None, xmax_n = None, ymin_n = None, ymax_n = None):
        if xmin_n is None:
            xmin_n = 0
        if xmax_n is None:
            xmax_n = self.xaxis.size
        if ymin_n is None:
            ymin_n = 0
        if ymax_n is None:
            ymax_n = self.yaxis.size
  
        xmin_n,xmax_n = (min(xmin_n,xmax_n),max(xmin_n,xmax_n))
        ymin_n,ymax_n = (min(ymin_n,ymax_n),max(ymin_n,ymax_n))
            
        newobj = AaData2D('', datatype = 'empty', name = self.name, note = self.note)
        newobj.reshape((xmax_n-xmin_n, ymax_n-ymin_n))
        newobj.data = np.array(self.data[xmin_n:xmax_n,ymin_n:ymax_n])
        newobj.p = self.p.copy()

        newobj.setXScale(self.xstart + self.xdelta*xmin_n, self.xdelta)
        newobj.setYScale(self.ystart + self.ydelta*ymin_n, self.ydelta)
        
        return newobj    

        
    def rebin(self, binning = (1,1)):
        
        shape = self.data.shape  
        #reduce to nearest integer multiple
        redshape = (int(shape[0]//binning[0])*binning[0],int(shape[1]//binning[1])*binning[1])

        self.reshape(redshape)

        newshape = (int(redshape[0]//binning[0]),int(redshape[1]//binning[1]))
        compression_pairs = [(d, int(c//d)) for d,c in zip(newshape,redshape)]


        flattened = [l for p in compression_pairs for l in p]

        newdata = self.data.reshape(flattened)
        for i in range(len(newshape)):
            newdata = newdata.mean(-1*(i+1))
        
        newobj = AaData2D('', datatype = 'empty', shape = newshape, name = self.name, note = self.note)
        newobj.p = self.p.copy()
        newobj.setXScale(self.xstart,binning[0]*self.xdelta)
        newobj.setYScale(self.ystart,binning[1]*self.ydelta)
        
        newobj.data = newdata
        
        return newobj
                
                
                
    def setangles(self, theta0 = None, phi0 = None, alpha0 = None, theta_m = None, phi_m = None):
        if not (theta0 == None):
            self.p['theta0'] = theta0
        if not (phi0 == None):
            self.p['phi0'] = phi0
        if not (alpha0 == None):
            self.p['alpha0'] = alpha0
        if not (theta_m == None):
            self.p['theta_m'] = theta_m
        if not (phi_m == None):
            self.p['phi_m'] = phi_m
        print('Set angles:')
        self.printangles()
        
    def getangles(self):
        return (self.p['theta0'],self.p['phi0'],self.p['alpha0'], self.p['theta_m'],self.p['phi_m'])
        
    def printangles(self):
        print('theta offset: %s' % self.p['theta0'])
        print('phi offset: %s' % self.p['phi0'])
        print('rotation: %s' % self.p['alpha0'])
        print('theta manipulator: %s' % self.p['theta_m'])
        print('phi manipulator: %s' % self.p['phi_m'])        
        
           
    # from arpes.plotters import plot2D
    # from arpes.plotters import kplot2D
    # from arpes.plotters import plotEDC2D
    # from arpes.plotters import plotMDC2D
    from plotters import plot2D
    from plotters import kplot2D
    from plotters import plotEDC2D
    from plotters import plotMDC2D


    def show(self, mode=None, **kwargs):
        """General plotfunction
        Args:
            mode: string indicating which plotfunction to run. Run as 'help' for available modes
            kwargs: keyword args necessary to run the chosen plot function
        
        """
        
        warn("show is deprecated, use the direct functions prepended with plot")


        if mode is None:
            if(self.p['type'] == 'rawdata') or (self.p['type'] == 'angledata'):
                mode = 'col'
            else:
                mode = 'k'


        plotfuncdict = {
                        'angle' : self.plot2D,
                        'a' : self.plot2D,
                        'col' : self.plot2D,
                        'k': self.kplot2D,
                        'theta' : self.plot2D,
                        'EDC' : self.plotEDC2D,
                        'MDC' : self.plotMDC2D
                        }
                        
        if mode == 'help':
            for k,v in plotfuncdict.items():
                print('Available mode \'{}\', can be called directly as obj.{}, callsign and docstring: '.format(k,v.__name__))
                help(v)
                print('\n\n\n')
            return
        
                        
        if(not mode in plotfuncdict):
            raise ValueError('Unkown plotfunction mode')
            

        ax = plotfuncdict[mode](**kwargs) #run the plot function

        return ax
    
    def get_kax(self, mode = 'simple', output = 'abs'):
        """Calculate the kx, ky and kabs for the datagrid
        Args:
            mode: 'edep' or 'simple', 'edep' calculates an energy dependent grid, simple just a single k-correction at EF
            output: 'abs' will give a single dimensional kaxis along a curve, (only in simple) 
                        'kxky' will give kx and ky coords.
                        'kmesh' will give a tuple containing the kmesh and Emesh"""
        
        if self.p['type'] == 'rawdata':
            raise RuntimeError('Only analyzed data')
            
        if(mode == 'simple'):
            ## phi is set to zero, because manipulator offset is already captured by phi_m
            kx,ky = kfromangles(self.yaxis, 0., self.p['EkatEF'],Aaobj = None) 
            dkx = np.zeros_like(kx)
            dky = np.zeros_like(ky)
            dkx[1:] = kx[1:]-kx[:-1]
            dky[1:] = ky[1:]-ky[:-1]
            dk = np.sqrt(dkx*dkx + dky*dky)
            kax = np.cumsum(dk)
            argzero = np.argmin(np.sqrt(kx**2 + ky**2))
            kax -= kax[argzero]
            kmesh,Emesh = np.meshgrid(kax,self.xaxis)
            if output == 'abs':
                return kax
        elif(mode == 'edep'):
            thmesh,Emesh = np.meshgrid(self.yaxis,self.xaxis+self.p['EkatEF'])
            kx,ky = kfromangles(thmesh, 0., Emesh)
            dkx = np.zeros_like(kx).astype('float64') #this is to prevent a bug in np that raises an error taking the sqrt of 0
            dky = np.zeros_like(ky).astype('float64')
            dkx[:,1:] = kx[:,1:]-kx[:,:-1]
            dky[:,1:] = ky[:,1:]-ky[:,:-1]            
            dk = np.sqrt(dkx*dkx + dky*dky)
            kmesh = np.cumsum(dk, axis = 1)
            
            argzero = np.argmin(np.sqrt(kx**2 + ky**2), axis = 1)
            
            for i in range(kmesh.shape[0]):
                kmesh[i] -= kmesh[i,argzero[i]]
        if output == 'kxky':
            return kx,ky
        elif output == 'kmesh':
            return kmesh,Emesh
        else:
            raise RuntimeError('Uknown or invalid output mode')
            
            

    
    def getEDCs(self, path, pathtype = 'a', NEDC = 10, avg = 1, Estart = None, Eend = None):
        
        if not Estart:
            Estart = self.xaxis[0]
        if not Eend:
            Eend = self.xaxis[-1]

        if len(path) == 2:
            path = np.linspace(path[0],path[1],NEDC)
        else:
            raise ValueError('Invalid path')
        
        if pathtype == 'a':
            atbins = np.digitize(path,self.yaxis)
        elif pathtype == 'k':
            kax = self.get_kax(mode = 'simple', output = 'abs')
            atbins = np.digitize(path,kax)
        elif pathtype == 'kcor':
            raise NotImplementedError('kcor is not implemented')
        else:
            raise ValueError('Invalid path')
        
        Estartn,Eendn = np.digitize([Estart, Eend], self.xaxis)
        
        EDCs = []
        for i in range(NEDC):            
            startbin = atbins[i] - int(avg/2)
            endbin = startbin + avg
            EDCs.append(np.mean(self.data[Estartn:Eendn, max(startbin,0):min(endbin, self.shape[1]-1)], axis = 1))
                    
        Escale = self.xaxis[Estartn:Eendn]
        
        return EDCs, Escale, path
    

     
    def CastGold(self):
        if(not (self.p['type'] == 'rawdata')):
            raise RuntimeError('Only raw data can be cast to gold')

        self.__class__ = Gold
        self.p['analyzed'] = False
        self.p['slitcorrected'] = False


    
    def CastSlits(self, th_center = 0.):
        if(not (self.p['type'] == 'rawdata')):
            raise RuntimeError('Only raw data can be cast to slits')

        self.__class__ = Slits
        self.p['th_center'] = th_center
        self.p['analyzed'] = False


    def QuickCorrect(self, EF):
        """Quick method of correcting so k-plotting is available, only use offset, no calibration
        Args:
            EF: energy to which to set the Fermi level"""

        self.setXScale(self.xaxis[0] - EF, self.xaxis[1] - self.xaxis[0])
        self.p['type'] = 'corrected'
        self.p['EkatEF'] = EF
        
    def Correct(self, slitobj = None, goldobj = None):
 

        
        if slitobj: #correct slits
            if(not (self.p['type'] == 'rawdata')):
                raise RuntimeError('Only raw data can corrected')
        
            if(not slitobj.p['analyzed']):
                raise RuntimeError('Analyze slits first')
            
            #slit correction
            slitcorrected = SlitCorrect(self.data, slitobj.p['Tslit'])
            self.reshape(slitcorrected.shape)
            np.copyto(self.data, slitcorrected)
            self.setYScale(slitobj.p['Tyaxis'][0] + self.p['cryo_theta'], slitobj.p['Tyaxis'][1]-slitobj.p['Tyaxis'][0])       
            self.p['type'] = 'angledata'
            
        if goldobj:
            #gold correction        
            if(not goldobj.p['analyzed']):
                raise RuntimeError('Analyze gold first')
            goldcorrected, newxaxis = GoldCorrect(self.data, self.xaxis, goldobj.p['EFcurve'])
            self.setXScale(newxaxis[0], newxaxis[1]-newxaxis[0])
    
            #copy and finish up
            self.reshape(goldcorrected.shape)
            np.copyto(self.data,goldcorrected)   
            self.p['type'] = 'corrected'
            self.p['EkatEF'] = goldobj.p['EkatEF']

                    

            
    def save(self, path = '', objtype = 'AaData2D'):
        if path == '':
            if not 'savepath' in self.p:
                raise RuntimeError('No previously known savepath')
            else:
                path = self.p['savepath']
        else:
            if((path.find('.h5') == -1) and (path.find('.hdf5') == -1)):
                path += '.h5'
            self.p['savepath'] = path

                
        with h5py.File(path,'w') as hf:
            #header:
            hf.clear()
            hf.create_dataset('h/ver', data =  np.array('1.0'.encode('utf-8')))
            hf.create_dataset('h/obj', data =  np.array(objtype.encode('utf-8')))
            #data:    
            hf.create_dataset('d/data', data = self.data)
            hf.create_dataset('d/datalabel', data = np.array(self.datalabel.encode('utf-8')))
            hf.create_dataset('d/xstart', data = np.array(self.xstart))
            hf.create_dataset('d/xdelta', data = np.array(self.xdelta))
            hf.create_dataset('d/ystart', data = np.array(self.ystart))
            hf.create_dataset('d/ydelta', data = np.array(self.ydelta))
            hf.create_dataset('d/xlabel', data = np.array(self.xlabel.encode('utf-8')))
            hf.create_dataset('d/ylabel', data = np.array(self.ylabel.encode('utf-8')))
            hf.create_dataset('d/name', data = np.array(self.name.encode('utf-8')))
            hf.create_dataset('d/note', data = np.array(self.note.encode('utf-8')))                
            #parameters:            
            for k,v in self.p.items():
                if (type(v) == str) or (type(v) == np.str_):
                    hf.create_dataset('p/str/' + k, data = np.array(v.encode('utf-8')))
                elif v is None:                   
                    hf.create_dataset('p/none/' + k, data =  np.array(''.encode('utf-8'))) #save None's as an empty string
                elif type(v) is bool:                   
                    hf.create_dataset('p/bool/' + k, data =  np.array([v])) #save bools separately
                else:
                    hf.create_dataset('p/num/' + k, data =  np.array([v]))
                
                
    
            
        


class AaData3D(Data3D):
    def __init__(self, path, datatype = 'load', shape = (1376,1024, 1), name = '', note = '', zaxis = 'polar', **kwargs):
        super(AaData3D, self).__init__(shape = shape, name = name, note = note)
                


        if datatype == 'load':
            self.p = dl.loadH5_3D(self,path)
        elif datatype == 'empty':
            self.p = {}
        elif datatype == 'UBC':
            self.p = dl.loadData3DUBC(self,path = path)
            
        elif datatype == 'SLS':
            self.p = dl.loadData3DSLS(self,path = path)
            
        elif datatype == 'ALS':
            self.p = dl.loadData3DALS(self,path = path, zaxis = zaxis)    
        elif datatype == 'Maestro':
            self.p = dl.loadMaestro_3D(self,path = path)    
        elif datatype == 'APE':
            if zaxis == 'hv':
                self.p = dl.loadData3DAPE_hv(self,paths = path)
            elif zaxis == 'theta_y':
                self.p = dl.loadData3DAPE_FS(self,path)
        elif datatype == 'CLS':
            self.p = dl.loadData3DCLS(self, path = path,  **kwargs)
        else:
            raise RuntimeError('invalid datatype')


        
    def setangles(self, theta0 = None, phi0 = None, alpha0 = None, theta_m = None, phi_m = None):
        if not (theta0 == None):
            self.p['theta0'] = theta0
        if not (phi0 == None):
            self.p['phi0'] = phi0
        if not (alpha0 == None):
            self.p['alpha0'] = alpha0
        if not (theta_m == None):
            self.p['theta_m'] = theta_m
        if not (phi_m == None):
            self.p['phi_m'] = phi_m
        print('Set angles:')
        self.printangles()
        
    def getangles(self):
        return (self.p['theta0'],self.p['phi0'],self.p['alpha0'], self.p['theta_m'],self.p['phi_m'])
        
    def printangles(self):
        print('theta offset: %s' % self.p['theta0'])
        print('phi offset: %s' % self.p['phi0'])
        print('rotation: %s' % self.p['alpha0'])
        print('theta manipulator: %s' % self.p['theta_m'])
        print('phi manipulator: %s' % self.p['phi_m'])     
        

    ## append the plot functions defined in arpes.plotters        
    # from arpes.plotters import plotCE3D
    # from arpes.plotters import plotCEk3D
    # from arpes.plotters import plotcol3D
    # from arpes.plotters import plotlay3D
    # from arpes.plotters import plotkpath3D
    # from arpes.plotters import plotkpath2_3D
    # from arpes.plotters import plotEDC3D
    # from arpes.plotters import plotkz3D
    from plotters import plotCE3D
    from plotters import plotCEk3D
    from plotters import plotcol3D
    from plotters import plotlay3D
    from plotters import plotkpath3D
    from plotters import plotkpath2_3D
    from plotters import plotEDC3D
    from plotters import plotkz3D
            
    def show(self, mode=None, **kwargs):
        """General plotfunction
        Args:
            mode: string indicating which plotfunction to run. Run as 'help' for available modes
            kwargs: keyword args necessary to run the chosen plot function
        
        """
        
        warn("show is deprecated, use the direct functions prepended with plot")
        
        if mode is None:
            if(self.p['type'] == 'corrected'):
                mode = 'CEk'
            else:
                mode = 'row'
            if not ('val' in kwargs or 'num' in kwargs):
                if(self.p['type'] == 'corrected'):
                    kwargs['val'] = 0.0 # plot at EF
                else:
                    kwargs['num'] = int(0.7*self.xaxis.shape[0]) #plot at 70% of x
        
        plotfuncdict = {
                        'CE' : self.plotCE3D,
                        'CEk' : self.plotCEk3D,
                        'theta': self.plotcol3D,
                        'phi' : self.plotlay3D,
                        'kpath' : self.plotkpath3D,
                        'k' : self.plotkpath3D,
                        'k2' : self.plotkpath2_3D,
                        'EDC' : self.plotEDC3D,
                        'row' : self.plotCE3D,
                        'col' : self.plotcol3D,
                        'lay' : self.plotlay3D,
                        'kz' : self.plotkz3D
                        }
        
        if mode == 'help':
            for k,v in plotfuncdict.items():
                print('Available mode \'{}\', can be called directly as obj.{}, callsign and docstring: '.format(k,v.__name__))
                help(v)
                print('\n\n\n')
            return
        
                        
        if(not mode in plotfuncdict):
            raise ValueError('Unkown plotfunction mode')
            

        ax = plotfuncdict[mode](**kwargs) #run the plot function
            
        return ax

    def makeMDC(self, path, Eb, pathtype = 'k', pts = 200, Eint = None):
        """path should be [kx,ky] or [kx0,ky0,kx1,ky1]"""
        

        if Eint is None:
            Ebs = [Eb]
        else:
            Estart, Eend = round((Eb-Eint/2-self.xstart)/self.xdelta)*self.xdelta + self.xstart,round((Eb+Eint/2-self.xstart)/self.xdelta)*self.xdelta + self.xstart
            Ebs = np.arange(Estart,Eend+self.xdelta,self.xdelta)

        if pathtype == 'a':
            pass
        elif pathtype == 'k':
            if len(path) == 4: #
                kx = np.linspace(path[0],path[2],pts)
                ky = np.linspace(path[1],path[3],pts)
            else:
                kx = path[0]
                ky = path[1]
        MDC = np.zeros(kx.shape)
        for Eb in Ebs:
            Eb = np.linspace(Eb,Eb,kx.shape[0])
            MDC += self.atk(Eb,kx,ky, makegrid = False)
        return MDC
        
            

    def getEDCs(self,path,pathtype = 'k', NEDC = 10, avg = 1, avgtype = 'l', Estart = None, Eend = None):
        """Return EDCs in a list, together with the angles they are taken at and the energy scaling"""
        
        if not Estart:
            Estart = self.xaxis[0]
        if not Eend:
            Eend = self.xaxis[-1]
        
        if len(path) == 4:
            thpath = np.linspace(path[0],path[2],NEDC)
            phipath = np.linspace(path[1],path[3],NEDC)
        elif len(path) == 2:
            thpath = np.linspace(path[0][0],path[0][-1],NEDC)
            phipath = np.linspace(path[1][0],path[1][-1],NEDC)
        else:
            raise ValueError('Invalid path')
        
        if pathtype == 'a':
            pass
        elif pathtype == 'k':
            thpath, phipath = anglesfromk(thpath,phipath,self.p['EkatEF'], Aaobj = self)
        elif pathtype == 'phipath':
            phipath = np.linspace(phipath[0],phipath[0],NEDC)
        elif pathtype == 'thpath':
            thpath = np.linspace(thpath[0],thpath[0],NEDC)
        elif pathtype == 'kcor':
            raise NotImplementedError('kcor is not implemented')
        else:
            raise ValueError('Uknown path type')
        
        Estartn,Eendn = int((Estart-self.xstart)/self.xdelta),int((Eend-self.xstart)/self.xdelta)

        
        EDCs = []
        offsetavg = 0.0
        for i in range(NEDC):
            th,phi = thpath[i],phipath[i]
            nth,nphi = int((th-self.ystart)/self.ydelta),int((phi-self.zstart)/self.zdelta)
            if avgtype == 'sq': #square
                sqsize = int(np.ceil(np.sqrt(avg)))
                sqth = np.arange(sqsize) - int(sqsize/2) + nth
                sqphi = np.arange(sqsize) - int(sqsize/2) + nphi
                thpts,phipts = np.meshgrid(sqth,sqphi)
                thpts = thpts.flatten()
                phipts = phipts.flatten()
            elif avgtype == 'c': #circle
                dsize = int(np.ceil(np.sqrt(4*avg/np.pi)))
                circx = np.arange(dsize) - int(dsize/2)
                circy = np.arange(dsize) - int(dsize/2)
                circx,circy = np.meshgrid(circx,circy)
                circfilter = np.argwhere((circx*circx + circy*circy) < int(round((dsize+2)*(dsize+2)/4)))
                phipts = circy[circfilter[:,0],circfilter[:,1]] + nphi
                thpts = circx[circfilter[:,0],circfilter[:,1]] + nth
            elif avgtype == 'l': #line
                pts = np.arange(avg)-int(avg/2)
                thpts = (pts + nth) #make pts for th
                phipts = np.ones(pts.shape,dtype = 'int32')*nphi #make pts for phi
            filtpts = np.argwhere((thpts >= 0) & (thpts < self.yaxis.size) & (phipts >= 0) & (phipts < self.zaxis.size))
            if filtpts.size > avg:
                filtpts = filtpts[:avg]
            EDCs.append(np.sum(self.data[Estartn:Eendn,thpts[filtpts],phipts[filtpts]],axis=(1,2))/avg)
            offsetavg += np.sum(EDCs[i][:10])
            
        Escale = self.xaxis[Estartn:Eendn]
        return EDCs,Escale,thpath,phipath

    def QuickCorrect(self, EF):
        """Quick method of correcting so k-plotting is available, only use offset, no calibration
        Args:
            EF: energy to which to set the Fermi level"""

        self.setXScale(self.xaxis[0] - EF, self.xaxis[1] - self.xaxis[0])
        self.p['type'] = 'corrected'
        self.p['EkatEF'] = EF
        
    def Correct(self, slitobj = None, goldobj = None):
        
        
        
        if slitobj: #angle correction
            if(not slitobj.p['analyzed']):
                raise RuntimeError('Analyze slits first')
                
            if(not (self.p['type'] == 'rawdata')):
                raise RuntimeError('Only raw data can be corrected')
        

            self.setYScale(slitobj.p['Tyaxis'][0] + self.p['cryo_theta'][0], slitobj.p['Tyaxis'][1]-slitobj.p['Tyaxis'][0])    
        
            tmpdata = np.zeros((*slitobj.p['Tslit'].shape,self.data.shape[2]))
            for ii in range(self.data.shape[2]):  
                tmpdata[:,:,ii] = SlitCorrect(self.data[:,:,ii], slitobj.p['Tslit'])
            self.reshape(tmpdata.shape)
            np.copyto(self.data,tmpdata)
            self.p['type'] = 'angledata'


        if goldobj: #gold correction            
            if(not goldobj.p['analyzed']):
                raise RuntimeError('Analyze gold first')
                
            for ii in range(self.data.shape[2]):    
                self.data[:,:,ii], newxaxis = GoldCorrect(self.data[:,:,ii], self.xaxis, goldobj.p['EFcurve'])                
            self.setXScale(newxaxis[0], newxaxis[1]-newxaxis[0])
            self.p['type'] = 'corrected'
            self.p['EkatEF'] = goldobj.p['EkatEF']

                            


    def glue(self, obj, axis = 'lay'):
        """glue obj to self along axis, self will determine the final scaling"""

        if((axis == 'col') or (axis == 'kx') or (axis == 'y') or (axis == 'th') or (axis == 'pol')):
            axis = 'col'
        elif((axis == 'lay') or (axis == 'ky') or (axis == 'z') or (axis == 'phi') or (axis == 'azi')):
            axis = 'lay'
        else:
            raise ValueError('Invalid axis')

        #check if obj is the right size:
        if(not((self.shape[0] == obj.shape[0]) and (((axis == 'lay') and (self.shape[1] == obj.shape[1])) or ((axis == 'col') and (self.shape[2] == obj.shape[2]))))):
            raise RuntimeError('Incompatible array shapes')

        #determine new size
        if(axis == 'lay'):
            self.reshape((self.shape[0], self.shape[1], self.shape[2]+obj.shape[2]))
            self.data[:,:,self.shape[2]-obj.shape[2]:] = obj.data
        else:
            self.reshape((self.shape[0], self.shape[1]+obj.shape[1], self.shape[2]))
            self.data[:,self.shape[1]-obj.shape[1]:,:] = obj.data
        
            

    
    
    def stitch(self, obj, val, axis = 'col', ratio = 1):
        """Stitch data together from two blocks on axis at val, ratio is for normalization"""
        
        if(not((self.shape[0] == obj.shape[0]) and (((axis == 'lay') and (self.shape[1] == obj.shape[1])) or ((axis == 'col') and (self.shape[2] == obj.shape[2]))))):
            raise RuntimeError('Incompatible array shapes')

        if(not ((axis == 'col') or (axis == 'lay'))):
            raise ValueError('Invalid axis')

        #find crossover point
        if(axis == 'col'):
            numself = np.argwhere(self.yaxis < val)[-1][0]
            numobj = np.argwhere(obj.yaxis >= val)[0][0]
            axisnum = 1
        else:
            numself = np.argwhere(self.zaxis < val)[-1][0]
            numobj = np.argwhere(obj.zaxis >= val)[0][0]
            axisnum = 2

        
        if((numself >= self.shape[axisnum]) or (numself < 0) or (numobj < 0) or (numobj >= obj.shape[axisnum])):
            raise ValueError('val out of range')
        
        if(axis == 'col'):
            self.reshape((self.shape[0], numself + obj.shape[1] - numobj + 1, self.shape[2]))
            self.data[:,numself+1:,:] = obj.data[:,numobj:,:]
            self.yaxis[numself+1:] = obj.yaxis[numobj:]
        else:
            self.reshape((self.shape[0],self.shape[1], numself +  obj.shape[2] - numobj + 1))
            self.data[:,:,numself+1:] = obj.data[:,:,numobj:]
            self.zaxis[numself+1:] = obj.zaxis[numobj:]


    def getZSlice(self, z, xmin = None, xmax = None, ymin = None, ymax = None):
        """ get slice in the layer (z) direction
        
        Args:
            z: z value to slice at
            xmin: crop window minimum for x
            xmax: crop window maximum for x
            ymin: etc
            
            
        Returns:
            cropped AaData2D object at z, with parameters at z
        
        """
        if xmin is None:
            xmin = self.xstart
        if xmax is None:
            xmax = self.xstart + self.xaxis.size * self.xdelta
        if ymin is None:
            ymin = self.ystart
        if ymax is None:
            ymax = self.ystart + self.yaxis.size * self.ydelta
            
        #get min_n:
        xmin_n = int((xmin-self.xstart)/self.xdelta)
        xmax_n = int((xmax-self.xstart)/self.xdelta)
        ymin_n = int((ymin-self.ystart)/self.ydelta)
        ymax_n = int((ymax-self.ystart)/self.ydelta)
        
        zn = np.absolute(self.zaxis-z).argmin()
                        
        newobj = AaData2D('', datatype = 'empty', name = self.name, note = self.note)
        newobj.reshape((xmax_n-xmin_n, ymax_n-ymin_n))
        newobj.data = self.data[xmin_n:xmax_n,ymin_n:ymax_n,zn]
        newobj.p = self.p.copy()
        
        for k,v in newobj.p.items():
            if type(v) == np.ndarray:
                if v.size == self.shape[2]:
                    newobj.p[k] = v[zn]
                else:
                    newobj.p[k] = v
            else:
                newobj.p[k] = v

        newobj.setXScale(self.xaxis[xmin_n], self.xdelta)
        newobj.setYScale(self.yaxis[ymin_n], self.ydelta)

        return newobj

                                
    def crop(self, xmin = None, xmax = None, ymin = None, ymax = None, zmin = None, zmax = None):
        """ crop AaData with set maxima
        
        Args:
            xmin: crop window minimum for x
            xmax: crop window maximum for x
            ymin: etc
            
        Returns:
            cropped AaData3D object with the same parameters
        
        """
        if xmin is None:
            xmin = self.xstart
        if xmax is None:
            xmax = self.xstart + self.xaxis.size * self.xdelta
        if ymin is None:
            ymin = self.ystart
        if ymax is None:
            ymax = self.ystart + self.yaxis.size * self.ydelta
        if zmin is None:
            zmin = self.zstart
        if zmax is None:
            zmax = self.zstart + self.zaxis.size * self.zdelta
            
        #get min_n:
        xmin_n = int((xmin-self.xstart)/self.xdelta)
        xmax_n = int((xmax-self.xstart)/self.xdelta)
        ymin_n = int((ymin-self.ystart)/self.ydelta)
        ymax_n = int((ymax-self.ystart)/self.ydelta)
        zmin_n = int((zmin-self.zstart)/self.zdelta)
        zmax_n = int((zmax-self.zstart)/self.zdelta)
        
        print(xmin_n,xmax_n,ymin_n,ymax_n,zmin_n,zmax_n)
        
        return self.cropN(xmin_n,xmax_n,ymin_n,ymax_n,zmin_n,zmax_n)
        
    
    def cropN(self, xmin_n = None, xmax_n = None, ymin_n = None, ymax_n = None, zmin_n = None, zmax_n = None):
        if xmin_n is None:
            xmin_n = 0
        if xmax_n is None:
            xmax_n = self.xaxis.size
        if ymin_n is None:
            ymin_n = 0
        if ymax_n is None:
            ymax_n = self.yaxis.size
        if zmin_n is None:
            zmin_n = 0
        if zmax_n is None:
            zmax_n = self.zaxis.size
            
        xmin_n,xmax_n = (min(xmin_n,xmax_n),max(xmin_n,xmax_n))
        ymin_n,ymax_n = (min(ymin_n,ymax_n),max(ymin_n,ymax_n))
        zmin_n,zmax_n = (min(zmin_n,zmax_n),max(zmin_n,zmax_n))
            
        newobj = AaData3D('', datatype = 'empty', shape = (xmax_n-xmin_n, ymax_n-ymin_n, zmax_n-zmin_n), name = self.name, note = self.note)
        newobj.data = np.array(self.data[xmin_n:xmax_n,ymin_n:ymax_n,zmin_n:zmax_n])
        newobj.p = self.p.copy()
        
        newobj.xaxis = self.xaxis[xmin_n:xmax_n]
        newobj.yaxis = self.yaxis[ymin_n:ymax_n]
        newobj.zaxis = self.zaxis[zmin_n:zmax_n]
        newobj.xstart = newobj.xaxis[0]
        newobj.ystart = newobj.yaxis[0]
        newobj.zstart = newobj.zaxis[0]
        newobj.xdelta = newobj.xaxis[1] - newobj.xaxis[0]
        newobj.ydelta = newobj.yaxis[1] - newobj.yaxis[0]
        newobj.zdelta = newobj.zaxis[1] - newobj.zaxis[0]

        
        return newobj
        
    def norm(self, axis = 0, fitpoly = None):
        axes = [0,1,2]
        axes.remove(axis)
        normfac = 1
        for ax in axes:
            normfac *= self.data.shape[ax]
        normvec = np.sum(self.data, axis = tuple(axes))/normfac
        
        if fitpoly:
            #fit poly of order fitpoly and return
            pass
        
        shape = self.data.shape
        
        newobj = AaData3D('', datatype = 'empty', shape = shape, name = self.name, note = self.note)
        newobj.p = self.p
        newobj.setXScale(self.xstart,self.xdelta)
        newobj.setYScale(self.ystart,self.ydelta)
        newobj.setZScale(self.zstart,self.zdelta)
        
        newobj.data = (self.data.reshape((shape[axes[0]],shape[axes[1]],shape[axis]))/normvec).reshape((shape[0],shape[1],shape[2]))
        
        return newobj
        
        
    def rebin(self, binning = (1,1,1)):
        
        shape = self.data.shape  
        #reduce to nearest integer multiple
        redshape = (int(shape[0]//binning[0])*binning[0],int(shape[1]//binning[1])*binning[1],int(shape[2]//binning[2])*binning[2])

        self.reshape(redshape)

        newshape = (int(redshape[0]//binning[0]),int(redshape[1]//binning[1]),int(redshape[2]//binning[2]))
        compression_pairs = [(d, int(c//d)) for d,c in zip(newshape,redshape)]

        flattened = [l for p in compression_pairs for l in p]

        newdata = self.data.reshape(flattened)
        for i in range(len(newshape)):
            newdata = newdata.mean(-1*(i+1))
        
        newobj = AaData3D('', datatype = 'empty', shape = newshape, name = self.name, note = self.note)
        newobj.p = self.p.copy()
        newobj.setXScale(self.xstart,binning[0]*self.xdelta)
        newobj.setYScale(self.ystart,binning[1]*self.ydelta)
        newobj.setZScale(self.zstart,binning[2]*self.zdelta)        
        
        newobj.data = newdata
        
        return newobj
    
    
    
    def atkCE(self, kx, ky, Eb):
        """gives the intensity at kx,ky,Eb. kx and ky can be a number vector and number
        or a meshgrid, when mesh is set to True a meshgrid is created from kx and ky"""
        if(not (self.p['type'] == 'corrected')):
            raise RuntimeError('k interpolation is only available for corrected data')
        
        Enum = int(round((Eb-self.xstart)/self.xdelta))
        
        if((Enum < 0 ) | (Enum >= self.xaxis.shape[0])):
            raise ValueError('Energy value out of range')
    

        if(not (type(kx) == np.ndarray)):
            kx = np.array([kx])
        if(not (type(ky) == np.ndarray)):
            ky = np.array([ky])
        
        
        th,phi = anglesfromk(kx,ky,self.p['EkatEF']+Eb,Aaobj = self)
    
        thmin = min(self.yaxis[0],self.yaxis[-1])
        thmax = max(self.yaxis[0],self.yaxis[-1])
        phimin = min(self.zaxis[0],self.zaxis[-1])
        phimax = max(self.zaxis[0],self.zaxis[-1])
    
        
    
        #interpolator:
        interp = RectBivariateSpline(self.yaxis, self.zaxis, self.data[Enum])
        
        dataout = interp(kx, ky, grid = False)
        dataout[np.where((th < thmin) | (th > thmax) | (phi < phimin) | (phi > phimax))] = 0

        return dataout

        
    def atk(self, E, kx, ky, makegrid = True):
        """return k corrected slice of the data with coordinates given by Eb,kx,ky (can be numbers, paths, meshgrids)"""

        if(not (self.p['type'] == 'corrected')):
            raise RuntimeError('k interpolation is only available for corrected data')
        
        if makegrid:
            if(not (type(kx) == np.ndarray)):
                kx = np.array([kx])
            if(not (type(ky) == np.ndarray)):
                ky = np.array([ky])
            if(not (type(E) == np.ndarray)):
                E = np.array([E])
            gridE,gridkx,gridky = np.meshgrid(E,kx,ky)
        else: #if providing own grid
            gridE = E
            gridkx = kx
            gridky = ky
            
        gridth,gridphi = anglesfromk(gridkx,gridky,self.p['EkatEF']+gridE,Aaobj = self)

        if len(gridE.shape) == 3:
            intgrid = np.array([gridE,gridth,gridphi]).transpose((2,1,3,0))
        elif len(gridE.shape) == 2:
            intgrid = np.array([gridE,gridth,gridphi]).transpose((1,2,0))
        elif len(gridE.shape) == 1:
            intgrid = np.array([gridE,gridth,gridphi]).transpose((1,0))

        return interpn((self.xaxis,self.yaxis,self.zaxis),self.data, intgrid, bounds_error = False, fill_value = 0.)
        

    def atkz(self, E, kx, kz, makegrid = True):
        """return k corrected slice of the data with coordinates given by Eb,kx,kz (can be numbers, paths, meshgrids)
        assumes zaxis is kz data"""
        
        if(not (self.p['type'] == 'corrected')):
            raise RuntimeError('k interpolation is only available for corrected data')
        
        if makegrid:
            if(not (type(kx) == np.ndarray)):
                kx = np.array([kx])
            if(not (type(kz) == np.ndarray)):
                kz = np.array([kz])
            if(not (type(E) == np.ndarray)):
                E = np.array([E])
            gridE,gridkx,gridkz = np.meshgrid(E,kx,kz)
        else: #if providing own grid
            gridE = E
            gridkx = kx
            gridkz = kz      
            
        if not 'V0' in self.p:
            self.p['V0'] = 10.
            print('V0 not found, added std value of 10.')
        if not 'WF' in self.p:
            self.p['WF'] = 4.3
            print('WF not found, added std value of 4.3')
            
            
        gridth,gridhv = anglesfromkz(gridkx,gridkz,gridE, Aaobj = self)
        
        if len(gridE.shape) == 3:
            intgrid = np.array([gridE,gridth,gridhv]).transpose((2,1,3,0))
        elif len(gridE.shape) == 2:
            intgrid = np.array([gridE,gridth,gridhv]).transpose((1,2,0))
        elif len(gridE.shape) == 1:
            intgrid = np.array([gridE,gridth,gridhv]).transpose((1,0))
        ##check order of zaxis:
        if self.zaxis[1] > self.zaxis[0]:
            return interpn((self.xaxis,self.yaxis,self.zaxis),self.data, intgrid, bounds_error = False, fill_value = 0.)
        else:
            return interpn((self.xaxis,self.yaxis,self.zaxis[::-1]),self.data[:,:,::-1], intgrid, bounds_error = False, fill_value = 0.)
        
        
                    
    def save(self, path = ''):
        if path == '':
            if not 'savepath' in self.p:
                raise RuntimeError('No previous savepath found')
            else:
                path = self.p['savepath']
        else:
            if((path.find('.h5') == -1) and (path.find('.hdf5') == -1)):
                path += '.h5'
            self.p['savepath'] = path

                
        with h5py.File(path,'w') as hf:
            #header:
            hf.clear()
            hf.create_dataset('h/ver', data =  np.array('1.0'.encode()))
            hf.create_dataset('h/obj', data =  np.array('AaData3D'.encode()))
            #data:    
            hf.create_dataset('d/data', data = self.data)
            hf.create_dataset('d/datalabel', data = np.array(self.datalabel.encode()))
            hf.create_dataset('d/xstart', data = np.array(self.xstart))
            hf.create_dataset('d/xdelta', data = np.array(self.xdelta))
            hf.create_dataset('d/ystart', data = np.array(self.ystart))
            hf.create_dataset('d/ydelta', data = np.array(self.ydelta))
            hf.create_dataset('d/zstart', data = np.array(self.zstart))
            hf.create_dataset('d/zdelta', data = np.array(self.zdelta))
            hf.create_dataset('d/xlabel', data = np.array(self.xlabel.encode()))
            hf.create_dataset('d/ylabel', data = np.array(self.ylabel.encode()))
            hf.create_dataset('d/zlabel', data = np.array(self.ylabel.encode()))
            hf.create_dataset('d/name', data = np.array(self.name.encode()))
            hf.create_dataset('d/note', data = np.array(self.note.encode()))                
            #parameters:            
            for k,v in self.p.items():
                if (type(v) == str) or (type(v) == np.str_):
                    hf.create_dataset('p/str/' + k, data = np.array(v.encode('utf-8')))
                elif v is None:                   
                    hf.create_dataset('p/none/' + k, data =  np.array(''.encode('utf-8'))) #save None's as an empty string
                elif type(v) is bool:                   
                    hf.create_dataset('p/bool/' + k, data =  np.array([v])) #save bools separately
                else:
                    hf.create_dataset('p/num/' + k, data =  np.array([v]))
                
                
            
            


class Slits(AaData2D):
    def __init__(self, path = '', datatype = None, shape = (1376,1024), name = '', note = ''):        

        if datatype is None:
            if ((path[-3:] == '.h5') or (path[-5:] == '.hdf5')): #if h5 in name, default to load saved data
                datatype = 'load'
            else:
                datatype = 'UBC' #else default to UBC load

        super(Slits, self).__init__(path = path, objtype = 'Slits', datatype = datatype, shape = shape, name = name, note = note)
            
            
        if datatype == 'UBC':
            self.p['Tslit'] = np.zeros(shape) #transformation matrix
            self.p['Txaxis'] = np.array(self.xaxis)
            self.p['Tyaxis'] = np.zeros(shape[1])
            self.p['curveparams'] = None
            self.p['slitcurves'] = None
            self.p['th_center'] = 0. # center of the transformation is zero
            self.p['analyzed'] = False
            self.p['mask'] = None # list of (row,col) tuples that define the mask
            self.p['usemask'] = False
        

        
    def setmask(self, mask = None, defaultmask = False): # inc standard shaped mask
        if defaultmask:
            mask = [(1375,250), (200,190), (200,850),(1375,794)]
        self.p['mask'] = np.array(mask)
        self.p['usemask'] = (mask is not None)

    
    def getmask(self):
        return [(v1,v2) for (v1,v2) in self.p['mask']]
    
    def makemask(self):        
        #define the polygon
        p = Path(self.p['mask'])
        
        #generate all array indices:
        idx_x = np.arange(self.data.shape[0])
        idx_y = np.arange(self.data.shape[1])
        IDX_y,IDX_x = np.meshgrid(idx_y,idx_x)
        idx = np.array([IDX_x.flatten(),IDX_y.flatten()]).T     
        
        mask = p.contains_points(idx).reshape(self.data.shape)
        
        return mask
    
    
    def showmask(self):
        
        if not self.p['mask']:
            raise RuntimeError('no mask set')

        if self.p['usemask']:
            plotdata = ma.filled(ma.array(self.data, mask = self.makemask()), fill_value = 0.)
        else: 
            plotdata = self.data
        
        #plot mask
        maskline = np.empty((len(self.p['mask'])+1,2))
        maskline[0:-1] = self.p['mask']
        maskline[-1] = maskline[0]

        fig, ax = plt.subplots(figsize=(6,6))
        ax.plot(maskline.T[1],maskline.T[0],color='red')
        ax.imshow(np.flipud(plotdata), cmap='Spectral_r', interpolation = 'none')
        plt.show()   



    def Analyze(self, N = None):
        """Analyze Slits with N lines"""
        
        if self.p['usemask']:
            slitdata = ma.filled(ma.array(self.data, mask = self.makemask()), fill_value = 0.)
        else: 
            slitdata = self.data

        self.p['Tslit'], self.p['Txaxis'], self.p['Tyaxis'], self.p['curveparams'], self.p['slitcurves'] = AnalyzeSlits(slitdata, N = N)
        self.p['Tyaxis'] += self.p['th_center']
        self.p['analyzed'] = True
        
    
    def showslits(self):
        if(not self.p['analyzed']):
            raise RuntimeError('Analyze slits first')
            
        ystart = self.ystart
        yend = self.ystart+(self.data.shape[1]-1)*self.ydelta
        xstart = self.xstart
        xend = self.xstart+(self.data.shape[0]-1)*self.xdelta
        asp = self.data.shape[0]/float(self.data.shape[1]) # actual aspect ratio that is shown for image asp = vertical/horizontal
        ratio = (yend-ystart)/(xend-xstart)*asp
        fig, ax = plt.subplots(figsize=(6,6))
        for ii in range(self.p['slitcurves'].shape[0]):
            ax.plot(self.p['slitcurves'][ii], self.xaxis, color = 'red')
        ax.imshow(np.flipud(self.data), cmap='Spectral_r', interpolation = 'none', extent = [ystart,yend,xstart,xend], aspect=ratio)
        
        plt.show()        
        

    def Correct(self, datain):
        corrected = AaData2D(shape = self.shape)
        corrected.Ek = datain.Ek
        corrected.Ep = datain.Ep
        corrected.paramdict = datain.paramdict
        corrected.lens_type = datain.lens_type
        corrected.slit_entrance = datain.slit_entrance

        corrected.setXScale(datain.xaxis[0], datain.xaxis[1]-datain.xaxis[0])
        corrected.setYScale(self.p['Tyaxis'][0], self.p['Tyaxis'][1]-self.p['Tyaxis'][0])

        np.copyto(corrected,SlitCorrect(datain, self.p['Tslit']))        
                    
        return corrected
        
    def save(self, path = ''):
        #can be removed by impletenting comparing to the obj type in AaData2D
        super(Slits, self).save(path = path, objtype = 'Slits')

        
def UBCSlits(*args, **kwargs):
    warn("UBCSlits is deprecated, use Slits instead")
    return Slits(*args, **kwargs)
    

        
        
        
        
class Gold(AaData2D):
    def __init__(self, path, objtype = 'Gold', datatype = None, name = '', note = '', mode = 'avg', pin = dl.std_p()):

        if datatype is None:
            if ((path[-3:] == '.h5') or (path[-5:] == '.hdf5')): #if h5 in name, default to load saved data
                datatype = 'load'
            else:
                datatype = 'UBC' #else default to UBC load
        
                 
        super(Gold, self).__init__(path = path, objtype = objtype, datatype = datatype, name = name, note = note, mode = mode, pin = pin)

 
        if datatype == 'UBC':
            self.p['EFparams'] = np.array([()])
            self.p['EFcurve'] = np.zeros(self.shape[1]) # center of the transformation is zero
            self.p['analyzed'] = False
            self.p['slitcorrected'] = False
        elif datatype == 'SLS' or datatype == 'ALS':
            self.p['EFparams'] = np.array([()])
            self.p['EFcurve'] = np.zeros(self.shape[1])
            self.p['analyzed'] = False
            self.p['slitcorrected'] = True


    def Correct(self, slitobj): #correction that produces gold object from raw gold and slits

        if(not slitobj.p['analyzed']):
            raise RuntimeError('Analyze slits first')
        
        if(self.p['slitcorrected']):
            raise RuntimeError('Already slit corrected')

        self.data = SlitCorrect(self.data, slitobj.p['Tslit'])  


        self.p['type'] = 'corrected'
        self.p['slitcorrected'] = True

        self.setYScale(slitobj.p['Tyaxis'][0], slitobj.p['Tyaxis'][1]-slitobj.p['Tyaxis'][0])

                            


    def Analyze(self, T = None, Np = 2, region = None, xmin = None, xmax = None, ymin = None, ymax = None): #T is temp, Np is polynomial order
        """Analyze Gold: intensity and energy correction
        Args:
            T: temperature in K
            Np: polynomial order to fit to gold edge
            region: tuple of region indices (xmin_n, xmax_n, ymin_n,ymax_n)
            xmin, xmax, ymin, ymax: region values"""
        if region is None:
            region = []
            region.append(int(np.digitize(xmin, self.xaxis)) if xmin is not None else 0)
            region.append(int(np.digitize(xmax, self.xaxis)) if xmax is not None else self.shape[0])
            region.append(int(np.digitize(ymin, self.yaxis)) if ymin is not None else 0)
            region.append(int(np.digitize(ymax, self.yaxis)) if ymax is not None else self.shape[1])
            region = tuple(region)
                
        if(self.p['slitcorrected'] == False):
            raise RuntimeError('Can only analyze slit corrected data')
        if(T == None):
            T = self.p['T_sample']
        fitdict, self.p['EFparams'], self.p['EFcurve'] = AnalyzeGold(self.xaxis, self.data, T, Np = Np, region = region)
        self.p['EFparams'] = np.array([self.p['EFparams']])
        self.p['analyzed'] = True
        self.p['EkatEF'] = np.mean(self.p['EFcurve'])
        for key in fitdict.keys(): #expand fitdict into p
            self.p['fitres_' + key] = fitdict[key]


    def showgold(self):
        if(not self.p['analyzed']):
            raise RuntimeError('Can only show analyzed data')
        ystart = self.ystart
        yend = self.ystart+(self.data.shape[1]-1)*self.ydelta
        xstart = self.xstart
        xend = self.xstart+(self.data.shape[0]-1)*self.xdelta
        asp = self.data.shape[0]/float(self.data.shape[1]) # actual aspect ratio that is shown for image asp = vertical/horizontal
        ratio = (yend-ystart)/(xend-xstart)*asp
        fig, ax = plt.subplots(figsize=(6,6))

        #add the found EF with resolution as errorbars in black squares
        ax.imshow(np.flipud(self.data), cmap='Spectral_r', interpolation = 'none', extent = [ystart,yend,xstart,xend], aspect=ratio)
        ax.plot(self.yaxis[np.array(self.p['fitres_col'], dtype = 'int64')], (self.p['fitres_EF']-self.p['fitres_sig']), color = 'black', linestyle = 'none', marker = '_')
        ax.plot(self.yaxis[np.array(self.p['fitres_col'], dtype = 'int64')], (self.p['fitres_EF']+self.p['fitres_sig']), color = 'black', linestyle = 'none', marker = '_')
        ax.plot(self.yaxis, self.p['EFcurve'], color = 'red')
        
        plt.show()     

    def save(self, path = ''):
        #can be removed by impletenting comparing to the obj type in AaData2D
        super(Gold, self).save(path = path, objtype = 'Gold')
        
        
    def getres(self, region = None, T = None):
        if(not self.p['analyzed']):
            raise RuntimeError('Can only get resolution of analyzed data')
        if region is None:
            region = (0,self.xaxis.shape[0],int(0.4*self.shape[1]),int(0.6*self.shape[1])) #middle fifth
        if T is None:
            T = self.p['T_sample']
        goldcorrected, newxaxis = GoldCorrect(self.data, self.xaxis, self.p['EFcurve'])
        
        print(region)
        fitdata = np.mean(goldcorrected[region[0]:region[1],region[2]:region[3]], axis = 1)
        
        
        
        fitres = FDexpfit(newxaxis[region[0]:region[1]], fitdata, T)
        
        print('Resolution is %0.1f mV' % (1000*fitres[2]))
        
        self.p['Res Fit'] = fitres

        fig, ax = plt.subplots(figsize=(6,3.5))
        ax.plot(newxaxis[region[0]:region[1]],fitdata, color = 'black', linestyle = 'none', marker = '.')
        ax.plot(newxaxis[region[0]:region[1]],FDexp(newxaxis[region[0]:region[1]], *fitres), color = 'red', label ='fit')
        ax.vlines(0.0,-2*np.max(fitdata),2*np.max(fitdata), color = 'red')
        ax.vlines(fitres[2]/2,-2*np.max(fitdata),2*np.max(fitdata), color = 'green')
        ax.vlines(-fitres[2]/2,-2*np.max(fitdata),2*np.max(fitdata), color = 'green', label = 'exp fwhm')
        ax.axis([newxaxis[region[0]],newxaxis[region[1]],0.0,1.2*np.max(fitdata)])
        ax.set_xlabel('Energy[eV]', fontsize = 12)
        ax.set_ylabel('Intensity []', fontsize = 12)
        ax.text(fitres[2],0.8*np.max(fitdata), ('Resolution %0.1f meV'  % (1000*fitres[2])))
        plt.show()

    
def UBCGold(*args, **kwargs):
    warn("UBCGold is deprecated, use Gold instead")
    return Gold(*args, **kwargs)
    
            
class Gold1D(AaData1D):
    pass 

    def __init__(self, path = '', datatype = 'load', name = '', note = ''):

        
        super(Gold1D, self).__init__(path = path, objtype = 'Gold1D', datatype = datatype, name = name, note = note)
 
        self.p['analyzed'] = False

                            


    def Analyze(self, T = None, Np = 2, region = None): #T is temp, Np is polynomial order
        """Analyze Gold: intensity and energy correction"""

        if(T == None):
            T = self.p['T_sample']        
      
        self.p['fitres'] = AnalyzeGold1D(self.xaxis, self.data, T, region)
        self.p['EkatEF'] = self.p['fitres'][0]
        self.p['analyzed'] = True

 

    def save(self, path = ''):
        #can be removed by impletenting comparing to the obj type in AaData2D
        super(Gold, self).save(path = path, objtype = 'Gold')
        
        
    def showgold(self, region = None, T = None):
        if(not self.p['analyzed']):
            raise RuntimeError('Can only show analyzed data')

        print('Resolution is %s mV' % (1000*self.p['fitres'][2]))
        
        fig, ax = plt.subplots(figsize=(6,3.5))
        ax.plot(self.xaxis,self.data, color = 'black', linestyle = 'none', marker = '.')
        ax.plot(self.xaxis,FDexp(self.xaxis, *self.p['fitres']), color = 'red')
        ax.vlines(0.0,-2*np.max(self.data),2*np.max(self.data), color = 'red')
        ax.vlines(self.p['fitres'][2]/2,-2*np.max(self.data),2*np.max(self.data), color = 'green')
        ax.vlines(-self.p['fitres'][2]/2,-2*np.max(self.data),2*np.max(self.data), color = 'green')
        ax.axis([self.xaxis[0],self.xaxis[-1],0.0,1.2*np.max(self.data)])
        ax.set_xlabel('Energy[eV]', fontsize = 12)
        ax.set_ylabel('Intensity []', fontsize = 12)    
        plt.show()
        
        
        
        
            
class GoldDummy():
    def __init__(self, Ek, size = 1024):
        self.p = {}

        self.p['EkatEF'] = Ek
        self.p['EFcurve'] = np.linspace(Ek,Ek, size)
        self.p['analyzed'] = True
    

class SpinPol():
    def __init__(self, up, dn, sens_up = 1.0, sens_dn = 1.0):
        self.data = {'up' : up, 'dn' : dn}
        self.p = {'sens_up' : sens_up, 'sens_dn' : sens_dn} ## possibly change this to avoid confusion

    # from arpes.plotters import plotspinpol
    # from arpes.plotters import plotspin
    from plotters import plotspinpol
    from plotters import plotspin

    def show(self, mode='spin', **kwargs):
        """General plotfunction
        Args:
            mode: string indicating which plotfunction to run. Run as 'help' for available modes
            kwargs: keyword args necessary to run the chosen plot function
        
        """
        
        warn("show is deprecated, use the direct functions prepended with plot")
                
        plotfuncdict = {
                        'pol' : self.plotspinpol,
                        'spin' : self.plotspin,
                        }
        
        if mode == 'help':
            for k,v in plotfuncdict.items():
                print('Available mode \'{}\', can be called directly as obj.{}, callsign and docstring: '.format(k,v.__name__))
                help(v)
                print('\n\n\n')
            return
        
                        
        if(not mode in plotfuncdict):
            raise ValueError('Unkown plotfunction mode')
            

        ax = plotfuncdict[mode](**kwargs) #run the plot function
            
        return ax
    
    
    
class CPSA():
    def __init__(self, CP_up, CP_dn, CM_up, CM_dn, sens_up = 1.0, sens_dn = 1.0, sens_CP = 1.0, sens_CM = 1.0):
        self.data = {'CP_up' : CP_up, 'CP_dn' : CP_dn, 'CM_up' : CM_up, 'CM_dn' : CM_dn}
        self.p = {'sens_up' : sens_up, 'sens_dn' : sens_dn, 'sens_CP' : sens_CP, 'sens_CM' : sens_CM} ## possibly change this to avoid confusion

    # from arpes.plotters import plotcpspol
    # from arpes.plotters import plotcoan
    from plotters import plotcpspol
    from plotters import plotcoan

    def show(self, mode='cps', **kwargs):
        """General plotfunction
        Args:
            mode: string indicating which plotfunction to run. Run as 'help' for available modes
            kwargs: keyword args necessary to run the chosen plot function
        
        """
        
        warn("show is deprecated, use the direct functions prepended with plot")
                
        plotfuncdict = {
                        'cps' : self.plotcpspol,
                        'coan' : self.plotcoan,
                        }
                        
        
        if mode == 'help':
            for k,v in plotfuncdict.items():
                print('Available mode \'{}\', can be called directly as obj.{}, callsign and docstring: '.format(k,v.__name__))
                help(v)
                print('\n\n\n')
            return
        
                        
        if(not mode in plotfuncdict):
            raise ValueError('Unkown plotfunction mode')
            

        ax = plotfuncdict[mode](**kwargs) #run the plot function
            
        return ax

