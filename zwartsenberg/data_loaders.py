# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 13:10:51 2017

@author: Berend
"""

#loaders for different beamlines

import numpy as np
import os
import h5py
import time
import datetime
import igor.binarywave
import struct


def check_keys(p):
    """Checks if keys in p contain forward slashes, causing issues saving h5 data
    args:
        p : the parameter dictionary to be checked
    returns:
        None
    """
    
    for k in p.keys():
        if k.find('/') > 0:
            print('Forward slash found in key \'%s\'. Saving H5 data will give problems' % k)


def std_p():
    p = {}
    p['theta0'] = 0. #sample surface offset angles
    p['phi0'] = 0.
    p['alpha0'] = 0. #sample rotation
    p['theta_m'] = 0. #manipulator offset
    p['phi_m'] = 0.
    p['EkatEF'] = None #kinetic energy of the gold at EF
    p['bz'] = 'square' #support more than one?
    p['alatt'] = 3.85 # in angstrom
    p['sample'] = ''
    p['compound'] = ''
    return p      



######################### 1D loaders ##########################################
def loadData1DAPE(Aa1D,path,pin = std_p()):
    p = pin.copy()

    loaddict = igor.binarywave.load(path)
    #loadpxp: 
    #import igor.packed
    #igor.packed.load
    wave = loaddict['wave']
    p['version'] = loaddict['version'] #save version number
    
    #incorporate loading different regions from the same file
    rows,cols = wave['wave_header']['nDim'][0:2]

    Aa1D.reshape((rows,))
    Aa1D.data = np.array(wave['wData'], dtype = 'float32')
    deltas,offsets = wave['wave_header']['sfA'],wave['wave_header']['sfB']
    Aa1D.setXScale(offsets[0],deltas[0])
    parameters = wave['note'].decode('utf-8').split('\r')
    for parameter in parameters:
        parameterset = parameter.split('=')
        if len(parameterset) == 2:
            try:
                p[parameterset[0]] = float(parameterset[1])
            except(ValueError):
                p[parameterset[0]] = parameterset[1] #no way to convert string parameters
                    

    tstamp_s = p['Date'] + ' ' + p['Time']
    p['tstamp'] = time.mktime(datetime.datetime.strptime(tstamp_s,"%Y-%m-%d %H:%M:%S").timetuple())
    
    p['date'] = time.ctime(p['tstamp'])
    p['loadpath'] = path
    p['loadname'] = path[path.rfind('/')+1:]
    p['type'] = 'angledata'

    check_keys(p)

    return p            


######################### 2D loaders ##########################################
    
def loadData2DAPE(Aa2D,path, pin = std_p()):
    p = pin.copy()

    loaddict = igor.binarywave.load(path)
    #loadpxp: 
    #import igor.packed
    #igor.packed.load
    wave = loaddict['wave']
    p['version'] = loaddict['version'] #save version number
    
    #incorporate loading different regions from the same file
    rows,cols = wave['wave_header']['nDim'][0:2]

    Aa2D.reshape((rows,cols))
    Aa2D.data = np.array(wave['wData'], dtype = 'float32')
    deltas,offsets = wave['wave_header']['sfA'],wave['wave_header']['sfB']
    Aa2D.setXScale(offsets[0],deltas[0])
    Aa2D.setYScale(offsets[1],deltas[1])
    parameters = wave['note'].decode('utf-8').split('\r')
    for parameter in parameters:
        parameterset = parameter.split('=')
        if len(parameterset) == 2:
            try:
                p[parameterset[0]] = float(parameterset[1])
            except(ValueError):
                p[parameterset[0]] = parameterset[1] #no way to convert string parameters
                    

    tstamp_s = p['Date'] + ' ' + p['Time']
    p['tstamp'] = time.mktime(datetime.datetime.strptime(tstamp_s,"%Y-%m-%d %H:%M:%S").timetuple())
    p['type'] = 'angledata'  

    p['date'] = time.ctime(p['tstamp'])
    p['loadpath'] = path
    p['loadname'] = path[path.rfind('/')+1:]

    check_keys(p)
    return p
    
        
def loadData2DSLS(Aa2D, path, pin = std_p()):
    p = pin.copy()

    
    with h5py.File(path,'r') as hf:
        #data
        Aa2D.reshape(hf['Electron Analyzer/Image Data'].shape)
        Aa2D.data = hf['Electron Analyzer/Image Data'].value
        Aa2D.datalabel = hf['Electron Analyzer/Image Data'].attrs['Intensity Units'].decode()
        Aa2D.setXScale(*hf['Electron Analyzer/Image Data'].attrs['Axis0.Scale'])
        Aa2D.setYScale(*hf['Electron Analyzer/Image Data'].attrs['Axis1.Scale'])
        Aa2D.xlabel = hf['Electron Analyzer/Image Data'].attrs['Axis0.Description'].decode()
        Aa2D.ylabel = hf['Electron Analyzer/Image Data'].attrs['Axis1.Description'].decode()

        #extract image attrs
        for k,v in hf['Electron Analyzer/Image Data'].attrs.items():
            if type(v) == bytes:
                p[k] = v.decode()
            else:
                p[k] = v

        for k in hf['Other Instruments'].keys():
            if hf['Other Instruments'][k].shape[0] <= 1:
                p[k] = hf['Other Instruments'][k].value[0]
            else:
                p[k] = hf['Other Instruments'][k].value
        p['hv mode'] = hf['Other Instruments']['hv'].attrs['Mode'].decode()


    tstamp_s = p['Date Created'] + ' ' + p['Time Created']
    p['tstamp'] = time.mktime(datetime.datetime.strptime(tstamp_s,"%m/%d/%Y %I:%M:%S %p").timetuple())            
    p['date'] = time.ctime(p['tstamp'])
    p['type'] = 'angledata'
    p['Ep'] = p['Pass Energy (eV)']
    p['lens_type'] = p['Lens Mode']
    p['T_sample'] = p['Temperature B']

    p['date'] = time.ctime(p['tstamp'])
    p['loadpath'] = path
    p['loadname'] = path[path.rfind('/')+1:]

    #add theta angle to yscale:
    Aa2D.setYScale(p['Theta'] + Aa2D.ystart, Aa2D.ydelta)

    check_keys(p)

    return p
 

    
def loadData2DUBC(Aa2D, path, mode = 'avg', pin = std_p()):
    p = pin.copy()

    if(os.path.isdir(path)):
        dataext = '.detector_data'
        datafiles = [f for f in os.listdir(path) if (os.path.isfile(os.path.join(path, f)) and (f[-14:] == dataext))]
        parext = '.parms'
        parmfiles = [f for f in os.listdir(path) if (os.path.isfile(os.path.join(path, f)) and (f[-6:] == parext))]
    else:
        if(path.find('.detector_data') != -1):
            #if filename given with extension, strip and make new for parms
            datafiles = [path]
            parmfiles = [path[0:-14] + '.parms']
            path = ''
        else:
            #add extension
            datafiles = [path + '.detector_data']
            parmfiles = [path + '.parms']
            path = ''
            
    if(mode == 'avg'):
        pass
    elif(mode == 'first'):
        datafiles = [datafiles[0]]
        parmfiles = [parmfiles[0]]
    elif(mode[0] == 'f'): #use as f0, f1, f2 etc
        if(int(mode[1:]) < len(datafiles)):
            datafiles = [datafiles[int(mode[1:])]]
            parmfiles = [parmfiles[int(mode[1:])]]
        else:
            raise RuntimeError('Not enough data files')
    else:
        raise ValueError('Invalid mode')
        
    numfiles = float(len(datafiles))
 
   
    #load data
    for ii in range(int(numfiles)):        
        with open(os.path.join(path,datafiles[ii]), 'rb') as f:
            bytestring = f.read()
        
        size = 1376*1024 #for UBC data
        
        Aa2D.data += np.array(struct.unpack('%sf' % size,bytestring), dtype = 'float32').reshape((1024,1376)).T/numfiles


        with open(os.path.join(path,parmfiles[ii]), 'r') as f:
            filestring = f.read()
        
        parameters = filestring.split('\n')
    
    
        for parameter in parameters:
            parameterset = parameter.split(' = ')
            if len(parameterset) == 2:
                try:
                    value = float(parameterset[1])
                    try:
                        p[parameterset[0]] += value/numfiles
                    except KeyError:
                        p[parameterset[0]] = value/numfiles #create the key if it doesn't exist                         
                except(ValueError):
                    p[parameterset[0]] = parameterset[1] #no way to average string parameters
    p['eVpxEp'] = p['eV/px/Ep']
    del(p['eV/px/Ep'])
    p['T_sample'] = p['T_bullet']


    p['date'] = time.ctime(p['tstamp'])
    p['loadpath'] = path
    p['loadname'] = path[path.rfind('/')+1:]

    if(p['loadname'].find('spol') >= 0):
        p['pol'] = 'spol'
        p['gamma'] = 0.
    elif(p['loadname'].find('ppol') >= 0):
        p['pol'] = 'ppol'
        p['gamma'] = -90.
    elif(p['loadname'].find('_pol_') >= 0):
        p['pol'] = 'mix'
        p['gamma'] = float(p['loadname'][p['loadname'].find('_pol_')+5:p['loadname'].find('_pol_')+8])
    else:
        p['pol'] = 'unknown' 
        p['gamma'] = None

    p['type'] = 'rawdata'
    #set energy scaling
    xdelta = p['Ep']*p['eVpxEp']
    xstart = p['Ek'] - xdelta*(Aa2D.data.shape[0]-1)/2. 
    Aa2D.setXScale(xstart, xdelta)   

    check_keys(p)
    return p 

    
def loadData2DIgor(Aa2D,path, mode = 'avg', pin = std_p()):
    p = pin.copy()

    loaddict = igor.binarywave.load(path)
    #loadpxp: 
    #import igor.packed
    #igor.packed.load
    wave = loaddict['wave']
    p['version'] = loaddict['version'] #save version number
    
    #incorporate loading different regions from the same file
    rows,cols = wave['wave_header']['nDim'][0:2]

    Aa2D.reshape((rows,cols))
    Aa2D.data = wave['wData']
    deltas,offsets = wave['wave_header']['sfA'],wave['wave_header']['sfB']
    Aa2D.setXScale(offsets[0],deltas[0])
    Aa2D.setYScale(offsets[1],deltas[1])
    parameters = wave['note'].decode('utf-8').split('\r')
    for parameter in parameters:
        parameterset = parameter.split('=')
        if len(parameterset) == 2:
            try:
                p[parameterset[0]] = float(parameterset[1])
            except(ValueError):
                p[parameterset[0]] = parameterset[1] #no way to average string parameters
                    

#    tstamp_s = p['Date'] + ' ' + p['Time']
#    p['tstamp'] = time.mktime(datetime.datetime.strptime(tstamp_s,"%m/%d/%Y %I:%M:%S %p").timetuple())

#    p['date'] = time.ctime(p['tstamp'])
    p['loadpath'] = path
    p['loadname'] = path[path.rfind('/')+1:]

    p['type'] = 'angledata'
#    p['Ek'] = p['Center Energy']
#    p['Ep'] = p['Pass Energy']
#    p['lens_type'] = p['Lens Mode']

    check_keys(p)            
    return p
    

    
def loadData2DALS(Aa2D,path, mode = 'avg', pin = std_p()):
    p = pin.copy()
    #avg is not incorporated yet
    #load data dat has been converted form pxt to IBW
    loaddict = igor.binarywave.load(path)
    #loadpxp: 
    #import igor.packed
    #igor.packed.load
    wave = loaddict['wave']
    p['version'] = loaddict['version'] #save version number
    
    #incorporate loading different regions from the same file
    rows,cols = wave['wave_header']['nDim'][0:2]

    Aa2D.reshape((rows,cols))
    Aa2D.data = wave['wData']
    deltas,offsets = wave['wave_header']['sfA'],wave['wave_header']['sfB']
    Aa2D.setXScale(offsets[0],deltas[0])
    Aa2D.setYScale(offsets[1],deltas[1])
    parameters = wave['note'].decode('utf-8').replace('\r','\n').split('\n')
    for parameter in parameters:
        parameterset = parameter.split('=')
        if len(parameterset) == 2:
            try:
                p[parameterset[0]] = float(parameterset[1])
            except(ValueError):
                p[parameterset[0]] = parameterset[1]

    try:
        if p['Energy Scale'] == 'Binding':
            Aa2D.setXScale(-Aa2D.xstart,-Aa2D.xdelta)
                        
    
        tstamp_s = p['Date'] + ' ' + p['Time']
        p['tstamp'] = time.mktime(datetime.datetime.strptime(tstamp_s,"%m/%d/%Y %I:%M:%S %p").timetuple())
        p['date'] = time.ctime(p['tstamp'])
        p['type'] = 'angledata'  
    
        p['date'] = time.ctime(p['tstamp'])
        p['loadpath'] = path
        p['loadname'] = path[path.rfind('/')+1:]
    
    
        p['Ek'] = p['Center Energy']
        p['Ep'] = p['Pass Energy']
        p['lens_type'] = p['Lens Mode']
        p['slit_entrance'] = p['Slit Plate']
    except KeyError:
        pass

    p['type'] = 'angledata'

    check_keys(p)  
    return p




def loadData2DCLS(Aa2D, path = '', pin = std_p()): 
    """Args:
        path: should be the full path to datafile, without underscore
        cycles: can be a number (load only that cycle) or 'avg', which averages, or a range"""
    
    
    paramfile = path + '_param.txt'
    

    data_slice = np.loadtxt(path)
    Aa2D.reshape(data_slice.shape)
    Aa2D.data = data_slice
    
    with open(paramfile,'r') as f:
        flines = f.read().split('\n')
        
    angle_off, angle_delta = [float(s.split()[2]) for s in (flines.pop(-1)).split(';')]
    energy_off, energy_delta = [float(s.split()[2]) for s in (flines.pop(-1)).split(';')]
    
    Aa2D.setXScale(energy_off, energy_delta)
    Aa2D.setYScale(angle_off, angle_delta)
    
    p = pin.copy()

    #load parameters

    for fline in flines:
        key, v = fline.split(': ')
        try:
            v = float(v) if '.' in v else int(v)
        except ValueError:
            v = v
        p[key] = v

    
    print('No time in this datatype, using creation time of parameter file')
    
    p['tstamp'] = os.path.getmtime(paramfile)        
    p['date'] = time.ctime(p['tstamp'])
    
    
    
    p['type'] = 'rawdata'
    p['Ep'] = p['Pass energy']
    p['lens_type'] = p['Lens mode']

    p['loadpath'] = path
    p['loadname'] = path[path.rfind('/')+1:]


    return p
    



    
def loadH5_2D(Aa2D, path, objtype = 'AaData2D'):     
    p = {}        
    with h5py.File(path,'r') as hf:
        #header
        if(hf['h/ver'].value.decode() != '1.0'):
            raise RuntimeError('This loader is only for version 1.0 data')
        if(hf['h/obj'].value.decode() != objtype):
            raise RuntimeError('Wrong obj type, trying to load {} data with in a {} container'.format(hf['h/obj'].value.decode(), objtype))                
        #data
        Aa2D.reshape(hf['d/data'].shape)
        Aa2D.data = hf['d/data'].value
        Aa2D.datalabel = hf['d/datalabel'].value.decode()
        Aa2D.setXScale(hf['d/xstart'].value, hf['d/xdelta'].value)
        Aa2D.setYScale(hf['d/ystart'].value, hf['d/ydelta'].value)
        Aa2D.xlabel = hf['d/xlabel'].value.decode()
        Aa2D.ylabel = hf['d/ylabel'].value.decode()
        Aa2D.name = hf['d/name'].value.decode()
        Aa2D.note = hf['d/note'].value.decode()
        #parameters:
        if 'p/str' in hf:
            for k in hf['p/str']:#str
                val = hf['p/str'][k].value.decode()
                p[k] = val
        if 'p/none' in hf:
            for k in hf['p/none']: #none
                p[k] = None
        if 'p/bool' in hf:                    
            for k in hf['p/bool']:
                if hf['p/bool'][k].shape[0] > 1:
                    p[k] = hf['p/bool'][k].value
                else:
                    p[k] = hf['p/bool'][k].value[0]
        if 'p/num' in hf:
            for k in hf['p/num']:#num
                if hf['p/num'][k].shape[0] > 1:
                    p[k] = hf['p/num'][k].value
                else:
                    p[k] = hf['p/num'][k].value[0]        
    check_keys(p)
    return p    
    

######################### 3D loaders ##########################################
        
def loadData3DSLS(Aa3D, path = '', pin = std_p()):    
    p = pin.copy()
    
    with h5py.File(path,'r') as hf:
        #data
        Aa3D.reshape(hf['Electron Analyzer/Image Data'].shape)
        Aa3D.data = hf['Electron Analyzer/Image Data'].value
        Aa3D.datalabel = hf['Electron Analyzer/Image Data'].attrs['Intensity Units'].decode()
        Aa3D.setXScale(*hf['Electron Analyzer/Image Data'].attrs['Axis0.Scale'])
        Aa3D.setYScale(*hf['Electron Analyzer/Image Data'].attrs['Axis1.Scale'])
        Aa3D.setZScale(*hf['Electron Analyzer/Image Data'].attrs['Axis2.Scale'])
        Aa3D.xlabel = hf['Electron Analyzer/Image Data'].attrs['Axis0.Description'].decode()
        Aa3D.ylabel = hf['Electron Analyzer/Image Data'].attrs['Axis1.Description'].decode()
        Aa3D.zlabel = hf['Electron Analyzer/Image Data'].attrs['Axis2.Description'].decode()

        #extract image attrs
        for k,v in hf['Electron Analyzer/Image Data'].attrs.items():
            if type(v) == bytes:
                p[k] = v.decode()
            else:
                p[k] = v

        for k in hf['Other Instruments'].keys():
            if hf['Other Instruments'][k].shape[0] <= 1:
                p[k] = hf['Other Instruments'][k].value[0]
            else:
                p[k] = hf['Other Instruments'][k].value
        p['hv mode'] = hf['Other Instruments']['hv'].attrs['Mode'].decode()


    tstamp_s = p['Date Created'] + ' ' + p['Time Created']
    p['tstamp'] = time.mktime(datetime.datetime.strptime(tstamp_s,"%m/%d/%Y %I:%M:%S %p").timetuple())            
    p['date'] = time.ctime(p['tstamp'])
    p['type'] = 'angledata'
    p['Ep'] = p['Pass Energy (eV)']
    p['lens_type'] = p['Lens Mode']

    p['loadpath'] = path
    p['loadname'] = path[path.rfind('/')+1:]


    #add theta angle to yscale:
    Aa3D.setYScale(p['Theta'] + Aa3D.ystart, Aa3D.ydelta)

    return p
    
def loadData3DALS(Aa3D,path,zaxis = 'polar', pin = std_p()):
    p = pin.copy()
    #path is basename, give as 'path/sri_024'
    
    #get files:
    files = []
    i = 1
    while True:
        if os.path.exists(path + '_S%03d.ibw' % i):
            files.append(path + '_S%03d.ibw' % i)
        else:
            break
        i += 1
    
    nfiles = len(files)
    if nfiles == 0:
        raise FileNotFoundError('No matching files found')

    for i in range(nfiles):   
        
        loaddict = igor.binarywave.load(files[i])
        #loadpxp: 
        #import igor.packed
        #igor.packed.load
        wave = loaddict['wave']
        parameters = wave['note'].decode('utf-8').replace('\r','\n').split('\n')
        
        if i == 0:
            rows,cols = wave['wave_header']['nDim'][0:2]
            Aa3D.reshape((rows,cols,nfiles))
            deltas,offsets = wave['wave_header']['sfA'],wave['wave_header']['sfB']
            Aa3D.setXScale(offsets[0],deltas[0])
            Aa3D.setYScale(offsets[1],deltas[1])
        
            for parameter in parameters:
                parameterset = parameter.split('=')
                if len(parameterset) == 2:
                    try:
                        float(parameterset[1])
                        p[parameterset[0]] = np.zeros((nfiles), dtype = 'float32')
                    except(ValueError):
                        p[parameterset[0]] = np.zeros((nfiles),dtype = '<U50')
            p['tstamp'] = np.zeros((nfiles),dtype = 'float64')

        

        for parameter in parameters:
            parameterset = parameter.split('=')
            if len(parameterset) == 2:
                try:
                    p[parameterset[0]][i] = float(parameterset[1])
                except(ValueError):
                    p[parameterset[0]][i] = parameterset[1]
        tstamp_s = p['Date'][i] + ' ' + p['Time'][i]
        p['tstamp'][i] = time.mktime(datetime.datetime.strptime(tstamp_s,"%m/%d/%Y %I:%M:%S %p").timetuple())
        Aa3D.data[:,:,i] = wave['wData']


    for k in p.keys():
        if not((k == 'tstamp') or (k == 'Azimuth') or (k == 'BL Energy') or (k == 'Beam Current') or (k == 'Entrance Slit') \
        or (k == 'Exit Slit') or (k == 'Polar') or (k == 'Pressure') or (k == 'Temperature Sensor A') \
        or (k == 'Temperature Sensor B') or (k == 'Temperature Sensor C') or (k == 'Temperature Sensor D') or (k == 'Tilt') \
        or (k == 'Sample X') or (k == 'Sample Y') or (k == 'Sample Z') or (k == 'Mesh Current')):
            if type(p[k]) == np.ndarray:
                p[k] = p[k][0]


    p['version'] = loaddict['version'] #save version number
    
    if zaxis == 'polar':
        Aa3D.setZScale(p['Polar'][0],p['Polar'][1]-p['Polar'][0])
    elif zaxis == 'hv':
        Aa3D.setZScale(p['BL Energy'][0],p['BL Energy'][1]-p['BL Energy'][0])   
    elif zaxis == 'T':
        Aa3D.setZScale(p['Temperature Sensor B'][0],p['Temperature Sensor B'][1]-p['Temperature Sensor B'][0])   
        Aa3D.zaxis = np.copy(p['Temperature Sensor B'])   

    if p['Energy Scale'] == 'Binding':
        Aa3D.setXScale(-Aa3D.xstart,-Aa3D.xdelta)
    p['date'] = time.ctime(p['tstamp'][0])
    p['type'] = 'angledata'  
    p['loadpath'] = path
    p['loadname'] = path[path.rfind('/')+1:]
              
    return p
    
    
    
def loadData3DAPE_hv(Aa3D,paths, pin = std_p()):
    p = pin.copy()
    #path is basename, give as 'path/sri_024'
    
    #get files:
    files = []
    
    if type(paths) == str:
        ##not supported yet:
        raise NotImplementedError('Not supported yet')
        i = 1
        while True:
            if os.path.exists(paths + '_S%03d.ibw' % i):
                files.append(paths + '_S%03d.ibw' % i)
            else:
                break
            i += 1
    elif type(paths) == list:
        files = paths
    else:
        raise TypeError('path needs to be basepath or list of paths')
    
    nfiles = len(files)
    if nfiles == 0:
        raise FileNotFoundError('no files found')
    for i in range(nfiles):   
        
        loaddict = igor.binarywave.load(files[i])
        #loadpxp: 
        #import igor.packed
        #igor.packed.load
        wave = loaddict['wave']
        parameters = wave['note'].decode('utf-8').replace('\r','\n').split('\n')
        
        if i == 0:
            rows,cols = wave['wave_header']['nDim'][0:2]
            Aa3D.reshape((rows,cols,nfiles))
            deltas,offsets = wave['wave_header']['sfA'],wave['wave_header']['sfB']
            Aa3D.setXScale(offsets[0],deltas[0])
            Aa3D.setYScale(offsets[1],deltas[1])
        
            for parameter in parameters:
                parameterset = parameter.split('=')
                if len(parameterset) == 2:
                    try:
                        float(parameterset[1])
                        p[parameterset[0]] = np.zeros((nfiles), dtype = 'float32')
                    except(ValueError):
                        p[parameterset[0]] = np.zeros((nfiles),dtype = '<U50')
            p['tstamp'] = np.zeros((nfiles),dtype = 'float64') 

        

        for parameter in parameters:
            parameterset = parameter.split('=')
            if len(parameterset) == 2:
                try:
                    p[parameterset[0]][i] = float(parameterset[1])
                except(ValueError):
                    p[parameterset[0]][i] = parameterset[1]
        tstamp_s = p['Date'][i] + ' ' + p['Time'][i]
        p['tstamp'][i] = time.mktime(datetime.datetime.strptime(tstamp_s,"%Y-%m-%d %H:%M:%S").timetuple())

        Aa3D.data[:,:,i] = wave['wData']


    for k in p.keys():
        if not((k == 'tstamp') or (k == 'Excitation Energy') or (k == 'Low Energy')\
               or (k == 'Center Energy') or (k == 'High Energy')):
            try:
                p[k] = p[k][0]
            except (TypeError,IndexError):
                pass


    p['version'] = loaddict['version'] #save version number
    
    Aa3D.setZScale(p['Excitation Energy'][0],p['Excitation Energy'][1]-p['Excitation Energy'][0])            


    p['date'] = time.ctime(p['tstamp'][0])
    p['type'] = 'angledata'  
    p['loadpath'] = files[0]
    p['loadname'] = files[0][files[0].rfind('/')+1:]
              
    return p
 
    

def loadData3DAPE_FS(Aa3D, path = '', pin = std_p()):    
    """Load FS for APE data, path should refer to the unzipped folder with bin,
    and ini files"""
    p = pin.copy()
    
    files = os.listdir(path)
    binfile = [name for name in files if '.bin' in name][0] #there should only be one bin per FS, except if you do multiple regions, in that case, separate into folder manually
    specini = binfile[:binfile.find('.bin')]+'.ini' #second file we need is the Spectrum_name
    nameini = specini[len('Spectrum_'):] # last one we need is name.ini, viewers are only viewer settings
    
    #load paramters
    with open(path + '/' + specini,'r') as f:
        parameters = f.read().split('\n') #get separate lines
    with open(path + '/' + nameini, 'r') as f:
        parameters += f.read().split('\n')
        
    for parameter in parameters:
        parameterset = parameter.split('=')
        if len(parameterset) == 2:
            try:
                p[parameterset[0]] = float(parameterset[1])
            except(ValueError):
                p[parameterset[0]] = parameterset[1]    
    ## first find ini files to load meta data, get parameters
    
    rows,cols,lays = int(p['width']),int(p['height']),int(p['depth'])

    Aa3D.reshape((rows,cols,lays))
    Aa3D.setXScale(p['widthoffset'],p['widthdelta'])
    Aa3D.setYScale(p['heightoffset'],p['heightdelta'])
    Aa3D.setZScale(p['depthoffset'],p['depthdelta'])
    Aa3D.xlabel,Aa3D.ylabel,Aa3D.zlabel = p['widthlabel'],p['heightlabel'],p['depthlabel']
    
    Aa3D.data = np.fromfile(path + '/' + binfile,dtype = 'float32').reshape((lays,cols,rows)).transpose((2,1,0))

    
    tstamp_s = p['Date'] + ' ' + p['Time']
    p['tstamp'] = time.mktime(datetime.datetime.strptime(tstamp_s,"%Y-%m-%d %H:%M:%S").timetuple())
    p['date'] = time.ctime(p['tstamp'])
    p['type'] = 'angledata'
    p['Ep'] = p['Pass Energy']
    p['lens_type'] = p['Lens Mode']

    p['loadpath'] = path
    p['loadname'] = path[path.rfind('/')+1:]

    return p    
    
 

   
def loadData3DUBC(Aa3D, path = '', pin = std_p()):
    p = pin.copy()

        
    #list .detector_data files
    dataext = '.detector_data'
    datafiles = [f for f in os.listdir(path) if (os.path.isfile(os.path.join(path, f)) and (f[-14:] == dataext))]
    parext = '.parms'
    parmfiles = [f for f in os.listdir(path) if (os.path.isfile(os.path.join(path, f)) and (f[-6:] == parext))]
    
    if not(len(parmfiles) == len(datafiles)):
        print('Number of parameter files does not match number of data files')
        
    numfiles = len(datafiles)
    
    
    Aa3D.reshape((1376,1024,numfiles))
    
    for kk in range(numfiles):
        #load data
        with open(os.path.join(path,datafiles[kk]), 'rb') as f:
            bytestring = f.read()
        
        size = 1376*1024 #for UBC data
        
        Aa3D.data[:,:,kk] = np.array(struct.unpack('%sf' % size,bytestring), dtype = 'float32').reshape((1024,1376)).T
        

        with open(os.path.join(path,parmfiles[kk]), 'r') as f:
            filestring = f.read()
        
        parameters = filestring.split('\n')
        if(kk == 0):        
            for parameter in parameters:
                parameterset = parameter.split(' = ')
                if len(parameterset) == 2:
                    try:
                        value = float(parameterset[1])
                        p[parameterset[0]] = np.zeros((numfiles),dtype = 'float32') #make array
                        p[parameterset[0]][kk] = value
                    except(ValueError):
                        p[parameterset[0]] = np.zeros((numfiles),dtype = '<U50') #make array
                        p[parameterset[0]][kk] = parameterset[1]
        else:
            for parameter in parameters:
                parameterset = parameter.split(' = ')
                if len(parameterset) == 2:
                    try:
                        value = float(parameterset[1])
                        p[parameterset[0]][kk] = value
                    except(ValueError):
                        p[parameterset[0]][kk] = parameterset[1]                
    for k in p.keys():
        if not((k == 'tstamp') or (k == 'cryo_phi') or (k == 'T_bullet') or (k == 'P_lower_chamber') or (k == 'cryo_z') \
        or (k == 'T_cold_shield') or (k == 'P_mono') or (k == 'cryo_theta') or (k == 'cryo_y') \
        or (k == 'P_upper_chamber') or (k == 'cryo_x') or (k == 'T_reservoir') or (k == 'P_lamp')):
            p[k] = p[k][0]
    p['eVpxEp'] = p['eV/px/Ep']
    del(p['eV/px/Ep'])
    p['type'] = 'rawdata'
    #set energy scaling
    xdelta = p['Ep']*p['eVpxEp']
    xstart = p['Ek'] - xdelta*(Aa3D.data.shape[0]-1)/2. 
    Aa3D.setXScale(xstart, xdelta)  
    Aa3D.setZScale(p['cryo_phi'][0],p['cryo_phi'][1]-p['cryo_phi'][0])
    p['date'] = time.ctime(p['tstamp'][0])    

    p['loadpath'] = path
    p['loadname'] = path[path.rfind('/')+1:]

    if(p['loadname'].find('spol') >= 0):
        p['pol'] = 'spol'
        p['gamma'] = 0.
    elif(p['loadname'].find('ppol') >= 0):
        p['pol'] = 'ppol'
        p['gamma'] = -90.
    elif(p['loadname'].find('_pol_') >= 0):
        p['pol'] = 'mix'
        p['gamma'] = float(p['loadname'][p['loadname'].find('_pol_')+5:p['loadname'].find('_pol_')+8])
    else:
        p['pol'] = 'unknown' 
        p['gamma'] = None


    
    return p



def loadData3DCLS(Aa3D, path = '', pin = std_p(), cycles = 'avg'): 
    """Args:
        path: should be the full path to the base filename without the underscore
        cycles: can be a number (load only that cycle) or 'avg', which averages, or a range"""
    
    folder,name = os.path.split(path)

    
    files = [folder + '/' + file for file in os.listdir(folder) if (name == file[0:len(name)])]
    
    paramfile = files.pop(files.index(folder + '/' + name + '_param.txt'))
    
    #get max cycles and steps:
    n_cycles = -1
    n_steps = -1
    for f in files:
        splits = f.replace('.','_').split('_')
        cycle = int(splits[splits.index('Cycle')+1])
        if cycle > n_cycles:
            n_cycles = cycle
        step = int(splits[splits.index('Step')+1])
        if step > n_steps:
            n_steps = step
    n_cycles += 1
    n_steps += 1    
    if n_cycles == 0:
        raise RuntimeError('No cycles found with path')
    if n_steps == 0:
        raise RuntimeError('No steps found with path')
        
    if cycles == 'avg':
        cycles = range(n_cycles)
    elif type(cycles) == int:
        cycles = [cycles]

    for i,c in enumerate(cycles):
        for s in range(n_steps):
            print('Now loading cycle {} and step {}'.format(c,s))
            file = folder + '/' + name + '_Cycle_{}_Step_{}.txt'.format(c,s)
            data_slice = np.loadtxt(file)
            if i == 0 and s == 0:
                shape = data_slice.shape + (n_steps,)
                Aa3D.reshape(shape)
            Aa3D.data[:,:,s] = data_slice
    Aa3D.data /= len(cycles)
    
    
    with open(paramfile,'r') as f:
        flines = f.read().split('\n')
    
    zscale = np.fromstring(flines.pop(-1), sep = ' ')
    zdiff = zscale[1:] - zscale[:-1]
    if not np.all(zdiff[0] == zdiff):
        raise RuntimeError('Code is not setup for non-uniform scaling')
    Aa3D.setZScale(zscale[0], zdiff[0])
    
    angle_off, angle_delta = [float(s.split()[2]) for s in (flines.pop(-1)).split(';')]
    energy_off, energy_delta = [float(s.split()[2]) for s in (flines.pop(-1)).split(';')]
    
    Aa3D.setXScale(energy_off, energy_delta)
    Aa3D.setYScale(angle_off, angle_delta)
    Aa3D.setZScale(zscale[0], zdiff[0])
    
    p = pin.copy()

    #load parameters

    for fline in flines:
        key, v = fline.split(': ')
        try:
            v = float(v) if '.' in v else int(v)
        except ValueError:
            v = v
        p[key] = v

    
    print('No time in this datatype, using creation time of parameter file')
    
    p['tstamp'] = os.path.getmtime(paramfile)        
    p['date'] = time.ctime(p['tstamp'])
    
    
    
    p['type'] = 'rawdata'
    p['Ep'] = p['Pass energy']
    p['lens_type'] = p['Lens mode']

    p['loadpath'] = path
    p['loadname'] = path[path.rfind('/')+1:]


    return p
    

    
def loadMaestro_3D(Aa3D, path):        
    p = std_p()        
    with h5py.File(path,'r') as hf:
        #0D data:
        for key in hf['0D_Data'].keys():
            p['0D_' + key] = np.array(hf['0D_Data'][key])

        #data
        inshape = hf['2D_Data']['Fixed_Spectra0'].shape
        Aa3D.reshape((inshape[1],inshape[0],inshape[2]))
        Aa3D.data = np.array(hf['2D_Data']['Fixed_Spectra0']).transpose((1,0,2))
        
        offset = hf['2D_Data']['Fixed_Spectra0'].attrs['scaleOffset']
        delta = hf['2D_Data']['Fixed_Spectra0'].attrs['scaleDelta']
        units = hf['2D_Data']['Fixed_Spectra0'].attrs['unitNames']
        
        
        angle_delta = 0.045
        
        offset_angle = (inshape[0] * angle_delta)/2
        
        Aa3D.setXScale(offset[0], delta[0])
        Aa3D.setYScale(offset_angle, angle_delta)
        Aa3D.xlabel =units[0]
        Aa3D.ylabel = units[1]
        Aa3D.datalabel = hf['2D_Data']['Fixed_Spectra0'].attrs['dataUnitName']
        
        
        
        
        for k,v in zip(hf['Comments']['PreScan'].dtype.names, hf['Comments']['PreScan'][0]):
            p[k] = v.decode()
            
            
        for key in hf['Headers'].keys():
            dset = hf['Headers'][key]
            for i in range(len(dset)):
                keyname = dset[i]['longname'].decode() + ' ' + dset[i]['comment'].decode()
                value = dset[i]['value'].decode()
                try:
                    p[keyname] = float(value)
                except ValueError:
                    p[keyname] = value
    p['type'] = 'angledata'
                    

    return p

    
def loadH5_3D(Aa3D, path):        
    p = {}        
    with h5py.File(path,'r') as hf:
        #header
        if(hf['h/ver'].value.decode() != '1.0'):
            raise RuntimeError('This loader is only for version 1.0 data')
        if(hf['h/obj'].value.decode() != 'AaData3D'):
            raise RuntimeError('Wrong obj type, trying to load {} data with in a 3D container'.format(hf['h/obj'].value.decode()))                
        
        #data
        Aa3D.reshape(hf['d/data'].shape)
        Aa3D.data = hf['d/data'].value
        Aa3D.datalabel = hf['d/datalabel'].value.decode()
        Aa3D.setXScale(hf['d/xstart'].value, hf['d/xdelta'].value)
        Aa3D.setYScale(hf['d/ystart'].value, hf['d/ydelta'].value)
        Aa3D.setZScale(hf['d/zstart'].value, hf['d/zdelta'].value)
        Aa3D.xlabel = hf['d/xlabel'].value.decode()
        Aa3D.ylabel = hf['d/ylabel'].value.decode()
        Aa3D.zlabel = hf['d/zlabel'].value.decode()
        Aa3D.name = hf['d/name'].value.decode()
        Aa3D.note = hf['d/note'].value.decode()
        #parameters:
        if 'p/str' in hf:
            for k in hf['p/str']:#str
                val = hf['p/str'][k].value.decode()
                p[k] = val
        if 'p/none' in hf:
            for k in hf['p/none']: #none
                p[k] = None
        if 'p/bool' in hf:                    
            for k in hf['p/bool']:
                if hf['p/bool'][k].shape[0] > 1:
                    p[k] = hf['p/bool'][k].value
                else:
                    p[k] = hf['p/bool'][k].value[0]
        if 'p/num' in hf:
            for k in hf['p/num']:#num
                if hf['p/num'][k].shape[0] > 1:
                    p[k] = hf['p/num'][k].value
                else:
                    p[k] = hf['p/num'][k].value[0]
    return p
 

       

