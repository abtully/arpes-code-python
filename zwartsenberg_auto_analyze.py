# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 17:02:59 2016

@author: berend
"""

import os
from arpes.data import UBCSlits, UBCGold, AaData2D, AaData3D
import numpy as np
import time



def loadc(path, subpath, slits, gold, objtype = 'Data3D', savepath = None, angles = None, sample = None, compound = None):
    # find full path
    if os.path.isdir(path+'/'+subpath):
        loadpath = path+'/'+subpath
    else:
        folderlist = [d for d in os.listdir(path) if (os.path.isdir(os.path.join(path, d)) and d.find(subpath)>=0)]
        if len(folderlist) == 1:
            loadpath = path+'/'+folderlist[0]
        else:
            raise FileNotFoundError('path not found')

    funcdict = {'slits' : UBCSlits, 'gold' : UBCGold, 'Data2D' : AaData2D, 'Data3D' : AaData3D}            
            
    data = funcdict[objtype](loadpath)
    data.Correct(slits,gold)
    
    if angles:
        data.setangles(*angles)
    if sample:
        data.sample = sample
    if compound:
        data.compound = compound
    
    data.name = loadpath[loadpath.rfind('/')+1:loadpath.rfind('/')+4]
    
    if savepath:
        data.save((savepath+'/'+data.name))
    
    return data


    
    
# save info should be (fullpath, Ep,Ek,lens_type,gam,tstamp)
class slitdb:
    def __init__(self, path):
        self.path = path
        self.dbfilepath = self.path + '/db'
        self.db = []
        if not os.path.exists(path):
            raise FileNotFoundError('db path not found')
        
        if os.path.exists(self.dbfilepath):
            with open(self.dbfilepath, 'r') as f:
                dbstrlist = f.read().split('\n')
            for dbstr in dbstrlist:
                dbitem = dbstr.split('\t')
                self.db.append((dbitem[0], float(dbitem[1]), float(dbitem[2]), dbitem[3], float(dbitem[4]), float(dbitem[5])))
            
        
    def add(self,slitobj):
        if not slitobj.__class__ == UBCSlits:
            raise ValueError('Only UBCSlits objects')
            
        if not slitobj.analyzed:
            raise RuntimeError('Analyze slits first')
            

        
        for i in range(len(self.db)):
            if np.isclose(self.db[i][5]-slitobj.tstamp,0.0): #if tstamp difference is close to zero
                raise RuntimeError('Entry already in db')

        savename = slitobj.ID
        fullpath = self.path + '/' + savename + '.pickle'
        slitobj.save(fullpath[:-7])
        self.db.append((fullpath,slitobj.Ep, slitobj.Ek, slitobj.lens_type,slitobj.gamma, slitobj.tstamp))
        self.writedb()
        
    def remove(self,item = None):
        if not item: #default is last added item
            item = len(self.db)-1
        os.remove(self.db[item][0])
        self.db.pop(item)
        self.writedb()
    
    def writedb(self):
        with open(self.dbfilepath, 'w') as f:
            f.write('\n'.join(['\t'.join([str(elem) for elem in dbitem]) for dbitem in self.db]))
            
    def loadslit(self, item = None):
        if not item:
            item = -1
        return UBCSlits(self.db[item][0])
    
    def determineslit(self, Aaobj): #return itemno 
        shortlist = []
        Ek = Aaobj.Ek
        Ep = Aaobj.Ep
        lens_type = Aaobj.lens_type
        gamma = Aaobj.gamma
        for i in range(len(self.db)):
            dbitem = self.db[i]
            if ((np.isclose(dbitem[1],Ep)) and (np.isclose(dbitem[2],Ek)) and (dbitem[3] == lens_type) and (np.isclose(dbitem[4],gamma))):
                shortlist.append(i)    
        if len(shortlist) == 0:
            raise RuntimeError('No suitable correction found')
            return
        best = shortlist[0]
        for i in shortlist:
            if abs(Aaobj.tstamp - self.db[i][5]) < abs(Aaobj.tstamp - self.db[best][5]):
                best = i
        return self.loadslit(best)

    def getmostrecent(self, Ek, Ep, lens_type, gamma): #return itemno 
        shortlist = []

        for i in range(len(self.db)):
            dbitem = self.db[i]
            if ((np.isclose(dbitem[1],Ep)) and (np.isclose(dbitem[2],Ek)) and (dbitem[3] == lens_type) and (np.isclose(dbitem[4],gamma))):
                shortlist.append(i)    
        if len(shortlist) == 0:
            raise RuntimeError('No suitable correction found')
            return
        best = shortlist[0]
        for i in shortlist:
            if abs(time.time() - self.db[i][5]) < abs(time.time() - self.db[best][5]):
                best = i
        return self.loadslit(best)
        
    def getiteminfo(self, item = None): #return last item as default
        if not item:
            item = -1
        return self.db[item]
            
            

class golddb:
    def __init__(self, path):
        self.path = path
        self.dbfilepath = self.path + '/db'
        self.db = []
        if not os.path.exists(path):
            raise RuntimeError('db path not found')
            return
        
        if os.path.exists(self.dbfilepath):
            with open(self.dbfilepath, 'r') as f:
                dbstrlist = f.read().split('\n')
            for dbstr in dbstrlist:
                dbitem = dbstr.split('\t')
                self.db.append((dbitem[0], float(dbitem[1]), float(dbitem[2]), dbitem[3], float(dbitem[4]), float(dbitem[5])))
            
        
    def add(self,goldobj):
        if not goldobj.__class__ == UBCGold:
            raise ValueError('Only UBCGold objects')
            
        if not goldobj.analyzed:
            raise RuntimeError('Analyze slits first')
            
        for i in range(len(self.db)):
            if np.isclose(self.db[i][5]-goldobj.tstamp,0.0): #if tstamp difference is close to zero
                raise RuntimeError('Entry already in db')
        
        savename = goldobj.ID

        fullpath = self.path + '/' + savename + '.pickle'
        goldobj.save(fullpath[:-7])
        self.db.append((fullpath,goldobj.Ep, goldobj.Ek, goldobj.lens_type,goldobj.gamma, goldobj.tstamp))
        self.writedb()
        
    def remove(self,item = None):
        if not item: #default is last added item
            item = len(self.db)-1
        os.remove(self.db[item][0])
        self.db.pop(item)
        self.writedb()
    
    def writedb(self):
        with open(self.dbfilepath, 'w') as f:
            f.write('\n'.join(['\t'.join([str(elem) for elem in dbitem]) for dbitem in self.db]))
            
    def loadgold(self, item = None):
        if not item:
            item = -1
        return UBCGold(self.db[item][0])
    
    def determinegold(self, Aaobj): #return itemno 
        shortlist = []
        Ek = Aaobj.Ek
        Ep = Aaobj.Ep
        lens_type = Aaobj.lens_type
        gamma = Aaobj.gamma
        for i in range(len(self.db)):
            dbitem = self.db[i]
            if ((np.isclose(dbitem[1],Ep)) and (np.isclose(dbitem[2],Ek)) and (dbitem[3] == lens_type) and (np.isclose(dbitem[4],gamma))):
                shortlist.append(i)    
        if len(shortlist) == 0:
            raise RuntimeError('No suitable correction found')
            return
        best = shortlist[0]
        for i in shortlist:
            if abs(Aaobj.tstamp - self.db[i][5]) < abs(Aaobj.tstamp - self.db[best][5]):
                best = i
        return self.loadgold(best)


    def getmostrecent(self, Ek, Ep, lens_type, gamma): #return itemno 
        shortlist = []

        for i in range(len(self.db)):
            dbitem = self.db[i]
            if ((np.isclose(dbitem[1],Ep)) and (np.isclose(dbitem[2],Ek)) and (dbitem[3] == lens_type) and (np.isclose(dbitem[4],gamma))):
                shortlist.append(i)    
        if len(shortlist) == 0:
            raise RuntimeError('No suitable correction found')
            return
        best = shortlist[0]
        for i in shortlist:
            if abs(time.time() - self.db[i][5]) < abs(time.time() - self.db[best][5]):
                best = i
        return self.loadgold(best)

    def getiteminfo(self, item = None): #return last item as default
        if not item:
            item = -1
        return self.db[item]



def getslitsandgold(path):
    print(path)
    
    #for now only recipe data
    folderlist = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    slitpathlist = []
    goldpathlist = []

    
    #sort out slits, gold and 2D data and 3D data
    for d in folderlist:
        if(d.find('slits') >= 0):
            slitpathlist.append(path + '/' + d)
        elif(d.find('gold') >= 0):
            goldpathlist.append(path + '/' + d)

    return slitpathlist, goldpathlist


def loadlist(pathlist, objtype = 'slits'):
    funcdict = {'slits' : UBCSlits, 'gold' : UBCGold, 'data2D' : AaData2D, 'data3D' : AaData3D}
    outdict = {}
    for path in pathlist:
        #get f name:
        name = path[path.rfind('/')+1:path.rfind('/')+4]
        print(name)
        outdict[name] = funcdict[objtype](path)
        
    return outdict
    
    
def getdatafiles(path):
    print(path)
    
    #for now only recipe data
    folderlist = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    pathlist = []

    
    #sort out slits, gold and 2D data and 3D data
    for d in folderlist:
        if((d.find('slits') == -1) and (d.find('gold') == -1)):
            pathlist.append(path + '/' + d)


    return pathlist
  

def addslitstodb(path, sdb):
    
    print(path)
    folderlist = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    slitpathlist = []
    for d in folderlist:
        if(d.lower().find('slits') >= 0):
            slitpathlist.append(d)

    print(slitpathlist)
    notanalyzed = []
            
    for slitpath in slitpathlist:
        slits = UBCSlits(path = (path + '/' + slitpath))
        print(slits.ID)
        print(slits.loadpath)
        try:
            slits.Analyze() 
            sdb.add(slits)
            print('Added %s' % slits.ID)
        except:
            print('failed to analyze %s' % slits.ID)
            notanalyzed.append(slits)
    print('Done with slits')
    return notanalyzed
            
def addgoldtodb(path, gdb, sdb):
    
    print(path)
    folderlist = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    goldpathlist = []
    for d in folderlist:
        if((d.lower().find('gold') >= 0) or (d.lower().find('au') >= 0)):
            goldpathlist.append(d)

    print(goldpathlist)
    notanalyzed = []
            
    for goldpath in goldpathlist:
        gold = UBCGold(path = (path + '/' + goldpath))
        print(gold.ID)
        print(gold.loadpath)

        try:
            slits = sdb.loadslit(sdb.determineslit(gold))
            gold.Correct(slits)
            gold.Analyze() 
            gdb.add(gold)
            print('Added %s' % gold.ID)
        except:
            print('failed to analyze %s' % gold.ID)
            notanalyzed.append(gold)
    print('Done with gold')
    return notanalyzed    

        
def AnalyzeFolder(path):
    
    print(path)
    
    #for now only recipe data
    folderlist = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    slitpathlist = []
    goldpathlist = []
    data2Dpathlist = []
    data3Dpathlist = []
    
    #sort out slits, gold and 2D data and 3D data
    for d in folderlist:
        if(d.find('slits') >= 0):
            slitpathlist.append(d)
        elif(d.find('gold') >= 0):
            goldpathlist.append(d)
        elif((d.find('2D') >= 0) or (d.find('slice') >= 0) or (d.find('cut') >= 0)):
            data2Dpathlist.append(d)
        else:
            data3Dpathlist.append(d)
    slitlist = []
    goldlist = []
    data2Dpathlist = []
    data3Dpathlist = []
    
    #load all the slits, save and make a library
    for slitpath in slitpathlist:
        slits = UBCSlits(path = os.path.join(path, slitpath))
        try:
            slits.Analyze() 
            shortname = slitpath[0:3]
            slits.save(os.path.join(path, shortname))
            pol = getpol(slitpath)
            slitlist.append((os.path.join(path,(shortname + '.pickle')), slits.paramdict['Ek'], slits.paramdict['Ep'], slits.paramdict['lens_type'], pol, slits.paramdict['tstamp']))
        except:
            print('failed to analyze %s' % shortname)
    print('Done with slits')
    print(slitlist)
    
    
    #load all the gold and make library too
    for goldpath in goldpathlist:
        gold = UBCGold(path = os.path.join(path, goldpath))
        pol = getpol(goldpath)
        slits = UBCSlits(pickcorrection(gold, slitlist, getpol(goldpath)))
    
        gold = gold.Correct(slits)
        gold.Analyze()
        shortname = goldpath[0:3]

        gold.save(os.path.join(path, shortname))

        goldlist.append((os.path.join(path,(shortname + '.pickle')), gold.paramdict['Ek'], gold.paramdict['Ep'], gold.paramdict['lens_type'], pol, gold.paramdict['tstamp']))

    print('Done with gold')
    print(goldlist) 
       
    for data2Dpath in data2Dpathlist:
        data = AaData2D(path = os.path.join(path, data2Dpath)) #for now just average for std
        pol = getpol(data2Dpath)
        slits = UBCSlits(pickcorrection(data, slitlist, pol))
        gold = UBCGold(pickcorrection(data, goldlist, pol))
    
        data = data.Correct(slits, gold)
        shortname = data2Dpath[0:3]

        data.save(os.path.join(path, shortname))
        data.saveH5(os.path.join(path, shortname))


    print('Done with data2D')

    for data3Dpath in data3Dpathlist:
        data = AaData3D(path = os.path.join(path, data3Dpath)) #for now just average for std
        pol = getpol(data3Dpath)
        slits = UBCSlits(pickcorrection(data, slitlist, pol))
        gold = UBCGold(pickcorrection(data, goldlist, pol))
    
        data = data.Correct(slits, gold)
        shortname = data3Dpath[0:3]

        data.save(os.path.join(path, shortname))
        data.saveH5(os.path.join(path, shortname))
      
    print('Done with data3D')


    return goldlist, slitlist




def pickcorrection(data, correctionlist, pol):
    shortlist = []
    Ek = data.paramdict['Ek']
    Ep = data.paramdict['Ep']
    lens_type = data.paramdict['lens_type']
    pol = pol #should be changed to 
    tstamp = data.paramdict['tstamp']
    for corr in correctionlist:
        if ((corr[1] == Ek) and (corr[2] == Ep) and (corr[3] == lens_type) and (corr[4] == pol)):
            shortlist.append(corr)

    if(len(shortlist) == 0):
        raise RuntimeError('Could not find suitable slit')

    tdiff = None    
    
    for corr in shortlist:
        if((tdiff == None) or (tdiff > abs(corr[5] - tstamp))):
            tdiff = abs(corr[5] - tstamp)
            bestcorr = corr[0]
            
    return bestcorr
    
    
def dictpickcor(data, cdict):
    shortlist = []
    Ek = data.Ek
    Ep = data.Ep
    lens_type = data.lens_type
    gamma = data.gamma #should be changed to 
    tstamp = data.tstamp
    
    for key in cdict.keys():
        if ((np.isclose(cdict[key].Ek,Ek)) and (np.isclose(cdict[key].Ep,Ep)) and (cdict[key].lens_type == lens_type) and (np.isclose(cdict[key].gamma,gamma))):
            shortlist.append(key)    

    if(len(shortlist) == 0):
        raise RuntimeError('Could not find suitable slit')

    tdiff = None    
    
    for key in shortlist:
        if((tdiff == None) or (tdiff > abs(cdict[key].tstamp - tstamp))):
            tdiff = abs(cdict[key].tstamp - tstamp)
            bestcorr = cdict[key]
            
    return bestcorr            
            
            
def getpol(path):
    if(path.find('spol') >= 0):
        pol = 'spol'
    elif(path.find('ppol') >= 0):
        pol = 'ppol'
    else:
        pol = 'unknown'  
    return pol
    
    