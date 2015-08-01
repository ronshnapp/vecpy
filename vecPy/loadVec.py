#!/usr/bin/pythonw
""" LOADVEC - load a PIV .vec file onto a python oriented envierment 

addapted from Alex Liberzon's code.
extended by Ron Shnapp 24.5.15


"""
import os
import numpy as np
from numpy import *
import matplotlib.pylab as mpl
import matplotlib.pyplot as plt
from vecPy import vec
from string import upper, lower

def get_dt(fname,path):
    """given a .vec file this will return the delta t 
    from the file in micro seconds"""
    # os.chdir(path) BUG
    fname = os.path.join(os.path.abspath(path),fname) # just make a full path name 
    # new way of opening and closing the file
    with open(fname) as f:
        header = f.readline()
        
    ind1 = header.find('MicrosecondsPerDeltaT')
    dt = float(header[ind1:].split('"')[1])
    return dt

def get_data(fname,path):
    """this function gathers and retuens the data found in
    a single .vec file"""
    fname = os.path.join(os.path.abspath(path),fname) # just make a full path name 
    if fname.lower().endswith('.vec'):
        data = np.genfromtxt(fname,skip_header=1,delimiter=',',usecols=(0,1,2,3,4))
    else:
        raise 'Wrong file extension'
        
    return data
    
def get_data_openpiv(fname,path):
    """this function gathers and retuens the data found in
    a single .txt file"""
    fname = os.path.join(os.path.abspath(path),fname)
    if fname.endswith('.txt'): 
        data = np.genfromtxt(fname,usecols=(0,1,2,3,4))
    else:
        raise 'Wrong file extension'
    
    return data
	
def read_directory(dirname, ext='vec'):
    # list_files = os.listdir(dirname)
    list_files = [s for s in os.listdir(dirname) if s.rpartition('.')[2] in (lower(ext),upper(ext))]
    return list_files
	
def patternize(lst):
    """helper function for vecToMatrix"""
    lst = sorted(lst)
    n = [lst[0]]
    for i in range(len(lst)-1):
        if lst[i] != lst[i+1]:
            n.append(lst[i+1])    
    return n
    
def vecToMatrix(data):
    """this function takes vector form data and shifts it in
    to a matrix form. return is 4 matrices:
    X,Y - x and y axis position matrices in meshgrid form
    U,V - velocity of the flow"""
    x = patternize(data[:,0])
    y = patternize(data[:,1])
    X,Y = meshgrid(x,y)
    u1 = reshape(data[:,2],shape(X))
    v1 = reshape(data[:,3],shape(Y))*-1
    chc = reshape(data[:,4],shape(X))
    
    return (X,Y,u1,v1,chc)

def vecToVec(fname,path,lUnits='m',tUnits='s'):
    """ 
    using vecToMatrix and get_data
    generate directly the vec class object
    """
    X,Y,U,V,CHC = vecToMatrix(get_data(fname,path))
    dt = get_dt(fname,path)
    vector = vec(X,Y,U,V,CHC,dt,lUnits=lUnits,tUnits=tUnits)
    return vector
    
        
def getVecList(path, resolution=1, LUnits='mm',crop=False,rotate=False,Filter=True):
    """
    this function returns a list of vec instances
    for each .vec file in directory 'path'
    """
    fnames = os.listdir(path)  
    lst = []
    for n in fnames:
        if '.vec' in n:
            X,Y,U,V,CHC = vecToMatrix(get_data(n,path))
            dt = get_dt(n,path)
            vector = vec(X,Y,U,V,CHC,dt,lUnits=LUnits)
            vector.scale(resolution)
            if Filter: vector.filterVelocity('med',5)
            if rotate: vector.rotate(rotate)
            if crop: vector.crop(crop[0],crop[1],crop[2],crop[3])
            lst.append(vector)
    return lst
     
def readTimeStamp(fname,path):
    """reads an insight tstmp file and returns
    an array of the times at which photos were
    taken at relative to the begining of
    aquasition"""
    fname = os.path.join(os.path.abspath(path),fname)
    num_lines = sum(1 for line in open(fname))
    f = open(fname)
    for i in range(3):
        f.readline()
    strt = [f.readline().split()[1] for i in range(num_lines-4)]
    print strt[0:3]
    t = [float(i)/1000000 for i in strt]
    return t