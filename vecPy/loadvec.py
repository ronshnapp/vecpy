#!/usr/bin/pythonw
""" LOADVEC - load a PIV .vec file onto a python oriented envierment 

addapted from Alex Liberzon's code.
extended by Ron Shnapp 24.5.15

this module is dedicated to loading the data within a .vec file
and using it to create instances of a vec object.
"""
import os
import numpy as np
from numpy import *


def get_dt(fname,path):
    """given a .vec file this will return the delta t 
    from the file in micro seconds"""
    os.chdir(path)
    f = open(fname)
    header = f.readline()
    f.close()
    ind1 = header.find('MicrosecondsPerDeltaT')
    dt = float(header[ind1:].split('"')[1])
    return dt

def get_data(fname,path):
    """this function gathers and retuens the data found in
    a single .vec file"""
    os.chdir(path)
    data = np.genfromtxt(fname,skip_header=1,delimiter=',',usecols=(0,1,2,3,4))
    return data
	
def read_directory(dirname):
    list_files = os.listdir(dirname)
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

    