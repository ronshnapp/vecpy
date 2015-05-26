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

    
def genQuiver(vec):
    """
    this function will generate a quiver plot from a vec 
    object    
    """
    mpl.contourf(vec.X,vec.Y,sqrt(vec.U**2+vec.V**2),alpha=0.3)
    cbar = mpl.colorbar()
    cbar.set_label(r'Velocity [m $\cdot$ s$^{-1}$]')
    mpl.quiver(vec.X,vec.Y,vec.U,vec.V,scale=2000)
    mpl.xlabel('x ['+vec.length+']')
    mpl.ylabel('y ['+vec.length+']')
    return
    
def fit2LLK(vec):
        """this function is specialy made for Ron's 
        LLK PIV experiment"""
        vec.rotate(-90)
        return
    
def genVelHist(vec):
    """
    this function will plot a normalized histogram of
    the velocity data. the velocity in the histpgram 
    contains all the data set such that CHC == 1
    """
    u,v,chc = vec.U.flatten() , vec.V.flatten(), vec.CHC.flatten()
    print shape(u),shape(v),shape(chc)
    u1, v1 = [], []
    for i in range(size(chc)):
        if chc[i]==1:
            u1.append(u[i])
            v1.append(v[i])
    ax1 = plt.subplot2grid((2,1),(0,0))
    ax1.hist(u,bins=np.sqrt(len(u))*0.5,normed=1)
    ax1.set_xlabel('u [mm/sec]')
    ax2 = plt.subplot2grid((2,1),(1,0))
    ax2.hist(v,bins=np.sqrt(len(v)*0.5),normed=1)
    ax2.set_xlabel('v [mm/sec]')
    plt.tight_layout()
    return
    
def genVorticityMap(vec):
    dUy = gradient(vec.U)[0]*cos(vec.theta)-gradient(vec.U)[1]*sin(vec.theta)
    dVx = gradient(vec.V)[1]*cos(vec.theta)+gradient(vec.V)[0]*sin(vec.theta)
    dx = gradient(vec.X)[1]*cos(vec.theta)+gradient(vec.X)[0]*sin(vec.theta)
    dy = gradient(vec.Y)[0]*cos(vec.theta)-gradient(vec.Y)[1]*sin(vec.theta)
    print amax(dUy),amax(dVx),amax(dx),amax(dy)
    vorticity = dVx/dy-dUy/dx
    plt.contourf(vec.X,vec.Y,vorticity)
    plt.xlabel('x ['+vec.length+']')
    plt.ylabel('y ['+vec.length+']')
    cbar = plt.colorbar()
    cbar.set_label(r'Vorticity [s$^{-1}$]')
    return