# -*- coding: utf-8 -*-
"""
Created on Wed May 27 14:04:24 2015

@author: Ron

This module contains function that help generating
velocity attributes of the .vec files that were open.
use this code over vec objects.
"""

import numpy as np
from numpy import *
from matplotlib.pylab import *
import matplotlib.pylab as mpl
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from os import chdir, getcwd, listdir
plt.rcParams['animation.ffmpeg_path'] = 'C:/ffmpeg/bin/ffmpeg.exe'
from scipy.ndimage.filters import median_filter

def genQuiver(vec, arrScale = 25.0, threshold = None, nthArr = 1):
    """
    this function will generate a quiver plot from a vec 
    object   
    threshold - values above the threshold will be set equal to threshold
    arrScale - use to change arrow scales
    nthArr - use to plot only every nth arrow from the array 
    """
    u = vec.u
    v = vec.v
    if threshold != None:
        u = thresholdArray(u, threshold)
        v = thresholdArray(v, threshold)
        print 1
    S = sqrt(u**2+v**2)
    levels = np.linspace(0, amax(S), 30)
    mpl.contourf(vec.x,vec.y,S,alpha=0.8,
                 cmap=plt.cm.get_cmap("Blues"), 
                 levels=levels)
    cbar = mpl.colorbar()
    cbar.set_label(r'Velocity ['+vec.lUnits+' $\cdot$ '+vec.tUnits+'$^{-1}$]')
    n = nthArr
    mpl.quiver(vec.x[1::n,1::n],vec.y[1::n,1::n],
               u[1::n,1::n],v[1::n,1::n],units='width',
               scale=amax(S*arrScale),headwidth=2 )
    mpl.xlabel('x [' + vec.lUnits + ']')
    mpl.ylabel('y [' + vec.lUnits + ']')
    
    
def genFluctuationQuiver(vec):
    """
    generate a quiver plot of velocity fluctuation
    i.e. velocity-mean_velocity.    
    """
    vec.getVelStat()
    u,v = vec.u-vec.Umean,vec.v-vec.Vmean
    S = sqrt(u**2+v**2)
    levels = np.linspace(0, amax(S), 30)
    mpl.contourf(vec.x,vec.y,S,alpha=0.5,
                 cmap=plt.cm.get_cmap("Greens"),
                 levels=levels)
    cbar = mpl.colorbar()
    cbar.set_label(r'Velocity [m $\cdot$ s$^{-1}$]')
    mpl.quiver(vec.x, vec.y, u, v, units='width',scale=amax(sqrt(u**2+v**2))*25.0,headwidth=2 )
    mpl.xlabel('x [' + vec.lUnits + ']')
    mpl.ylabel('y [' + vec.lUnits + ']')
    
    
def genVelHist(vec):
    """
    this function will plot a normalized histogram of
    the velocity data.
    """
    u1, v1 = vec.u.flatten(), vec.v.flatten()
    ax1 = plt.subplot2grid((2,1),(0,0))
    ax1.hist(u1,bins=np.sqrt(len(u1))*0.5,normed=1)
    ax1.set_xlabel('u [mm/sec]')
    ax2 = plt.subplot2grid((2,1),(1,0))
    ax2.hist(v1,bins=np.sqrt(len(v1)*0.5),normed=1)
    ax2.set_xlabel('v'+vec.velUnits)
    plt.tight_layout()
    
    
def genVorticityMap(vec, threshold = None):
    """ why do we rotate the vector before taking derivative? """
    # BUG:
    dUy = gradient(vec.u)[0]*cos(vec.theta)-gradient(vec.u)[1]*sin(vec.theta)
    dVx = gradient(vec.v)[1]*cos(vec.theta)+gradient(vec.v)[0]*sin(vec.theta)
    dx = gradient(vec.x)[1]*cos(vec.theta)+gradient(vec.x)[0]*sin(vec.theta)
    dy = gradient(vec.y)[0]*cos(vec.theta)-gradient(vec.y)[1]*sin(vec.theta)
    vorticity = dVx/dy-dUy/dx
    if threshold != None:
        vorticity = thresholdArray(vorticity,threshold)
    m = amax(absolute(vorticity))
    levels = np.linspace(-m, m, 30)
    plt.contourf(vec.x,vec.y,vorticity, levels=levels,
                 cmap = plt.cm.get_cmap('RdYlBu'))
    plt.xlabel('x [' + vec.lUnits + ']')
    plt.ylabel('y [' + vec.lUnits + ']')
    cbar = plt.colorbar()
    cbar.set_label(r'Vorticity [s$^{-1}$]')


def genShearMap(vec, threshold = None):
    """this function plots a map of the xy strain e_xy"""
    dUy = gradient(vec.u)[0]*cos(vec.theta)-gradient(vec.u)[1]*sin(vec.theta)
    dVx = gradient(vec.v)[1]*cos(vec.theta)+gradient(vec.v)[0]*sin(vec.theta)
    dx = gradient(vec.x)[1]*cos(vec.theta)+gradient(vec.x)[0]*sin(vec.theta)
    dy = gradient(vec.y)[0]*cos(vec.theta)-gradient(vec.y)[1]*sin(vec.theta)
    strain = dVx/dy+dUy/dx
    if threshold != None:
        strain = thresholdArray(strain,threshold)
    m = amax(absolute(strain))
    levels = np.linspace(-m, m, 30)
    plt.contourf(vec.x,vec.y,strain, levels=levels,
                 cmap = plt.cm.get_cmap('PRGn'))
    plt.xlabel('x [' + vec.lUnits + ']')
    plt.ylabel('y [' + vec.lUnits + ']')
    cbar = plt.colorbar()
    cbar.set_label(r'strain [s$^{-1}$]')
    

def genFlowAcceleration(vec, arrScale = 25.0, threshold = None, nthArr = 1):
    """this function will plot a contour plot of
    the convective term of material derivative.
    i.e. it plots the magnitude of the vector 
    (u*dudx + v*dudy , u*dvdx + v*dvdy)"""  
    dUx = gradient(vec.u)[1]*cos(vec.theta)+gradient(vec.u)[0]*sin(vec.theta)
    dUy = gradient(vec.u)[0]*cos(vec.theta)-gradient(vec.u)[1]*sin(vec.theta)
    dVx = gradient(vec.v)[1]*cos(vec.theta)+gradient(vec.v)[0]*sin(vec.theta)
    dVy = gradient(vec.v)[0]*cos(vec.theta)-gradient(vec.v)[1]*sin(vec.theta)
    dx = gradient(vec.x)[1]*cos(vec.theta)+gradient(vec.x)[0]*sin(vec.theta)
    dy = gradient(vec.y)[0]*cos(vec.theta)-gradient(vec.y)[1]*sin(vec.theta)
    ax = median_filter(vec.u*dUx/dx + vec.v*dUy/dy , (3,3))
    ay = median_filter(vec.u*dVx/dx + vec.v*dVy/dy , (3,3))
    if threshold != None:
        ax = thresholdArray(ax,threshold)
        ay = thresholdArray(ay,threshold)
    S = sqrt(ax**2+ay**2)
    levels = np.linspace(0, amax(S), 30)
    mpl.contourf(vec.x,vec.y,S,alpha=0.5,
                 cmap=plt.cm.get_cmap("OrRd"), 
                 levels=levels)
    cbar = mpl.colorbar()
    cbar.set_label(r'Spatial Material Derivative ['+vec.lUnits+' $\cdot$ '+vec.tUnits+'$^{-2}$]')
    n = nthArr
    mpl.quiver(vec.x[1::n,1::n],vec.y[1::n,1::n],
               ax[1::n,1::n],ay[1::n,1::n],units='width',
               scale=amax(S)*arrScale,headwidth=2 )
    mpl.xlabel('x [' + vec.lUnits + ']')
    mpl.ylabel('y [' + vec.lUnits + ']')
    
    
def animateVecList(vecList, arrowscale=1, savepath=None):
    X, Y = vecList[0].x, vecList[0].y
    U, V = vecList[0].u, vecList[0].v
    fig, ax = plt.subplots(1,1)
    #Q = ax.quiver(X, Y, U, V, units='inches', scale=arrowscale)
    M = sqrt(pow(U, 2) + pow(V, 2))    
    Q = ax.quiver(X[::3,::3], Y[::3,::3], 
                  U[::3,::3], V[::3,::3], M[::3,::3],
                 units='inches', scale=arrowscale)
    cb = plt.colorbar(Q)
    cb.ax.set_ylabel('velocity ['+vecList[0].lUnits+'/'+vecList[0].tUnits+']')
    text = ax.text(0.2,1.05, '1/'+str(len(vecList)), ha='center', va='center', transform=ax.transAxes)
    def update_quiver(num,Q,vecList,text):
        U,V = vecList[num].u[::3,::3],vecList[num].v[::3,::3]
        M = sqrt(pow(U, 2) + pow(V, 2))   
        Q.set_UVC(U,V,M)
        #Q.set_UVC(vecList[num].u,vecList[num].v)
        text.set_text(str(num+1)+'/'+str(len(vecList)))
        return Q,
    anim = animation.FuncAnimation(fig, update_quiver, fargs=(Q,vecList,text),
                               frames = len(vecList), blit=False)
    mywriter = animation.FFMpegWriter()
    if savepath:
        p = getcwd()
        chdir(savepath)
        anim.save('im.mp4', writer=mywriter)
        chdir(p)
    else: anim.save('im.mp4', writer=mywriter)  
    
def thresholdArray(array, th):
    index = where(absolute(array)>th)
    for i in range(len(index[0])):
        array[index[0][i],index[1][i]] = th*sign(array[index[0][i],index[1][i]])
    return array
    