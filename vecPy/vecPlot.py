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


def genQuiver(vec):
    """
    this function will generate a quiver plot from a vec 
    object    
    """
    S = sqrt(vec.u**2+vec.v**2)
    levels = np.linspace(0, amax(S), 30)
    mpl.contourf(vec.x,vec.y,S,alpha=0.5,
                 cmap=plt.cm.get_cmap("Blues"), 
                 levels=levels)
    cbar = mpl.colorbar()
    cbar.set_label(r'Velocity [m $\cdot$ s$^{-1}$]')
    mpl.quiver(vec.x,vec.y,vec.u,vec.v,units='width',scale=amax(sqrt(vec.u**2+vec.v**2))*25.0,headwidth=2 )
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
    
    
def genVorticityMap(vec):
    """ why do we rotate the vector before taking derivative? """
    # BUG:
    dUy = gradient(vec.u)[0]*cos(vec.theta)-gradient(vec.u)[1]*sin(vec.theta)
    dVx = gradient(vec.v)[1]*cos(vec.theta)+gradient(vec.v)[0]*sin(vec.theta)
    dx = gradient(vec.x)[1]*cos(vec.theta)+gradient(vec.x)[0]*sin(vec.theta)
    dy = gradient(vec.y)[0]*cos(vec.theta)-gradient(vec.y)[1]*sin(vec.theta)
    vorticity = dVx/dy-dUy/dx
    m = amax(absolute(vorticity))
    levels = np.linspace(-m, m, 30)
    plt.contourf(vec.x,vec.y,vorticity, levels=levels,
                 cmap = plt.cm.get_cmap('RdYlBu'))
    plt.xlabel('x [' + vec.lUnits + ']')
    plt.ylabel('y [' + vec.lUnits + ']')
    cbar = plt.colorbar()
    cbar.set_label(r'Vorticity [s$^{-1}$]')


def genShearMap(vec):
    """this function plots a map of the xy strain e_xy"""
    dUy = gradient(vec.u)[0]*cos(vec.theta)-gradient(vec.u)[1]*sin(vec.theta)
    dVx = gradient(vec.v)[1]*cos(vec.theta)+gradient(vec.v)[0]*sin(vec.theta)
    dx = gradient(vec.x)[1]*cos(vec.theta)+gradient(vec.x)[0]*sin(vec.theta)
    dy = gradient(vec.y)[0]*cos(vec.theta)-gradient(vec.y)[1]*sin(vec.theta)
    strain = dVx/dy+dUy/dx
    m = amax(absolute(strain))
    levels = np.linspace(-m, m, 30)
    plt.contourf(vec.x,vec.y,strain, levels=levels,
                 cmap = plt.cm.get_cmap('PRGn'))
    plt.xlabel('x [' + vec.lUnits + ']')
    plt.ylabel('y [' + vec.lUnits + ']')
    cbar = plt.colorbar()
    cbar.set_label(r'strain [s$^{-1}$]')
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