# -*- coding: utf-8 -*-
"""
Created on Thu May 28 12:30:24 2015

@author: Ron
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
plt.rcParams['animation.ffmpeg_path'] = 'C:/ffmpeg/bin/ffmpeg.exe'

import sys
modulepath = 'C:/Users/Ron/Desktop/work/torbulece_flow_lab/vector_piv'
sys.path.append(modulepath)
from vecpy import *


def animate(vecList, arrowscale=1, frameWait = 100):
    
    X, Y = vecList[0].x, vecList[0].y
    U, V = vecList[0].u, vecList[0].v
    fig, ax = plt.subplots(1,1)
    Q = ax.quiver(X, Y, U, V, units='inches', scale=arrowscale)
    ax.text(amax(X)-10, amax(Y)+2, 'test')
    def update_quiver(num,Q,vecList):
        Q.set_UVC(vecList[num].u,vecList[num].v)
        return Q,
    anim = animation.FuncAnimation(fig, update_quiver, fargs=(Q,vecList),
                               frames = len(vecList), interval=frameWait, blit=False)
    mywriter = animation.FFMpegWriter()
    anim.save('im.mp4', writer=mywriter)


def getVecList(path, resolution=1, LUnits='mm',crop=False,rotate=False,Filter=True):
    """
    this function returns a list of vec instances
    for each .vec file in directory 'path'
    """
    fnames = os.listdir(path)  
    lst = []
    for n in fnames:
        X,Y,U,V,CHC = vecToMatrix(get_data(n,path))
        dt = get_dt(n,path)
        vector = vec(X,Y,U,V,CHC,dt,lUnits=LUnits)
        vector.scale(resolution)
        if Filter: vector.filterVelocity('med',5)
        if rotate: vector.rotate(rotate)
        if crop: vector.crop(crop[0],crop[1],crop[2],crop[3])
        lst.append(vector)
    return lst
    
    



