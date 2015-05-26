# -*- coding: utf-8 -*-
"""
Created on Sun May 24 22:08:51 2015

@author: Ron
"""
from numpy import *
import numpy as np
import matplotlib.pylab as mpl
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from loadvec import *
from vecClass import *


path = 'C:/Users/Ron/Desktop/Ron_LLK_PIV_data/Analysis'
fname = 'camera3_LLK001612.T000.D000.P052.H008.L.vec'

X,Y,U,V = vecToMatrix(get_data(fname,path))
dt = get_dt(fname,path)
vec = vec(X,Y,U,V,dt)
resolution = 1.0/71.96 #[mm/px]

vec.rotate(-90)
vec.scale(resolution)
vec.crop(-30,-5,15,30)
vec.getVelStat()
vec.filterVelocity('med')   # med / gauss
genQuiver(vec)
plt.figure()
genVelHist(vec)
plt.figure()
genVorticityMap(vec)