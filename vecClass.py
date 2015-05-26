# -*- coding: utf-8 -*-
"""
Created on Sun May 24 22:02:49 2015

@author: Ron
"""
from numpy import *
from scipy.stats import norm
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import median_filter

class vec:
    
    def __init__(self,X,Y,U,V,dt):
        self.X = X
        self.Y = Y
        self.U = U
        self.V = V
        self.dt = dt
        
    def rotate(self,theta):
        """ 
        use this method in order to rotate the data 
        by theta degrees in the clock wise direction
        """
        theta = theta/360.0*2*pi
        xi = self.X*cos(theta)+self.Y*sin(theta)
        eta = self.Y*cos(theta)-self.X*sin(theta)
        Uxi = self.U*cos(theta)+self.V*sin(theta)
        Ueta = self.V*cos(theta)-self.U*sin(theta)
        self.X = xi
        self.Y = eta
        self.U = Uxi
        self.V = Ueta
        
    def scale(self,resolution):
        """
        use this method to fux the resilution of the 
        vector from [px/frame] to [m/sec] or any similar
        - resolution should be in [length/px]
        - time is generated from the original file and
          it is in seconds
        """
        self.X=self.X*resolution
        self.Y=self.Y*resolution
        self.U=self.U*resolution/(self.dt*1e-6)
        self.V=self.V*resolution/(self.dt*1e-6)
        
    def move(self,dx,dy):
        """
        use this method to move the origin of the frame
        by dx and dy
        """
        self.X=self.X+dx
        self.Y=self.Y+dy
        
    def crop(self,xmin,xmax,ymin,ymax):
        """
        this method is used to crop a rectangular section 
        of the vector field difined as the region between 
        (xmin,ymin) and (xmax,ymax) 
        """
        x = self.X[:,0]
        y = self.Y[0,:]
        il,ih,jl,jh = -1,-1,-1,-1,
        for i in range(len(x)-1):
            if x[i]>xmin and x[i+1]<xmin:
                ih = i
            elif x[i]>xmax and x[i+1]<xmax:
                il = i
        for j in range(len(y)-1):
            if y[j]<ymin and y[j+1]>ymin:
                jl = j
            elif y[j]<ymax and y[j+1]>ymax:
                jh = j
        for i in [il,ih,jl,jh]: 
            if i==-1: 
                print 'error with limits the were given'
                return
        self.X = self.X[il:ih,jl:jh]
        self.Y = self.Y[il:ih,jl:jh]
        self.U = self.U[il:ih,jl:jh]
        self.V = self.V[il:ih,jl:jh]
        
    def getVelStat(self):
        """
        assuming normaly distributed values of U and V in the
        data, this method calculates its mean and standard
        deviation values and assigns them to new atribtes
        of the instance vec.
        """
        self.Umean, self.Ustd = norm.fit(self.U)
        self.Vmean, self.Vstd = norm.fit(self.V)
        
    def filterVelocity(self,filtr = 'med'):
        """
        this method passes the velocity vectors U and V
        through a either a 4X4 median filter or a gaussian
        filter with sigma = 1
        """
        if filtr == 'med':
            self.U = median_filter(self.U,size=(4,4))
            self.V = median_filter(self.V,size=(4,4))
        elif filtr == 'gauss':
            self.U = gaussian_filter(self.U,1)
            self.V = gaussian_filter(self.V,1)
        else: print "Bad choise of filter! - try again"
        