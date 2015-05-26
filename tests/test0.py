#!/usr/bin/python

from vecpy import loadvec
from vecpy.vecClass import vec
import matplotlib.pyplot as plt


test_dir = "tests/data"
lst = loadvec.read_directory(test_dir)
data = loadvec.get_data(lst[0],test_dir)
dt = loadvec.get_dt(lst[0],test_dir)
X,Y,U,V,chc = loadvec.vecToMatrix(data)
v = vec(X,Y,U,V,chc,dt,Lunits='mm')


resolution = 1.0/71.96 #[mm/px]
v.rotate(-90)
v.scale(resolution)

plt.figure()
loadvec.genQuiver(v)
plt.show()

plt.figure()
loadvec.genVorticityMap(v)
plt.show()

