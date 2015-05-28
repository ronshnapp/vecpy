#!/usr/bin/python

from vecPy import loadvec
from vecPy import vecPlot
from vecPy.vecClass import vec
import matplotlib.pyplot as plt


test_dir = "tests/data"
lst = loadvec.read_directory(test_dir)
data = loadvec.get_data(lst[3],test_dir)
dt = loadvec.get_dt(lst[3],test_dir)
x,y,u,v,chc = loadvec.vecToMatrix(data)
vec = vec(x,y,u,v,chc,dt,lUnits='mm',tUnits = 's')


resolution = 1.0/71.96 #[mm/px]
vec.rotate(-90)
vec.scale(resolution)

plt.figure()
vecPlot.genQuiver(vec)
plt.show()

plt.figure()
vecPlot.genVorticityMap(vec)
plt.show()

