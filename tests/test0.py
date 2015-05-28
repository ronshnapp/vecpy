#!/usr/bin/python

from vecpy import loadvec
from vecpy import vecplot
from vecpy.vecpy import vec
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
vecplot.genQuiver(vec)
plt.show()

plt.figure()
vecplot.genVorticityMap(vec)
plt.show()

