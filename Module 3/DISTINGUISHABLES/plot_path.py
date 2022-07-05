import numpy 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pylab
import math

n,x = pylab.loadtxt('path1.txt',unpack = True)
nn, y = pylab.loadtxt('path2.txt',unpack = True)

path = 50

x11 = numpy.zeros(path)
x12 =  numpy.zeros(path)
y11 = numpy.zeros(path)
y12 =  numpy.zeros(path)
pas = numpy.zeros(path)

for k in range(0,path,1):
    x11[k] = x[k]
    x12[k] = x[k + path]
    y11[k] = y[k]
    y12[k] = y[k + path]
    pas[k] = n[k]
    


pylab.figure(1)
plt.errorbar(x11,pas,linestyle='-',color='blue', label='1st cicle')
plt.errorbar(y11,pas,linestyle='-',color='red', label = '2nd cicle')
plt.rc('font',size=18)
plt.xlabel('x[n]')
plt.ylabel('n')
plt.legend(loc='upper right')
plt.minorticks_on

'''pylab.figure(2)
plt.errorbar(x12,pas,linestyle='-',color='darkslategray',marker='.')
plt.errorbar(y12,pas,linestyle='-',color='red',marker='.')
plt.rc('font',size=18)xx
plt.minorticks_on'''





pylab.show()
