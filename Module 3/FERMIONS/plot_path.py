import numpy 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pylab
import math

n,x = pylab.loadtxt('path1.txt',unpack = True)
nn, y = pylab.loadtxt('path2.txt',unpack = True)

path = 100

'''x11 = numpy.zeros(path)
x12 =  numpy.zeros(path)
x21 = numpy.zeros(path)
x22 =  numpy.zeros(path)
pas = numpy.zeros(path)'''
z1 = numpy.zeros(2*path)
z2 = numpy.zeros(2*path)

'''for k in range(0,path,1):
    x11[k] = x[k]
    x12[k] = x[k + path]
    x21[k] = y[k]
    x22[k] = y[k + path]
    pas[k] = n[k]'''
    
for k in range(0,path,1):
    z1[k]=x[k]
    z1[k+path]=y[k]
    z2[k]=x[k+path]
    z2[k+path]=y[k+path]
    

pylab.figure(1)
plt.errorbar(z1,n,linestyle='-',color='blue', label='path 1')
plt.errorbar(z2,n,linestyle='-',color='red', label = 'path 2')
plt.rc('font',size=18)
plt.xlabel('x[n]')
plt.ylabel('n')
plt.legend(loc='upper right')
plt.axhline(y=100, color='gray', linestyle='--')
'''plt.axvline(x=z1[0], color='powderblue', linestyle='--')
plt.axvline(x=z2[0], color='powderblue', linestyle='--')'''
plt.minorticks_on

'''pylab.figure(2)
plt.errorbar(x12,pas,linestyle='-',color='darkslategray',marker='.')
plt.errorbar(y22,pas,linestyle='-',color='red',marker='.')
plt.rc('font',size=18)xx
plt.minorticks_on'''





pylab.show()
