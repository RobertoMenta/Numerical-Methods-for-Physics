import numpy 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pylab
import math

#dati

j1,o1,do1=pylab.loadtxt('L20.txt',usecols=(0,3,4),unpack=True)
j2,o2,do2=pylab.loadtxt('L30.txt',usecols=(0,3,4),unpack=True)
j3,o3,do3=pylab.loadtxt('L40.txt',usecols=(0,3,4),unpack=True)
j4,o4,do4=pylab.loadtxt('L50.txt',usecols=(0,3,4),unpack=True)
j5,o5,do5=pylab.loadtxt('L60.txt',usecols=(0,3,4),unpack=True)
j6,o6,do6=pylab.loadtxt('L70.txt',usecols=(0,3,4),unpack=True)


b= 0.44
v = 0

L1=20
L2=30
L3=40
L4=50
L5=60
L6=70


    
    

    
pylab.figure(1)

plt.errorbar((j1-b)*L1,o1/((L1)**v),do1/((L1)**v),linestyle='',color='darkslategray',marker='.',fillstyle='none',label='L=20')
plt.errorbar((j2-b)*L2,o2/((L2)**v),do2/((L2)**v),linestyle='',color='slateblue',marker='.',fillstyle='none',label='L=30')
plt.errorbar((j3-b)*L3,o3/((L3)**v),do3/((L3)**v),linestyle='',color='turquoise',marker='.',fillstyle='none',label='L=40')
plt.errorbar((j4-b)*L4,o4/((L4)**v),do4/((L4)**v),linestyle='',color='gold',marker='.',fillstyle='none',label='L=50')
plt.errorbar((j5-b)*L5,o5/((L5)**v),do5/((L5)**v),linestyle='',color='lightgreen',marker='.',fillstyle='none',label='L=60')
plt.errorbar((j6-b)*L6,o6/((L6)**v),do6/((L6)**v),linestyle='',color='lightcoral',marker='.',fillstyle='none',label='L=70')


plt.rc('font',size=12)
plt.ylabel('$C \\vert/L^{\\alpha/\\nu}$ ')
plt.xlabel('$(\\beta - \\beta_c)L^{1/\\nu}$')
plt.minorticks_on
plt.legend(loc='best')

'''plt.errorbar(j1,o1,linestyle='',color='darkslategray',marker='.',fillstyle='none',label='L=20')
plt.errorbar(j2,o2,linestyle='',color='slateblue',marker='.',fillstyle='none',label='L=30')
plt.errorbar(j3,o3,linestyle='',color='turquoise',marker='.',fillstyle='none',label='L=40')
plt.errorbar(j4,o4,linestyle='',color='gold',marker='.',fillstyle='none',label='L=50')
plt.errorbar(j5,o5,linestyle='',color='lightgreen',marker='.',fillstyle='none',label='L=60')
plt.errorbar(j6,o6,linestyle='',color='lightcoral',marker='.',fillstyle='none',label='L=70')




plt.rc('font',size=12)
plt.ylabel('$\\langle \\vert M \\vert \\rangle$')
plt.xlabel('$\\beta$')
plt.minorticks_on
plt.legend(loc='best',)'''




pylab.show()
