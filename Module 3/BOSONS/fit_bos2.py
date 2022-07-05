import numpy 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pylab
import math

n,y,dy = pylab.loadtxt('bosoni_min.txt',unpack = True)



init=(0,1,0.05)
#funzione di fit
def f1(n,a,b,c):
    return a + b*((pylab.cosh((c*n)/2)*(pylab.sinh(c*n))**2 +          2*(pylab.cosh(c*n)*(pylab.sinh((c*n)/2))**3))/((pylab.sinh(c*n)*(pylab.sinh((c*n)/2)))*(pylab.sinh(c*n) + 2*(pylab.sinh((c*n)/2))**2)))
    

popt,pcov=curve_fit(f1,n,y,init,dy,absolute_sigma=False)
a,b,c=popt
da,db,dc= pylab.sqrt(pcov.diagonal())
#eta = c*2
#deta = dc*2
Chi=(((y-f1(n,a,b,c))/dy)**2).sum()

ndof=len(n)-len(init)

resnorm=(y-f1(n,a,b,c))/dy
print('a =%f +- %f' %(a, da))
print('b =%f +- %f ' %(b, db))
print('c =%f +- %f ' %(c, dc))
print('cov',pcov)
print('Chi,ndof:',Chi,ndof)


#grafico
xx=numpy.linspace(0,600,1000000)
pylab.figure(1)
plt.errorbar(n,y,dy,linestyle='',color='blue',marker='x')
plt.rc('font',size=18)
plt.xlabel('beta')
plt.ylabel('$<U>$ ')
plt.minorticks_on

'''plt.xlim(0,3)
plt.ylim(0,5)'''
plt.grid(color='gray',ls='--')
plt.subplots_adjust()
plt.plot(xx,f1(xx,a,b,c), color='black')


pylab.show()
