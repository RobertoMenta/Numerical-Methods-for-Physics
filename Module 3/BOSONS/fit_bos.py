import numpy 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pylab
import math

n,y,dy = pylab.loadtxt('bosoni_min.txt',unpack = True)

x = 1/(n*0.05)

init=(0,1)
#funzione di fit
def f1(x,a,b):
    return a + b*((pylab.cosh(1/(2*x))*(pylab.sinh(1/x))**2 + 2*(pylab.cosh(1/x)*(pylab.sinh(1/(2*x)))**3))/((pylab.sinh(1/x)*(pylab.sinh(1/(2*x))))*(pylab.sinh(1/x) + 2*(pylab.sinh(1/(2*x)))**2)))
    

popt,pcov=curve_fit(f1,x,y,init,dy,absolute_sigma=False)
a,b=popt
da,db= pylab.sqrt(pcov.diagonal())
#eta = c*2
#deta = dc*2
Chi=(((y-f1(x,a,b))/dy)**2).sum()

ndof=len(x)-len(init)

resnorm=(y-f1(x,a,b))/dy
print('a =%f +- %f' %(a, da))
print('b =%f +- %f ' %(b, db))

#print ('eta = %f +- %f'%(eta,deta))
print('cov',pcov)
print('Chi,ndof:',Chi,ndof)


'''m= len(x)
T = numpy.zeros(m)
for k in range(0,m,1):
    T[k]= 1/(eta*x[k])

j = numpy.zeros(1000000)

for k in range(0,1000000,1):
    j[k]= 1/(0.01 + k/(1000000))
  
p=0
s=1
u=0.0025'''

#grafico
xx=numpy.linspace(0,2,1000000)
pylab.figure(1)
plt.errorbar(x,y,dy,linestyle='',color='darkslategray',marker='.')
plt.rc('font',size=18)
plt.xlabel('T')
plt.ylabel('$<E>$ ')
plt.minorticks_on

'''plt.xlim(0,3)
plt.ylim(0,5)'''
plt.subplots_adjust()
plt.plot(xx,f1(xx,a,b), color='skyblue')
plt.axhline(y=1, color='powderblue', linestyle='--')



pylab.show()
