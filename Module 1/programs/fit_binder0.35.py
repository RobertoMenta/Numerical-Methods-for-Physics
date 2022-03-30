import numpy 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pylab
import math

#dati

j,o,do=pylab.loadtxt('binder0.35.txt',unpack=True)
#n serve a selezionare il range dei punti relativi al picco da cambiare di volta in volta a seconda del file e di L

n=10
#array di appoggio su cui inserire i punti selezionati
x=numpy.zeros(n-2)
y=numpy.zeros(n-2)
dy=numpy.zeros(n-2)
    

for k in range(2,n,1):
    x[k-2]=1/(j[k])
    y[k-2]=o[k]
    dy[k-2]=do[k]
    

    
init=(1,-1)
#funzione di fit
def f1(x,a,c):
    return a + c*x**2
    
popt,pcov=curve_fit(f1,x,y,init,dy,absolute_sigma=True)
a,c=popt
da,dc= pylab.sqrt(pcov.diagonal())

Chi=(((y-f1(x,a,c))/dy)**2).sum()

ndof=len(x)-len(init)

resnorm=(y-f1(x,a,c))/dy
covn=numpy.zeros(1)
covn=pcov[0][1]/(pylab.sqrt(pcov[0][0]*pcov[1][1]))
print('cumulante di binder =%f +- %f' %(a, da))
print('c =%f +- %f ' %(c, dc))
#print('cov',pcov)
print('Chi,ndof:',Chi,ndof)
#print('covarianza normalizzata',covn)


#grafico
xx=numpy.linspace(0,0.50,1000)

'''plt.subplot(2,1,1)'''
pylab.figure(1)
plt.errorbar(1/j,o,do,linestyle='',color='darkslategray',marker='.')
plt.rc('font',size=18)
plt.xlabel('$1/L$ ')
plt.ylabel('BC')
plt.minorticks_on

plt.xlim(0,0.21)
plt.ylim(0,3.5)
plt.subplots_adjust()
plt.plot(xx,f1(xx,a,c), color='skyblue')


'''plt.subplot(2,1,2)
plt.errorbar(x,resnorm,dy,linestyle='',color='red',marker='.')
plt.rc('font',size=18)
plt.ylim(-6,6)
plt.xlabel('$coordinata x$ [bo]')
plt.ylabel('Resnorm')
plt.minorticks_on
plt.grid(color='gray')
plt.plot(x,resnorm,dy,linestyle='--', color='blue',marker='o')'''

pylab.show()

