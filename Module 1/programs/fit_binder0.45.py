import numpy 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pylab


#dati

j,o,do=pylab.loadtxt('binder0.45.txt',unpack=True)
#n serve a selezionare il range dei punti relativi al picco da cambiare di volta in volta a seconda del file e di L

n=12
#array di appoggio su cui inserire i punti selezionati
x=numpy.zeros(n)

    

for k in range(0,n,1):
    x[k]=1/(j[k])

    
init=(1,-1)
#funzione di fit
def f1(x,a,c):
    return a + c*numpy.sqrt(x)
    
popt,pcov=curve_fit(f1,x,o,init,do,absolute_sigma=True)
a,c=popt
da,dc= pylab.sqrt(pcov.diagonal())

Chi=(((o-f1(x,a,c))/do)**2).sum()

ndof=len(x)-len(init)

resnorm=(o-f1(x,a,c))/do
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
plt.errorbar(x,o,do,linestyle='',color='darkslategray',marker='.')
plt.rc('font',size=18)
plt.xlabel('$1/L$ ')
plt.ylabel('BC')
plt.minorticks_on

plt.xlim(0,0.30)
plt.ylim(1,1.25)
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


