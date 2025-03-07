

I due codici scalar4d_simple.f e scalar2d_simple.f fanno la simulazione 
di un campo scalare neutro libero rispettivamente in 4 e 2 dimensioni,
per essere compilati necessitano di parameters_4d.f e parameters_2d.f.
Per essere eseguiti necessitano del file randomseed con dentro un numero 
intero da leggere e del file di input.

Ecco un benchmark per il caso in 2 dimensioni:

40 x 8    massa 0.04
1100000 (1.1M) misure con una spazzata di heatbath fra una misura e l'altra
scartiamo le prime 100K misure, otteniamo
il codice ci mette 22 secondi

< m^2 phi^2 >               = 0.00420(16)

< (phi(n + x) - phi(n))^2 > = 0.50441(10)

< (phi(n + t) - phi(n))^2 > = 0.49150(6)


se invece fra una misura e l'altra facciamo 
1 spazzata di heat-bath + 4 spazzate di overrelaxation
facendo lo stesso numero di misure e di scarto
il codice ci mette 27 secondi 
(quindi OR costa circa 20 volte in meno rispetto a HB)


< m^2 phi^2 >               = 0.004243(12)

< (phi(n + x) - phi(n))^2 > = 0.50423(7)

< (phi(n + t) - phi(n))^2 > = 0.49148(6)


il guadagno sull'errore per m^2 phi^2 e` un'ordine di grandezza


=====================================

CASO 4D

reticolo 8^3 x 4 con massa 0.1
110K misure con 5 spazzate (1 HB + 4 OR) fra una misura e l'altra
tempo totale 23 secondi:  

< m^2 phi^2 >                 = 0.002011(5)

< (phi(n + x) - phi(n))^2 > +
< (phi(n + y) - phi(n))^2 > +
< (phi(n + z) - phi(n))^2 >   = 0.75034(9)

< (phi(n + t) - phi(n))^2 >   = 0.24774(4)
