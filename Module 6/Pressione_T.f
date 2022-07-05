      
	  
	  
	  program pressione
	  
	  implicit real (a-h,o-z)
	  implicit integer (i-n)
	  parameter(nx = 25, xmass = 0.07)
	  real r(nx),s(nx),t(nx),x(nx)
	  real rs(nx), xdt(nx)
	  
	  open(3,file='data_anomalia.txt',status='old')
	  
CC    numero a cui mi fermo nell'integrazione, nel mio file devo mettere gli Nt in ordine decrescente perchè devo integrare la T da 0 a T: il numero è è nx
      

CC   dati che nella prima colonna presentano Nt e nella seconda il valore (e - p)/T**2 e nella terza l'errore di questa quantità
      do j= 1, nx
	   read(3,*) rs(j),s(j),t(j)
	  enddo
	  
CC   esprimo tutto nei termini di T/m	 
	  do i= 1, nx
	   r(i) = 1/(rs(i)*xmass) 
	  enddo
	  xav = 0.0
	  xxav = 0.0
	  do i=1,nx-1
		 x(i) = s(i)/r(i)
		 x(i) = x(i) + s(i+1)/r(i+1)
		 x(i)=x(i)*(r(i+1) - r(i))
		 x(i) = x(i)/float(2)
		 
c		 x(i) = ((s(i)/(r(i))) + (s(i + 1)/(r(i+1))))*(r(i+1) - r(i))/float(2)
         xdt(i) = sqrt((t(i)/(r(i)))**2 + (t(i+1)/(r(i + 1)))**2)
		 xdt(i) = (((r(i+1) - r(i)))/float(2))*xdt(i)
		 
		 xav = xav + x(i)
		 xxav = xxav + xdt(i)**2
		 
	  enddo
	  xxav = sqrt(xxav)
	    
c	  xav = xav + s(nx)
c	  xxav = sqrt(xxav**2 + t(nx)**2)
	  
	  write(*,*) xav, xxav
	  close(3)
	  return
      end