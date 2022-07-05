CC####################################################################
      program scalar 2d
CC####################################################################
c programma per la simulazione del campo scalare in 2 D
C--------------------------------------------------------------------
      parameter(nx = 100,nt = 50 ,nvol = nx*nt)
      parameter(nlatt = nx)   !! must be the maximum dimension size
      parameter(pigr = 3.141592654)
      real mass,extfield,mass2,mass2p4  
      common/various/mass,extfield,mass2,mass2p4 !! parametri della simulazione
      real field
      common/lattice/field(nx,nt) !! matrice con le variabili di spin
                                        !! che rappresenta lo stato del sistema
      common/move/npp(nlatt,2),nmm(nlatt,2) 


C      include "parameters_2d.f"

      integer*8 seed,seed_old



CC    per far girare il codice dopo averlo compilato, e` necessario che 
CC    esista un file "input" contenente i 5 numeri (3 interi e 2 reali) letti
CC    poche righe piu` sotto, ed un file "randomseed" contenente lo stato
CC    precedente del generatore random oppure il nuovo seme 
CC    (un numero intero qualsiasi)

CC    inoltre se si mette come "iflag" di partenza un numero che non sia
CC    0 o 1, prova a leggere la configurazione di spin dal file "lattice"
CC    che quindi deve gia` esistere


      call ranstart                 !! initialize random number generator
                                    !! per chiamare un numero random usare
                                    !! la funzione ran2(), vedere gli esempi
                                    !! riportati piu` avanti per capire come
                                    !! si usa

cc apertura file da dove leggere i parametri della simulazione
      open(1,file='input.txt',status='old')

cc apertura file sul quale scrivere le misure della magnetizzazione
      open(2,file='data.txt',status='unknown')

CC======================================
CC LETTURA PARAMETRI DELLA SIMULAZIONE
CC====================================== 
      read(1,*) iflag          !!partenza caldo(1)/freddo(0)/precedente(altro) 
      read(1,*) measures       !!numero di misure
      read(1,*) i_decorrel     !!updating fra una misura e altra
      read(1,*) extfield       !!valore del campo esterno (non in uso)
      read(1,*) mass           !!valore della massa
	read(1,*) i_term 
CC======================================
      mass2   = mass*mass
      mass2p4 = mass*mass + 4.0

CC======================================
CC OPERAZIONI PRELIMINARI
CC======================================
      CALL geometry()                     !!inizializza condizioni al bordo
      CALL initialize_lattice(iflag,seed) !!inizializza configurazione iniziale
CC======================================
      do it = 1,i_term
            call update_heatbath()
            call update_overrelax()
            call update_overrelax()
            call update_overrelax()
            call update_overrelax() 
      enddo


      do iter = 1,measures

CC   AGGIORNAMENTO CONFIGURAZIONE: i_decorrel spazzate di tutto il reticolo
         do idec = 1,i_decorrel
            call update_heatbath()
            call update_overrelax()
            call update_overrelax()
            call update_overrelax()
            call update_overrelax()
         enddo

CC   MISURA DELLE OSSERVABILI FISICHE
         call energy(xene_mass,xene_spat,xene_temp)

CC   SCRIVO LE MISURE SU FILE PER POI EFFETTUARE L'ANALISI
         write(2,*) xene_mass,xene_spat,xene_temp

      enddo !! measures
CC=========TERMINE SIMULAZIONE MONTE-CARLO===========

CC SALVO CONFIGURAZIONE E STATO GEN. RANDOM PER 
CC POTER EVENTUALMENTE RIPARTIRE
      open(3,file='lattice.txt',status='unknown')
      write(3,*) field
      close(3)      
      call ranfinish

CC==============================================
      stop
      end
CC===========F   I   N   E======================
CC===M  A   I   N      P  R  O  G  R  A  M======   





CC################################################################
CC      INIZIO SUBROUTINES
CC################################################################


c*****************************************************************
      subroutine geometry()    
c*****************************************************************
c per ogni coordinata definisco il passo in avanti o indietro
c con le opportune condizioni al bordo
c=================================================================
      parameter(nx = 100, nt = 50,nvol = nx*nt)
      parameter(nlatt = nx)   !! must be the maximum dimension size
      parameter(pigr = 3.141592654)
      real mass,extfield,mass2,mass2p4  
      common/various/mass,extfield,mass2,mass2p4 !! parametri della simulazione
      real field
      common/lattice/field(nx,nt) !! matrice con le variabili di spin
                                        !! che rappresenta lo stato del sistema
      common/move/npp(nlatt,2),nmm(nlatt,2) 


c     include "parameters_2d.f"
      
      !! le funzioni npp ed nmm sono costruite come dei vettori
      !! di interi, in modo da non essere ricalcolate ogni volta
      !! e rendere tutto piu` efficiente, prendono in input una 
      !! coordinata e restituiscono la coordinata in avanti o
      !! indietro, tenendo conto delle opportune condizioni

      !! le direzioni possono avere lunghezze diverse quindi mi servono
      !! 2 tabelle diverse


      do i = 1,nx
         npp(i,1) = i + 1
         nmm(i,1) = i - 1
      enddo
      npp(nx,1) = 1             !! RIAGGIUSTO IL CALCOLO AI BORDI PER TENERE
      nmm(1,1) = nx             !! CONTO DELLE CONDIZIONI AL BORDO PERIODICHE

      do i = 1,nt
         npp(i,2) = i + 1
         nmm(i,2) = i - 1
      enddo
      npp(nt,2) = 1             !! RIAGGIUSTO IL CALCOLO AI BORDI PER TENERE
      nmm(1,2) = nt             !! CONTO DELLE CONDIZIONI AL BORDO PERIODICHE

      return
      end
c=================================================================

c*****************************************************************
      subroutine initialize_lattice(iflag,seed)
c ASSEGNO LA CONFIGURAZIONE DI PARTENZA DELLA CATENA DI MARKOV
c=================================================================
      parameter(nx = 100,nt = 50,nvol = nx*nt)
      parameter(nlatt = nx)   !! must be the maximum dimension size
      parameter(pigr = 3.141592654)
      real mass,extfield,mass2,mass2p4  
      common/various/mass,extfield,mass2,mass2p4 !! parametri della simulazione
      real field
      common/lattice/field(nx,nt) !! matrice con le variabili di spin
                                        !! che rappresenta lo stato del sistema
      common/move/npp(nlatt,2),nmm(nlatt,2) 

c      include "parameters_2d.f"

      integer*8 seed

CC  PARTENZA A FREDDO (tutto il campo a zero come se fosse T = 0)
      if (iflag.eq.0) then
         do ix = 1,nx                        !! loop su tutti i siti
            do it = 1,nt                     !! del reticolo
               field(ix,it) = 0.0
            enddo
         enddo
CC  ... A CALDO ... (campo random fra -1 ed 1) 
      elseif (iflag.eq.1) then
         do ix = 1,nx                        !! loop su tutti i siti
            do it = 1,nt                     !! del reticolo
               x = ran2()               !! ran2() random fra 0 e 1
               field(ix,it) = 1.0 - 2.0*x
            enddo
         enddo

CC  ... O DA DOVE ERO RIMASTO L'ULTIMA VOLTA
      else
         open(9,file='lattice.txt',status='old')
         read(9,*) field
         read(9,*) seed  
         close(9)
      endif

      return
      end
c=================================================================


c*****************************************************************
      subroutine energy(xene_mass,xene_spat,xene_temp)
c================================================================= 
      parameter(nx = 100,nt = 50,nvol = nx*nt)
      parameter(nlatt = nx)   !! must be the maximum dimension size
      parameter(pigr = 3.141592654)
      real mass,extfield,mass2,mass2p4  
      common/various/mass,extfield,mass2,mass2p4 !! parametri della simulazione
      real field
      common/lattice/field(nx,nt) !! matrice con le variabili di spin
                                        !! che rappresenta lo stato del sistema
      common/move/npp(nlatt,2),nmm(nlatt,2) 
  


c      include "parameters_2d.f"  

      xene_mass = 0.0                             !! inizializzo variabile
      xene_spat = 0.0                             !! inizializzo variabile
      xene_temp = 0.0                             !! inizializzo variabile
      do ix = 1,nx                        !! loop su tutti i siti
         do it = 1,nt                     !! del reticolo

                 phi     = field(ix,it)
                 force_s = field(npp(ix,1),it)
                 force_t = field(ix,npp(it,2))

                 xene_mass = xene_mass + mass2*phi**2 
                 xene_spat = xene_spat - 2.0*phi*force_s + 2.0*phi**2
                 xene_temp = xene_temp - 2.0*phi*force_t + 2.0*phi**2
 
          enddo
      enddo

      xene_mass = xene_mass/float(nvol)   !! normalizzo -> densita` di energia !!! O1
      xene_spat = xene_spat/float(nvol)   !! normalizzo -> densita` di energia !!! O2 
      xene_temp = xene_temp/float(nvol)   !! normalizzo -> densita` di energia !!! O3

      return
      end
c=================================================================

c*****************************************************************
      subroutine update_heatbath()
c*****************************************************************
      parameter(nx = 100,nt = 50,nvol = nx*nt)
      parameter(nlatt = nx)   !! must be the maximum dimension size
      parameter(pigr = 3.141592654)
      real mass,extfield,mass2,mass2p4  
      common/various/mass,extfield,mass2,mass2p4 !! parametri della simulazione
      real field
      common/lattice/field(nx,nt) !! matrice con le variabili di spin
                                        !! che rappresenta lo stato del sistema
      common/move/npp(nlatt,2),nmm(nlatt,2) 



c      include "parameters_2d.f"

      do ix = 1,nx                        !! loop su tutti i siti
            do it = 1,nt               !! del reticolo

                 force = 0.0 
                 phi     = field(ix,it)
                 force   = force   + field(npp(ix,1),it)
                 force   = force   + field(nmm(ix,1),it)
                 force   = force   + field(ix,npp(it,2))
                 force   = force   + field(ix,nmm(it,2))


                 sigma2 = 1.0/mass2p4  !! variance of the gaussian distribution
                                       !! is 1/(m^2 + 4)
                 aver = force*sigma2   !! average of the gaussian distribution
                                       !! is force/(m^2 + 4)
                 
c                 write(*,*) mass2p4, aver, force
                 x = sqrt(sigma2)*sqrt(-2.0*log(ran2()))!! 
                 y = x*cos(2.0*pigr*ran2()) + aver !!  BOX MULLER ALGORITHM
                 field(ix,it) = y !!

         enddo
      enddo


      return
      end
c=================================================================

c*****************************************************************
      subroutine update_overrelax()
c*****************************************************************
      parameter(nx = 100,nt = 50,nvol = nx*nt)
      parameter(nlatt = nx)   !! must be the maximum dimension size
      parameter(pigr = 3.141592654)
      real mass,extfield,mass2,mass2p4  
      common/various/mass,extfield,mass2,mass2p4 !! parametri della simulazione
      real field
      common/lattice/field(nx,nt) !! matrice con le variabili di spin
                                        !! che rappresenta lo stato del sistema
      common/move/npp(nlatt,2),nmm(nlatt,2) 



c      include "parameters_2d.f"

      do ix = 1,nx                        !! loop su tutti i siti
         do it = 1,nt           !! del reticolo

                 force = 0.0 
                 phi     = field(ix,it)
                 force   = force   + field(npp(ix,1),it)
                 force   = force   + field(nmm(ix,1),it)
                 force   = force   + field(ix,npp(it,2))
                 force   = force   + field(ix,nmm(it,2))


                 aver = force/mass2p4  !! average of the gaussian distribution
                                       !! is force/(m^2 + 4)
                 
                 field(ix,it) = 2.0*aver - phi

         enddo
      enddo


      return
      end
c=================================================================

c============================================================================
c  RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &          ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &          ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,
     &          rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
c      save iv,iy,idum2
c      data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end

c=============================================================================
      subroutine ranstart
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed.txt', status='unknown')
      read(23,*) idum
      read(23,*,end=117) idum2
      do i=1,32
         read(23,*) iv(i)
      enddo
      read(23,*) iy
      close(23)
      goto 118                          !!takes account of the first start
 117  if(idum.ge.0) idum = -idum -1     !!
      close(23)
 118  continue                          !!

      return
      end

c=============================================================================
      subroutine ranfinish
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed.txt', status='unknown')
      write(23,*) idum
      write(23,*) idum2
      do i=1,32
         write(23,*) iv(i)
      enddo
      write(23,*) iy
      close(23)

      return
      end
c=============================================================================





















