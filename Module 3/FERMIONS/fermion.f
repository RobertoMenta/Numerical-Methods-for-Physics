      program FERMIONI
      
      parameter (nlatt =100)
      parameter (nlatt2 = 2*nlatt)
      real lattice
      integer seed
      common/lattice/field(nlatt2)
      common/various/eta,d_metro,seed
      integer sigma
      common/sig/sigma

      open(1,file='input.txt',status='old')
      open(2,file='meas_out.txt',status='unknown')

CC====================================
CC LETTURA PARAMETRI DELLA SIMULAZIONE
CC==================================== 
      read(1,*) iflag          !!partenza caldo/freddo/precedente 
      read(1,*) measures       !!numero di misure
      read(1,*) i_decorrel     !!updating fra una misura e l altra
      read(1,*) i_term         !!passi di termalizzazione
c     read(1,*) d_metro        !!parametro del metropolis
      read(1,*) eta            !!valore del parametro eta = omega * a
CC=====================================
      d_metro = 2.0*sqrt(eta)

CC=====================================
CC OPERAZIONI PRELIMINARI
CC=====================================
      call ranstart                     !! initialize random number generator
      CALL initialize_lattice(iflag)    !! inizializza configurazione iniziale
      CALL geometry()                   !! inizializza condizioni al bordo
CC=====================================


CC=============================
CC TERMALIZZAZIONE
CC=============================
      do it = 1,i_term
         call update_metropolis()
      enddo
CC=============================


CC=============================================================
CC SESSIONE ALL'EQUILIBRIO CON MISURE
CC=============================================================
      do iter = 1,measures

CC   AGGIORNAMENTO CONFIGURAZIONE
         do idec = 1,i_decorrel
            call update_metropolis()
         enddo

         call measure()

      enddo 
      close(2)

c      open(7,file='path.txt', status='unknown')

c      do j=1,nlatt2
c         write(7,*) field(j),j
c      enddo
c      close(7)
      call tau()
CC=========TERMINE SIMULAZIONE MONTE-CARLO===========



CC==============================================
CC PRENDO L'ULTIMO STATO DEL GENERATORE RANDOM
CC==============================================
      call ranfinish

CC==============================================
CC SALVO CONFIGURAZIONE E STATO GEN. RANDOM PER POTER RIPARTIRE
CC==============================================
      open(3,file='lattice.txt',status='unknown')
      write(3,*) field
      close(3)      

CC==============================================
      stop
      end
CC===========F   I   N   E======================




CC      INIZIO DEFINIZIONE SUBROUTINES

c*****************************************************************
      subroutine geometry()    
c*****************************************************************
c DEFINISCO LA COORDINATA + e - 
c=================================================================
      parameter (nlatt = 100)
      parameter (nlatt2 = 2*nlatt)
      common/move/npp(nlatt2),nmm(nlatt2)
      integer sigma
      common/sig/sigma

      do i = 1,nlatt2
         npp(i) = i + 1
         nmm(i) = i - 1
      enddo

      npp(nlatt) = 1*(1-sigma) + (nlatt+1)*sigma    !!        CONDIZIONI AL BORDO
      npp(nlatt2) = (nlatt+1)*(1-sigma) + 1*sigma   !!        dipendenti da sigma
      nmm(1) = nlatt*(1-sigma) + nlatt2*sigma       !!        sigma  = 0 - no incrocio      
      nmm(nlatt+1) = nlatt2*(1-sigma) + nlatt*sigma !!        sigma  = 1 - incrocio

      return
      end
c=================================================================


c*****************************************************************
      subroutine initialize_lattice(iflag)
c*****************************************************************
c ASSEGNO LA CONFIGURAZIONE DI PARTENZA DELLA CATENA DI MARKOV
c=================================================================
      parameter (nlatt = 100)
      parameter (nlatt2 = 2*nlatt)
      integer sigma
      common/sig/sigma

      common/lattice/field(nlatt2)
      common/various/eta,d_metro,seed

      sigma = 0

CC  PARTENZA A FREDDO ...
      if (iflag.eq.0) then
         do i = 1,nlatt2                  !! loop su tutti i siti
               field(i) = 0.0
         enddo
CC  ... A CALDO ...
      elseif (iflag.eq.1) then
         do i = 1,nlatt2                  !! loop su tutti i siti
               x = 1.0 - 2.*ran2()       !! frand() random fra -1 e 1
                  field(i) = x
         enddo
CC  ... O DA DOVE ERO RIMASTO L'ULTIMA VOLTA
      else
         open(9,file='lattice.txt',status='old')
         read(9,*) field
c         read(9,*) seed  
         close(9)
      endif

      return
      end
c=================================================================


c*****************************************************************
      subroutine update_metropolis()
c*****************************************************************
 
      parameter (nlatt = 100)
      parameter (nlatt2 = 2*nlatt)

      common/lattice/field(nlatt2)
      common/move/npp(nlatt2),nmm(nlatt2)
      common/various/eta,d_metro,seed                  
      integer sigma
      common/sig/sigma

      c1 = 1./eta
      c2 = (1./eta + eta/2.)
      
C      open(7,file='path.txt', status='unknown')
      do i = 1,nlatt2            !! loop su tutti i siti, qui il sito 
                                !! non e` scelto a caso ma faccio una spazzata 
                                !! iterativa su tutti i siti, si puo` dimostrare
                                !! che va bene lo stesso per il bilancio, ma meno banale da provare 
         
         ip = npp(i)            !! calcolo le coordinate
         im = nmm(i)            !! dei due primi vicini
         
         force = field(ip) + field(im) !! costruisco la forza
         
         phi =  field(i)        !! phi = valore attuale del campo. 
         phi_prova = phi + 2.*d_metro*(0.5-ran2())                  

         p_rat = c1 * phi_prova * force - c2 * phi_prova**2
         p_rat = p_rat - c1 * phi * force + c2 * phi**2
         
         
         x = log(ran2())                      !! METRO-TEST! x = random (0,1)
                                              !! x < p_rat verifica anche caso 
         if (x.lt.p_rat) field(i) = phi_prova !! p_rat > 1 -> se si accetto
         write(7,*) i,field(i)
      enddo                     !!  chiudo il loop

         i = nlatt
         j = nlatt2
         ip = npp(i)
         jp = npp(j) 
         p_rat = -c1*(field(i) - field(j))*(field(ip) - field(jp)) 
         x = log(ran2())                      !! METRO-TEST! x = random (0,1)
                                              !! x < p_rat verifica anche caso 
         if (x.lt.p_rat) sigma = 1 - sigma
         call geometry()
c         write(*,*) sigma

      return
      end
c=================================================================


c*****************************************************************
      subroutine measure()
c*****************************************************************

      parameter (nlatt = 100)
      parameter (nlatt2 = 2*nlatt)
      common/lattice/field(nlatt2)
      common/move/npp(nlatt2),nmm(nlatt2)
      integer sigma
      common/sig/sigma

      obs = 0.0
      do i = 1,nlatt2      
        obs = obs + field(i)**2
      enddo          

      obs = obs/float(nlatt) 
      if (sigma.eq.0) then  
            write(2,*) obs, 1
      elseif (sigma.eq.1) then
            write(2,*) obs, -1
      endif

      return
      end
c=================================================================



c*****************************************************************
      subroutine tau()
c*****************************************************************
      parameter (nlatt = 100)
      parameter (nlatt2 = 2*nlatt)
      parameter (measures = 1000000, nfake = 500)

      real q(measures),s(measures),a(measures),b(measures),o(nfake)
      common/lattice/field(nlatt2)
      common/move/npp(nlatt2),nmm(nlatt2)
      common/various/eta,d_metro,seed                  
      integer sigma
      common/sig/sigma

      open(2,file='meas_out.txt',status='old')
      open(7,file='bootstrap.txt',status='unknown')
      xav = 0.0
      sign = 0.0
      r = 0.0
      t = 0.0
      do i=1,measures
            read(2,*) q(i),s(i)
           
      enddo
      close(2)

      do i=1,nfake 
            do j=0,999
                  ival = int(ran2()*1000)
                  ival = ival*1000
                  do k=1,1000
                        b(j*1000 + k) = q(k+ival)
                        a(j*1000 + k) = s(k+ival)
                        xav = xav + b(j*1000 + k)*a(j*1000 + k)
                        segn = segn + a(j*1000 + k)
                  enddo
            enddo
         
            xav = xav/float(measures)
            segn = segn/float(measures)
            rav = xav/segn
            write(7,*) rav
      enddo
      close(7)

      zed = 0.0
      zig = 0.0

      open(7,file='bootstrap.txt',status='old')
      do m=1,nfake
	  	read(7,*) o(m)
	  	zed = zed + o(m)
            zig = zig + o(m)**2
	enddo 
	close(7)
    
	zed = zed/float(nfake)

c      do l=1,nfake
c            zig = zig + o(l)**2
c	 	zig = zig + (o(l) - zed)**2
c	enddo
      
      zig = zig/float(nfake)
      zig = sqrt(zig - zed**2)

      
      write(*,*) zed,zig,segn

      return 
      end

c============================================================================         
            


 

      

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


