       
	   
	   program costant
	   implicit real (a-h,o-z)
	   implicit integer (i-n)
       parameter (measures =1000000, nfake=1000)
c	   integer measures,nfake
	   integer*8 seed,seed_old
	   
	   
	   real r(1000000),s(1000000),t(1000000),xr(1000000),xs(1000000) 
	   real xu(1000000)
	   real xa(1000000)
	   real y(1000)
c		measures = 1000000
c		nfake= 1000
		
	   open(2,file='data.txt',status='unknown')
	   open(3,file ='data_tot.txt',status='unknown')
	    do i=1,measures 
		  read(2,*) r(i),s(i),t(i)
		enddo
        do j =1,measures
		   xr(j) = 0.0
		   xs(j) = 0.0
		   xu(j) = 0.0
		   xa(j) = 0.0
		enddo


		
		kh=1000
        	
         do i=1,nfake
		 xav = 0.0
	      do kj=0,999
	        ival=int(ran2()*kh)
		    ival=ival*1000
		    do j=1,1000                  
		      xr(kj*1000 +j)=r(ival + j)
			  xs(kj*1000 +j)=s(ival + j)
			  xu(kj*1000 +j)=t(ival + j)
			  
		      xa(kj*1000+j)=xr(kj*1000+j)+xs(kj*1000+j)-xu(kj*1000+j)
			  
		    enddo
		  enddo
		  do ik=1,measures
		     xav = xav + xa(ik)
		  enddo
		  xav = xav/float(measures)
		  write(3,*) xav
		 enddo
         close(2)
	  close(3)
      open(3, file='data_tot.txt',status = 'old')
	  
	  xmed = 0.0
	  x1= 0.0
	  var = 0.0
	  do iii=1,nfake
	     read(3,*) y(iii)
		 xmed = xmed + y(iii)
		 x1= x1 + y(iii)**2
	  enddo
	  xmed = xmed/float(nfake)
      x1 = x1/(float(nfake))
	  var= sqrt(x1 - (xmed)**2)
	  
	  write(*,*) xmed,var
      close(3)    
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