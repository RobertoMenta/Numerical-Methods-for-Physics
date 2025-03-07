      program SWITCH
      
      parameter (nlatt = 100)
      parameter (nlatt2 = 2*nlatt)
      real q(10000000),p(1000000),r(10000000)
      integer seed
      common/lattice/field(nlatt2)
      common/various/eta,d_metro,seed
      integer sigma
      common/sig/sigma

      open(1,file='path.txt',status='old')
      open(2,file='wave_function_right.txt',status='unknown')


      do i = 1,10000000
         read(1,*) r(i), q(i)
      enddo

      do k=1,1000000
         p(k) = q(k + 8000000)
      enddo

      do j=0,4999
            do jj=1,100
                write(2,*) p(2*j*100 + jj),p((2*j+1)*100 + jj)
            enddo
      enddo
      close(1)
      close(2)


      stop
      end
