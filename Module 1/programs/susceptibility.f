c====================
        program susceptibility 
c===================

        implicit real (a-h,o-z)
	    implicit integer (i-n)
		parameter(nfake=500)
	    real q(500), p(500)


        open(4,file='chi.txt',status='unknown')  
	    start = 0.0

	    sig = start 
	    zed = start 
		sigg = start 
	    zedd = start
        Dchi = start

        do i=1,nfake
	  	    read(4,*) q(i), p(i)
	  	    zed = zed + q(i)
			zedd = zedd + p(i)
	    enddo 
		close(4)
    
	    zed = zed/float(nfake)
		zedd = zedd/float(nfake)

   	    do n=1,nfake
	 	    sig = sig + (q(n) - zed)**2
			sigg = sigg + (p(n) - zedd)**2
	    enddo
        sig = sqrt(sig/float(nfake))
	    sigg = sqrt(sigg/float(nfake))

		write(*,*) 'Suscettibilità', zed, '+-', sig
		write(*,*) 'Calore Specifico', zedd, '+-', sigg
        end	




















