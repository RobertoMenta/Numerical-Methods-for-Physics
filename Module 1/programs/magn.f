c==========================================
		program magnetizzazione media 
c==========================================
	
		implicit real (a-h,o-z)
		implicit integer (i-n)
		real q(50000)

    
    	open(3, file='data',status='unknown')

		start = 0.0
		nstat = 50000
	
		con = start
		sig = start 
		tau = start
		zed = start 

    	do i=1,nstat
	  		read(3,*) q(i)
	  		zed = zed + q(i)
	 	enddo 
	  	close(3)
	  	zed = zed/nstat

   		do n=1,nstat
	 		sig = sig + (q(n) - zed)**2
	 	enddo
	
		sig = sig/(float(nstat)*(float(nstat) - 1))
		

		do k=1,5000
			do j=1,nstat-k
				con = con + (q(j) - zed)*(q(j+k) - zed)
		 	enddo
			con = con * (1/(float(nstat)-k))
			tau = tau + con
			con = start
	 	enddo
		

		var = sig*(1+2*tau)
		var = sqrt(var)

		write(*,*) 'magnetizzazione', zed, '+-', var
	
		stop
		end	




















