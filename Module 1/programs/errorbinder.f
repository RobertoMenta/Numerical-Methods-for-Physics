c================================
            program errorbinder
c================================
            implicit real (a-h,o-z)
	        implicit integer (i-n)
            parameter (nfake = 500, measures = 50000)
	        real o(500)

                
            
            open(3,file='binboost.txt',status='unknown')

           
                
            zed = start
            zig = start

            do ik=1,nfake
	  	        read(3,*) o(ik)
	  	        zed = zed + o(ik)
		    enddo 
		    close(3)
    
	        zed = zed/float(nfake)
		

            do l=1,nfake
	 	        zig = zig + (o(l) - zed)**2
		    enddo
            zig = sqrt(zig/float(nfake))
	      

            write(*,*) 'c-binder =', zed, '+-', zig
   

            stop
            end 