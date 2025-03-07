
      parameter(nx = 160,nt = 120,nvol = nx*nt)
      parameter(nlatt = nx)   
      parameter(pigr = 3.141592654)
      parameter(measures = 110000)
      real mass,extfield,mass2,mass2p4  
      common/various/mass,extfield,mass2,mass2p4 
      real field
      common/lattice/field(nx,nt) 
      common/move/npp(nlatt,2),nmm(nlatt,2) 
