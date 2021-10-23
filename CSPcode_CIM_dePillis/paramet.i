cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            THIS FILE CONTAINS IN FORTRAN FORM ALL THE PARAMETERS
c                 THAT ARE GOING TO BE USED FOR THE CODE TO RUN
c
c            EXCEPT FROM PROBLEM PARAMETERS IT ALSO HAS TIME AND OTHER STUFF
c
      DOUBLE PRECISION ai,bi,ei,fi,mi,alp,bet,gam,muI
      DOUBLE PRECISION ci,di,si,li
      DOUBLE PRECISION KT,KL,KN,KC
      DOUBLE PRECISION piI,giI
      DOUBLE PRECISION gi,hi,ji,ki,ri1,ri2
      DOUBLE PRECISION qi,pi,ui
      DOUBLE PRECISION Tc0,NK0,L0,C0
      DOUBLE PRECISION tend
C==============================================================C
C===============    Patient 9    Parameters   =================C
C==============================================================C

C==============================================================C
C=============== Growth and Death Parameters  =================C
C==============================================================C
      parameter(ai=4.31d-1)         ! 1/day       Normal value 4.31d-1
      parameter(bi=1.02d-9)         ! 1/cells
      parameter(ei=2.08d-7)         ! 1/day         
      parameter(fi=4.12d-2)         ! 1/day 
      parameter(mi=2.04d-1)         ! 1/day
      parameter(alp=7.5d+8)         ! cells/day
      parameter(bet=1.2d-2)         ! 1/day
      parameter(gam=9.0d-1)         ! 1/day  
      parameter(muI=1.0d+1)         ! 1/day 
C==============================================================C
C===============  cell-cell kill parameters   =================C
C==============================================================C
      parameter(ci=6.41d-11)        ! 1/(day cells)      
      parameter(di=2.34d0)          ! 1/day 
      parameter(si=8.39d-2)         !
      parameter(li=2.09d0)          !
C==============================================================C
C===============  cell-chemo kill parameters  =================C
C==============================================================C             
      parameter(KT=9.0d-1)          ! 1/day
      parameter(KL=6.0d-1)          ! 1/day
      parameter(KN=6.0d-1)          ! 1/day
      parameter(KC=6.0d-1)          ! 1/day
C==============================================================C
C=============== cell-immuno kill parameters  =================C
C==============================================================C 
      parameter(piI=1.25d-1)        ! 1/day
      parameter(giI=2.0d+7)         ! cells^2
C==============================================================C
C===============    recruitment parameters    =================C
C==============================================================C
      parameter(gi=1.25d-2)         ! 1/day       
      parameter(hi=2.02d+7)         ! cells^2    
      parameter(ji=2.49d-2)         ! 1/day       
      parameter(ki=3.66d+7)         ! cells^2       
      parameter(ri1=1.1d-7)         ! 1/(day cells)       
      parameter(ri2=6.5d-11*0.8d0)        ! 1/(day cells) 
C==============================================================C
C===============   inactivation parameters    =================C
C==============================================================C
      parameter(qi=1.42d-6)         ! 1/(day cells)
      parameter(pi=3.42d-6)         ! 1/(day cells)                  
      parameter(ui=3.0d-10)         ! 1/(day cells^2)        
C==============================================================C
C=============== Initial (non-zero) conditions ================C
C==============================================================C
      parameter(Tc0=3.19392d+5+1.0d0)         ! cells       
      parameter(NK0=1.0d+3)         ! cells        
      parameter(L0=1.0d+1)          ! cells 
      parameter(C0=6.0d+8)          ! cells

C=============== Tumor Persistence ================C
c      parameter(Tc0=1.0d+6)         ! cells       
c      parameter(NK0=1.0d+3)         ! cells        
c      parameter(L0=1.0d+1)          ! cells 
c      parameter(C0=6.0d+8)          ! cells
C=============== Tumor Remission  ================C
c      parameter(Tc0=1.0d+6)         ! cells       
c      parameter(NK0=1.0d+5)         ! cells        
c      parameter(L0=1.0d+2)          ! cells 
c      parameter(C0=6.0d+10)          ! cells            
C==============================================================C
C================    Simulation paremeters    =================C
C==============================================================C

      parameter(tend=2.0d+3)        ! day

