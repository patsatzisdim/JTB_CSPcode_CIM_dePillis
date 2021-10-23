cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            THIS FILE CONTAINS THE SUBROUTINES THAT DEFINE THE PROBLEM
c
c                     REACTION RATES (rates)
c                     STOICHIOMETRY (stoic)     
c                     GRADIENT OF REACTION RATES (gradR)
c                     TIME DERIVATIVE OF JACOBIAN (Djac_dt)                  
c                 
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd   SUBROUT used to generate the Rates
      subroutine rates(n,t,yy,k,RR)
      implicit none 
      integer, intent(in) :: n,k
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: RR(k)
      double precision :: y1,y2,y3,y4,y5,y6
      include 'paramet.i'

cpd   assign each variable to a vector
      y1=yy(1)
      y2=yy(2)
      y3=yy(3)
      y4=yy(4)
c      y5=yy(5)
c      y6=yy(6)   
cpd
cpd   Reaction Rates from mathematica
      RR(1)=ai*y1*(1.0d0-bi*y1)
      RR(2)=ei*y4         
      RR(3)=alp        
      RR(4)=fi*y2
      RR(5)=mi*y3
      RR(6)=bet*y4        
c      RR(7)=gam*y5
c      RR(8)=muI*y6
      RR(7)=ci*y1*y2
      RR(8)=di*y1*(1.0d0-si/(si+(y3/y1)**li))
c      RR(9)=KT*(1.0d0-exp(-y5))*y1
c      RR(10)=KN*(1.0d0-exp(-y5))*y2
c      RR(11)=KL*(1.0d0-exp(-y5))*y3
c      RR(12)=KC*(1.0d0-exp(-y5))*y4
c      RR(13)=(piI*y3*y6)/(giI+y6)
      RR(9)=(gi*y1**2*y2)/(hi+y1**2)        
      RR(10)=y3*(di**2*ji*y1**2*(y3/y1)**(2*li))/((si + (y3/y1)**li)**2*(ki + (di**2*y1**2*(y3/y1)**(2*li))/(si + (y3/y1)**li)**2))
      RR(11)=ri1*y1*y2
      RR(12)=ri2*y1*y4
      RR(13)=pi*y1*y2
      RR(14)=qi*y1*y3
      RR(15)=ui*y2*y3**2                
cpd          
      return
      end
cpd
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
c
c
c
c         Stoichiometry
c
c
c
c
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd   SUBROUT used for stoichiometry
      subroutine stoic(n,k,st)
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(out) :: st(n,k)
      include 'paramet.i'
cpd
cpd   initialize first      
      st(:,:)=0.0d0
cpd
cpd   values from mathematica
      st(1,1)=1.0d0
      st(1,7)=-1.0d0
      st(1,8)=-1.0d0
      
      st(2,2)=1.0d0
      st(2,4)=-1.0d0
      st(2,9)=1.0d0
      st(2,13)=-1.0d0
      
      st(3,5)=-1.0d0
      st(3,10)=1.0d0
      st(3,11)=1.0d0
      st(3,12)=1.0d0
      st(3,14)=-1.0d0
      st(3,15)=-1.0d0
      
      st(4,3)=1.0d0
      st(4,6)=-1.0d0

cpd      
      return
      end      
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
c
c
c
c
c
c
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd   SUBROUT used to calulate grad(R) analytically from mathematica
      subroutine gradR(n,t,yy,k,gR) 
      implicit none
      integer, intent(in) :: n,k
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: gR(k,n)
      double precision :: y1,y2,y3,y4,y5,y6,y7
      include 'paramet.i'
cpd   assign each variable to a vector
      y1=yy(1)
      y2=yy(2)
      y3=yy(3)
      y4=yy(4)
      y5=yy(5)
      y6=yy(6)
cpd
cpd   initialize first
      gR(:,:)=0.0d0

cpd
cpd   values grad(R) from mathematica  
c     
      gR(1,1)=ai*(1.0d0-2.0d0*bi*y1)
      gR(2,4)=ei
      gR(4,2)=fi
      gR(5,3)=mi
      gR(6,4)=bet
c      gR(7,5)=gam
c      gR(8,6)=muI
      gR(7,1)=ci*y2      
      gR(7,2)=ci*y1
      gR(8,1)=(di*(y3/y1)**li*(si-li*si+(y3/y1)**li))/(si+(y3/y1)**li)**2 
      gR(8,3)=(di*li*si*(y3/y1)**(-1.0d0+li))/(si+(y3/y1)**li)**2
c      gR(11,1)=KT*(1.0d0-exp(-y5))
c      gR(11,5)=KT*y1*exp(-y5)
c      gR(12,2)=KN*(1.0d0-exp(-y5))      
c      gR(12,5)=KN*y2*exp(-y5)
c      gR(13,3)=KL*(1.0d0-exp(-y5)) 
c      gR(13,5)=KL*y3*exp(-y5)
c      gR(14,4)=KC*(1.0d0-exp(-y5))
c      gR(14,5)=KC*y4*exp(-y5)
c      gR(15,3)=(piI*y6)/(giI+y6)      
c      gR(15,6)=(giI*pi*y3)/(giI+y6)**2
      gR(9,1)=(2.0d0*gi*hi*y1*y2)/(hi+y1**2)**2 
      gR(9,2)=(gi*y1**2)/(hi+y1**2)
      gR(10,1)=(2.0d0*di**2*ji*ki*y1*(y3/y1)**(2*li)*(si + (y3/y1)**li)*(si - li*si + (y3/y1)**li))/
     -  (di**2*y1**2*(y3/y1)**(2*li) + ki*(si + (y3/y1)**li)**2)**2
      gR(10,3)=(2.0d0*di**2*ji*ki*li*si*y1*(y3/y1)**(-1.0d0 + 2*li)*(si + (y3/y1)**li))/
     -  (di**2*y1**2*(y3/y1)**(2*li) + ki*(si + (y3/y1)**li)**2)**2
      gR(11,1)=ri1*y2     
      gR(11,2)=ri1*y1
      gR(12,1)=ri2*y4 
      gR(12,4)=ri2*y1
      gR(13,1)=pi*y2
      gR(13,2)=pi*y1
      gR(14,1)=qi*y3      
      gR(14,3)=qi*y1
      gR(15,2)=ui*y3**2 
      gR(15,3)=2.0d0*ui*y2*y3
c
      return
      end


ccpd-----------------------------------------------
ccpd
ccpd------------------------------------------------------------------
ccpd   SUBROUT to calculate dJac/dt
      subroutine Djac_dt(n,t,yy,dJacdt)
      integer, intent(in) :: n
      double precision, intent(in) :: t,yy(n)
      double precision, intent(out) :: dJacdt(n,n)
      double precision :: y1,y2,y3,y4,y5,y6,y7,ydot(n)
      double precision :: dJ1(n,n),dJ2(n,n),dJ3(n,n),dJ4(n,n),dJ5(n,n),dJ6(n,n),dJ7(n,n)
      include 'paramet.i'
c
c     Right hand side
      call FEX(n,t,yy,ydot)
c      
cpd   assign each variable to a vector
      y1=yy(1)
      y2=yy(2)
      y3=yy(3)
      y4=yy(4)
      y5=yy(5)
      y6=yy(6)
c
c
c
      dJ1(:,:)=0.0d0
      dJ1(:,:)=dJ1(:,:)*ydot(1)
cc
cc
      dJ2(:,:)=0.0d0
      dJ2(:,:)=dJ2(:,:)*ydot(2)
cc      
cc
      dJ3(:,:)=0.0d0
      dJ3(:,:)=dJ3(:,:)*ydot(3)
cc
cc     
      dJ4(:,:)=0.0d0 
      dJ4(:,:)=dJ4(:,:)*ydot(4)
cc
cc     
      dJ5(:,:)=0.0d0 
      dJ5(:,:)=dJ5(:,:)*ydot(5)
cc
cc     
      dJ6(:,:)=0.0d0 
      dJ6(:,:)=dJ6(:,:)*ydot(6)
c            
c
      dJacdt(:,:)=dJ1(:,:)+dJ2(:,:)+dJ3(:,:)+dJ4(:,:)+dJ5(:,:)+dJ6(:,:)  
c
c           
      return
      end
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------
cpd------------------------------------------------------------------