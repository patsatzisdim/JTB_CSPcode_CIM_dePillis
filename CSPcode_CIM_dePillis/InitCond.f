cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            THIS FILE CONTAINS THE INITIAL CONDITIONS FOR THE PROBLEM TO RUN
c
c
c            WATCH OUT! IT GIVES YOU THE OPPORTUNITY TO READ THE INITIAL CONDITIONS FROM
C                       THE END OF THE CURRENT SOLUTION FILE. THIS HELPS FOR MULTPLE CALLS
C                       FROM THE SCRIPT FILE!!!!!!
C
C            (COMMENTED OUT CURRENTLY)
c

      subroutine InitCond(n,t,y)
      implicit none
      include 'paramet.i'
      integer :: i
      integer, intent(in) :: n
      double precision, intent(out) :: y(n)
      double precision, intent(inout) :: t
      integer :: IOState,iflag


c      iflag=2
      do i=1,n
       y(i)=1.0d0
      enddo

      y(1)=Tc0                 
      y(2)=NK0              
      y(3)=L0
      y(4)=C0    
      y(5)=0.0d0
      y(6)=0.0d0                 

      return 
      end
cpd------------------------------------------------------------------
