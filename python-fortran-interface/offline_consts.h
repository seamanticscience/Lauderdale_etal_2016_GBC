C variable declarations etc for offline_routines_python.f
C Declare all parameters first
        real*8, parameter :: half=1.d0/2.d0
        real*8, parameter :: onesixth=1.d0/6.d0
        real*8, parameter :: thetamax = 1.d+20 
        real*8, parameter :: thetamin = 1.d-20 
        real*8, parameter :: eps = 1.d-20
        real*8, parameter :: rksign=-1.d0

C Array dimensions supplied as input arguments        
        integer*4, intent(in):: nx, ny, nr, nt

