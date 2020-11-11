C Arrays for fortran data outputs from computational routine
C Advective fluxes
        
        integer*4, intent(in):: advectionscheme
        
C Arrays for fortran data inputs to computational routine
        real*8, intent(in), dimension(nx, ny)          :: dy
        real*8, intent(in), dimension(nx, ny, nr)      :: ray
        real*8, intent(in), dimension(nx, ny, nr)      :: maskv
        real*8, intent(in), dimension(nx, ny, nr)      :: maskc
        real*8, intent(in), dimension(nx, ny, nr, nt)  :: tracer
        real*8, intent(in), dimension(nx, ny, nr, nt)  :: vvel
        real*8, intent(out), dimension(nx, ny, nr, nt) :: vflux
