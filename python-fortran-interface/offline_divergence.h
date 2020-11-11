        real*8, intent(in), dimension(nr)              :: dz    
        real*8, intent(in), dimension(nx, ny, nr)      :: volc
        real*8, intent(in), dimension(nx, ny, nr)      :: maskc
        real*8, intent(in), dimension(nx, ny, nr, nt)  :: fzon
        real*8, intent(in), dimension(nx, ny, nr, nt)  :: fmer
        real*8, intent(in), dimension(nx, ny, nr, nt)  :: fvert
        real*8, intent(in), dimension(nx, ny, nt)      :: klev
        real*8, intent(out), dimension(nx, ny, nt)     :: hfint
        real*8, intent(out), dimension(nx, ny, nt)     :: rfint