#include "fintrf.h"
C======================================================================
C     MIT_ADVECT_FLUX calls advection schemes for scalar transports
C     [uflux,vflux,wflux]=cadvect(timestep,advscheme,
C                                  dx,dy,dz,
C                                  rax,ray,rac,
C                                  umask,vmask,cmask,
C                                  uvel,vvel,wvel,
C                                  scalar);
C     
C======================================================================
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer*4 nlhs, nrhs

C     Function declarations:
      mwPointer mxCreateNumericArray
      integer*4 mxClassIDFromClassName
      integer mxIsDouble, mxIsNumeric
      mwPointer mxGetM, mxGetN
      mwSize mxGetNumberOfDimensions
      mwPointer mxGetDimensions
      mwPointer mxGetPr

C     Arguments for mxCreateNumericArray
      integer*4 classid, complexflag
      mwSize ndim, dims(4)
    
C     Pointers to input/output mxArrays:
      mwPointer deltat, advscheme
      mwPointer dx, dy, dz
      mwPointer rax, ray, rac
      mwPointer masku,maskv,maskc
      mwPointer uvel, vvel, wvel
      mwPointer scalar 
      mwPointer uflux, vflux, wflux

C     Array information:
      mwPointer mrows, ncols
      mwSize size

CC    For writing error to the screen
C     mexErrMsgIdAndTxt causes crash
      character*250 line
      integer*4 k
      integer*4 mexPrintf

CC     Arguments for computational routine:
      real*8 iAdvSch
      
C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 15) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nInput     ',
C     &                           'Thirteen inputs required. ')
         write(line,*)'ERROR: Fifteen inputs required.'
         k=mexPrintf(line//achar(13)) 
         return
      elseif(nlhs .gt. 4) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nOutput    ',
C     &                           'Too many output arguments.')
         write(line,*)'ERROR: Too many output arguments.'
         k=mexPrintf(line//achar(13)) 
         return
      endif

C     Validate inputs
C     Check that the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:NonNumeric',
C     &                           'Input must be a number.  ')
         write(line,*)'ERROR: Input must be numeric.'
         k=mexPrintf(line//achar(13)) 
         return     
      endif
      
C     Check data type of input argument.
      if (mxIsDouble(prhs(1)) .eq. 0) then
C         call mexErrMsgTxt('Input argument must be of type double.')
         write(line,*)'ERROR: Input must be double precision.'
         k=mexPrintf(line//achar(13)) 
         return     
      endif
      
C     Get the size of the input arrays using scalar as template.
C     Use mxGetDimensions to get dimensions into dims
      ndim = mxGetNumberOfDimensions(prhs(nrhs))
      call mxCopyPtrToInteger4(mxGetDimensions(prhs(nrhs)),  
     &           dims, mxGetNumberOfDimensions(prhs(nrhs)))
      
C      write(line,*)'ndims: ',ndim
C      k=mexPrintf(line//achar(13)) 
C      write(line,*)'dimensions: ',dims(1),'x',dims(2),
C     &                        'x',dims(3),'x',dims(4)
C      k=mexPrintf(line//achar(13)) 
          
C     Create Fortran arrays from the input arguments.
      deltat = mxGetPr(prhs(1))
      advscheme = mxGetPr(prhs(2))
      dx = mxGetPr(prhs(3))
      dy = mxGetPr(prhs(4))
      dz = mxGetPr(prhs(5))
      rax = mxGetPr(prhs(6))
      ray = mxGetPr(prhs(7))
      rac = mxGetPr(prhs(8))
      masku= mxGetPr(prhs(9))
      maskv= mxGetPr(prhs(10))
      maskc= mxGetPr(prhs(11))
      uvel = mxGetPr(prhs(12))
      vvel = mxGetPr(prhs(13))
      wvel = mxGetPr(prhs(14))
      scalar = mxGetPr(prhs(15))

C     Copy advection scheme to a variable
      mrows = mxGetM(prhs(2))
      ncols = mxGetN(prhs(2))
      size = mrows*ncols
      call mxCopyPtrToReal8(advscheme, iAdvSch, size)

C     Output will be uflux,vflux,wflux
C     Create matrix for the return argument.
      classid = mxClassIDFromClassName('double')
      complexflag = 0

      plhs(1) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)
      plhs(2) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)
      plhs(3) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)

      uflux = mxGetPr(plhs(1))
      vflux = mxGetPr(plhs(2))
      wflux = mxGetPr(plhs(3))
CC     Call the computational subroutines.

CC----------------ADVECTION IN X DIRECTION----------------
      if (iAdvSch .eq. 1) then
C        first order upwind      
         call gad_upwind_x(dims(1),dims(2),dims(3),dims(4),
     &           %val(rax),%val(maskc),%val(masku),
     &           %val(uvel),%val(scalar),%val(uflux))
      elseif (iAdvSch .eq. 2) then
C        centre second order      
         call gad_adv_x(dims(1),dims(2),dims(3),dims(4),
     &           %val(rax),%val(maskc),%val(masku),
     &           %val(uvel),%val(scalar),%val(uflux))
      elseif (iAdvSch .eq. 3) then
C        third order upwind biased
         call gad_u3_adv_x(dims(1),dims(2),dims(3),dims(4),
     &           %val(rax),%val(maskc),%val(masku),%val(deltaT),
     &           %val(dx),%val(uvel),%val(scalar),%val(uflux))
      elseif (iAdvSch .eq. 4) then
C        centre fourth order
         call gad_c4_adv_x(dims(1),dims(2),dims(3),dims(4),
     &           %val(rax),%val(maskc),%val(masku),
     &           %val(uvel),%val(scalar),%val(uflux))
      elseif (iAdvSch .eq. 7) then
C     7th Order One Step method with Monotonicity Preserving Limiter
         call gad_os7mp_adv_x(dims(1),dims(2),dims(3),dims(4),
     &           %val(rax),%val(maskc),%val(masku),%val(deltaT),
     &           %val(dx),%val(uvel),%val(scalar),%val(uflux))
      elseif (iAdvSch .eq. 20 .or. iAdvSch .eq. 21) then
C        second order direct-space-and-time, first upwind
         call gad_dst2u1_adv_x(dims(1),dims(2),dims(3),dims(4),
     &           %val(advscheme),    
     &           %val(rax),%val(maskc),%val(masku),%val(deltaT),
     &           %val(dx),%val(uvel),%val(scalar),%val(uflux))
      elseif (iAdvSch .eq. 30) then
C        third order direct-space-and-time
         call gad_dst3_adv_x(dims(1),dims(2),dims(3),dims(4),
     &           %val(rax),%val(maskc),%val(masku),%val(deltaT),
     &           %val(dx),%val(uvel),%val(scalar),%val(uflux))
      elseif (iAdvSch .eq. 33) then
C        third order direct-space-and-time (flux limited)    
         call gad_dst3fl_adv_x(dims(1),dims(2),dims(3),dims(4),
     &           %val(rax),%val(maskc),%val(masku),%val(deltaT),
     &           %val(dx),%val(uvel),%val(scalar),%val(uflux))
      elseif (iAdvSch .eq. 77) then
C        centre second order (flux liimited)    
         call gad_fluxlimit_adv_x(dims(1),dims(2),dims(3),dims(4),
     &           %val(rax),%val(maskc),%val(masku),%val(deltaT),
     &           %val(dx),%val(uvel),%val(scalar),%val(uflux))
      else
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nInput',
C     &                           'iAdvSch = 2,3,30,33,77.')
         write(line,*)'ERROR: iAdvSch = 1,2,3,30,33,77.'
         k=mexPrintf(line//achar(13)) 
         return     
      endif
     
CC----------------ADVECTION IN Y DIRECTION----------------      
      if (iAdvSch .eq. 1) then
C        first order upwind     
         call gad_upwind_y(dims(1),dims(2),dims(3),dims(4),
     &           %val(ray),%val(maskc),%val(maskv),
     &           %val(vvel),%val(scalar),%val(vflux))
      elseif (iAdvSch .eq. 2) then
C        centre second order      
         call gad_adv_y(dims(1),dims(2),dims(3),dims(4),
     &           %val(ray),%val(maskc),%val(maskv),
     &           %val(vvel),%val(scalar),%val(vflux))
      elseif (iAdvSch .eq. 3) then
C        third order upwind biased
         call gad_u3_adv_y(dims(1),dims(2),dims(3),dims(4),
     &           %val(ray),%val(maskc),%val(maskv),%val(deltaT),
     &           %val(dy),%val(vvel),%val(scalar),%val(vflux))
      elseif (iAdvSch .eq. 4) then
C        centre fourth order
         call gad_c4_adv_y(dims(1),dims(2),dims(3),dims(4),
     &           %val(ray),%val(maskc),%val(maskv),
     &           %val(vvel),%val(scalar),%val(vflux))
      elseif (iAdvSch .eq. 7) then
C     7th Order One Step method with Monotonicity Preserving Limiter
         call gad_os7mp_adv_y(dims(1),dims(2),dims(3),dims(4),
     &           %val(ray),%val(maskc),%val(maskv),%val(deltaT),
     &           %val(dy),%val(vvel),%val(scalar),%val(vflux))
      elseif (iAdvSch .eq. 20 .or. iAdvSch .eq. 21) then
C        second order direct-space-and-time, first upwind
         call gad_dst2u1_adv_y(dims(1),dims(2),dims(3),dims(4),
     &           %val(advscheme),    
     &           %val(ray),%val(maskc),%val(maskv),%val(deltaT),
     &           %val(dy),%val(vvel),%val(scalar),%val(vflux))
      elseif (iAdvSch .eq. 30) then
C        third order direct-space-and-time
         call gad_dst3_adv_y(dims(1),dims(2),dims(3),dims(4),
     &           %val(ray),%val(maskc),%val(maskv),%val(deltaT),
     &           %val(dy),%val(vvel),%val(scalar),%val(vflux))
      elseif (iAdvSch .eq. 33) then
C        third order direct-space-and-time (flux limited)    
         call gad_dst3fl_adv_y(dims(1),dims(2),dims(3),dims(4),
     &           %val(ray),%val(maskc),%val(maskv),%val(deltaT),
     &           %val(dy),%val(vvel),%val(scalar),%val(vflux))
      elseif (iAdvSch .eq. 77) then
C        centre second order (flux liimited)    
         call gad_fluxlimit_adv_y(dims(1),dims(2),dims(3),dims(4),
     &           %val(ray),%val(maskc),%val(maskv),%val(deltaT),
     &           %val(dy),%val(vvel),%val(scalar),%val(vflux))
      else
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nInput',
C     &                           'iAdvSch = 2,3,30,33,77.')
         write(line,*)'ERROR: iAdvSch = 1,2,3,30,33,77.'
         k=mexPrintf(line//achar(13)) 
         return     
      endif

CC----------------ADVECTION IN Z DIRECTION----------------
      if (iAdvSch .eq. 1) then
C        first order upwind      
         call gad_upwind_z(dims(1),dims(2),dims(3),dims(4),
     &           %val(rac),%val(maskc),
     &           %val(wvel),%val(scalar),%val(wflux))
      elseif (iAdvSch .eq. 2) then
C        centre second order      
         call gad_adv_z(dims(1),dims(2),dims(3),dims(4),
     &           %val(rac),%val(maskc),
     &           %val(wvel),%val(scalar),%val(wflux))
      elseif (iAdvSch .eq. 3) then
C        third order upwind biased
         call gad_u3_adv_z(dims(1),dims(2),dims(3),dims(4),
     &           %val(rac),%val(maskc),%val(deltaT),%val(dz),
     &           %val(wvel),%val(scalar),%val(wflux))
      elseif (iAdvSch .eq. 4) then
C        centre fourth order
         call gad_c4_adv_z(dims(1),dims(2),dims(3),dims(4),
     &           %val(rac),%val(maskc),
     &           %val(wvel),%val(scalar),%val(wflux))
      elseif (iAdvSch .eq. 7) then
C     7th Order One Step method with Monotonicity Preserving Limiter
         call gad_os7mp_adv_z(dims(1),dims(2),dims(3),dims(4),
     &           %val(rac),%val(maskc),%val(deltaT),%val(dz),
     &           %val(wvel),%val(scalar),%val(wflux))
      elseif (iAdvSch .eq. 20 .or. iAdvSch .eq. 21) then
C        second order direct-space-and-time, 1st upwind
         call gad_dst2u1_adv_z(dims(1),dims(2),dims(3),dims(4),
     &           %val(advscheme),    
     &           %val(rac),%val(maskc),%val(deltaT),%val(dz),
     &           %val(wvel),%val(scalar),%val(wflux))
      elseif (iAdvSch .eq. 30) then
C        third order direct-space-and-time
         call gad_dst3_adv_z(dims(1),dims(2),dims(3),dims(4),
     &           %val(rac),%val(maskc),%val(deltaT),%val(dz),
     &           %val(wvel),%val(scalar),%val(wflux))
      elseif (iAdvSch .eq. 33) then
C        third order direct-space-and-time (flux liimited)    
         call gad_dst3fl_adv_z(dims(1),dims(2),dims(3),dims(4),
     &           %val(rac),%val(maskc),%val(deltaT),%val(dz),
     &           %val(wvel),%val(scalar),%val(wflux))
      elseif (iAdvSch .eq. 77) then
C        centre second order (flux liimited)    
         call gad_fluxlimit_adv_z(dims(1),dims(2),dims(3),dims(4),
     &           %val(rac),%val(maskc),%val(deltaT),%val(dz),
     &           %val(wvel),%val(scalar),%val(wflux))
      else
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nInput',
C     &                           'iAdvSch = 2,3,30,33,77.')
         write(line,*)'ERROR: iAdvSch = 1,2,3,30,33,77.'
         k=mexPrintf(line//achar(13)) 
         return     
      endif
      
      return
      end