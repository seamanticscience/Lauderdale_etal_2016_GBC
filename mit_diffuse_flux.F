#include "fintrf.h"
C======================================================================
C     MIT_DIFFUSE_FLUX calculates explicit scalar diffusive transports
C     [uflux,vflux,wflux]=mit_diffuse_flux( dx ,dy ,dz ,
C                                           umask,vmask,cmask.
C                                           rax,ray,rac,
C                                           kux,kuz,kvy,kvz,
C                                           kwx,kwy,scalar);
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
      mwSize mxGetNumberOfDimensions
      mwPointer mxGetDimensions
      mwPointer mxGetPr

C     Arguments for mxCreateNumericArray
      integer*4 classid, complexflag
      mwSize ndim, dims(4)
    
C     Pointers to input/output mxArrays:
      mwPointer dx, dy, dz
      mwPointer rax, ray, rac
      mwPointer masku, maskv, maskc
      mwPointer kux,kuz,kvy,kvz,kwx,kwy
      mwPointer scalar 
      mwPointer uflux, vflux, wflux

CC    For writing error to the screen
C     mexErrMsgIdAndTxt causes crash
      character*250 line
      integer*4 k
      integer*4 mexPrintf

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 16) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nInput     ',
C     &                           'Fourteen inputs required. ')
         write(line,*)'ERROR: Sixteen inputs required.'
         k=mexPrintf(line//achar(13)) 
         return
      elseif(nlhs .gt. 3) then
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
      dx   = mxGetPr(prhs(1))
      dy   = mxGetPr(prhs(2))
      dz   = mxGetPr(prhs(3))
      masku= mxGetPr(prhs(4))
      maskv= mxGetPr(prhs(5))
      maskc= mxGetPr(prhs(6))
      rax  = mxGetPr(prhs(7))
      ray  = mxGetPr(prhs(8))
      rac  = mxGetPr(prhs(9))
      kux  = mxGetPr(prhs(10))
      kuz  = mxGetPr(prhs(11))
      kvy  = mxGetPr(prhs(12))
      kvz  = mxGetPr(prhs(13))
      kwx  = mxGetPr(prhs(14))
      kwy  = mxGetPr(prhs(15))
      scalar = mxGetPr(prhs(16))

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

CC----------------DIFFUSION IN X DIRECTION----------------

         call mit_diff_x(dims(1),dims(2),dims(3),dims(4),
     &           %val(maskc),%val(dz),%val(dx),%val(rax),
     &           %val(kux),%val(kuz),%val(scalar),%val(uflux))

CC----------------DIFFUSION IN Y DIRECTION----------------

         call mit_diff_y(dims(1),dims(2),dims(3),dims(4),
     &           %val(maskc),%val(dz),%val(dy),%val(ray),
     &           %val(kvy),%val(kvz),%val(scalar),%val(vflux))

CC----------------DIFFUSION IN Z DIRECTION----------------

         call mit_diff_z(dims(1),dims(2),dims(3),dims(4),
     &           %val(masku),%val(maskv),%val(maskc),
     &           %val(dx),%val(dy),%val(rac),%val(kwx),
     &           %val(kwy),%val(scalar),%val(wflux))
     
      return
      end