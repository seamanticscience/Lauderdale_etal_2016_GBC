#include "fintrf.h"
C======================================================================
C     CADVECT calls advection schemes to calculate scalar transports
C     [wflux,itertrac]=mit_difimpl_flux(timestep,rac,
C                               dzf,dzc,
C                               mask,
C                               diffkz,
C                               scalar);
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
      mwSize ndim, dims(4), nd, dimsd(5)
    
C     Pointers to input/output mxArrays:
      mwPointer deltat
      mwPointer iterret
      mwPointer dzf, dzc, rac
      mwPointer mask
      mwPointer diffkz
      mwPointer scalar 
      mwPointer wflux

C     Array information:
      mwPointer mrows, ncols
      mwSize size

CC    For writing error to the screen
C     mexErrMsgIdAndTxt causes crash
      character*250 line
      integer*4 k
      integer*4 mexPrintf
      integer*4 errCode
      real*8 iter
C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 8) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nInput     ',
C     &                           'Fourteen inputs required. ')
         write(line,*)'ERROR: Eight inputs required.'
         k=mexPrintf(line//achar(13)) 
         return
      elseif(nlhs .gt. 1) then
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
      iterret= mxGetPr(prhs(2))
      rac    = mxGetPr(prhs(3))
      dzf    = mxGetPr(prhs(4))
      dzc    = mxGetPr(prhs(5))
      mask   = mxGetPr(prhs(6))
      diffkz = mxGetPr(prhs(7))
      scalar = mxGetPr(prhs(8))
      
C     Copy iteration number to a variable
      mrows = mxGetM(prhs(2))
      ncols = mxGetN(prhs(2))
      size = mrows*ncols
      call mxCopyPtrToReal8(iterret, iter, size)
      
C     Output will be wflux
C     Create matrix for the return argument.
      classid = mxClassIDFromClassName('double')
      complexflag = 0

      plhs(1) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)
     
C     nd = 5
C     dimsd(1)=dims(1)
C     dimsd(2)=dims(2)
C     dimsd(3)=dims(3)
C     dimsd(4)=dims(4)
C     dimsd(5)=15
C
C     plhs(1) = mxCreateNumericArray(nd, dimsd, 
C    &          classid, 
C    &          complexflag)

C      plhs(2) = mxCreateNumericArray(nd, dimsd, 
C     &          classid, 
C     &          complexflag)


      wflux = mxGetPr(plhs(1))
C      itertrac = mxGetPr(plhs(2))
       write(line,*)'Implicit diffusion iterations:',int(iter) 
          k=mexPrintf(line//achar(13))
         
CC     Call the computational subroutines.
         call mit_implicit_rmix(dims(1),dims(2),dims(3),dims(4),
     &           %val(deltat), %val(iterret), %val(mask), %val(rac),
     &           %val(dzf), %val(dzc),%val(diffkz),%val(scalar),
     &           %val(wflux),errCode)
     
C     %val(itertrac),errCode)
     
CC     Check data type of input argument.
      if (errCode .ne. 0) then
CC         call mexErrMsgTxt('Input argument must be of type double.')
         write(line,*)'ERROR:SOLVE_TRIDIAGONAL:Single matrix.'
         k=mexPrintf(line//achar(13)) 
         write(line,*)'errCode=',errCode
         k=mexPrintf(line//achar(13)) 
         return     
      endif
      
      return
      end