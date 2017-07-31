#include "fintrf.h"
C======================================================================
C     MIT_CALC_DIVERGENCE calculates transport divergence and 
C                         integrates with depth.
C     [rint]=mit_integrate_vert(dz,kave,mask,fld);
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

CC     Array information:
C      mwPointer mxGetM, mxGetN
C      mwPointer mrows, ncols
C      mwSize size
      
C     Arguments for mxCreateNumericArray
      integer*4 classid, complexflag
      mwSize ndim, dims(4), fdims(3)
    
C     Pointers to input/output mxArrays:
      mwPointer dz, mask
      mwPointer kave
      mwPointer fld
      mwPointer rint
      
C     integer*4, Dimension(:,:,:), Allocatable :: kave

C    For writing error to the screen
C     mexErrMsgIdAndTxt causes crash
      character*250 line
      integer*4 k
      integer*4 mexPrintf

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 4) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nInput     ',
C     &                           'Fourteen inputs required. ')
         write(line,*)'ERROR: Four inputs required.'
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
      ndim = mxGetNumberOfDimensions(prhs(4))
      call mxCopyPtrToInteger4(mxGetDimensions(prhs(4)),  
     &           dims, mxGetNumberOfDimensions(prhs(4)))
      
C      write(line,*)'ndims: ',ndim
C      k=mexPrintf(line//achar(13)) 
C      write(line,*)'dimensions: ',dims(1),'x',dims(2),
C     &                        'x',dims(3),'x',dims(4)
C      k=mexPrintf(line//achar(13)) 
    
C     Create Fortran arrays from the input arguments.
      dz    = mxGetPr(prhs(1))
      kave  = mxGetPr(prhs(2))
      mask  = mxGetPr(prhs(3))
      fld   = mxGetPr(prhs(4))

CC     Copy kave to an integer variable
C      mrows = mxGetM(prhs(3))
C      ncols = mxGetN(prhs(3))
C      size = mrows*ncols
C      call mxCopyPtrToInteger4(kave_ptr, kave, size)

C     Output will be uflux,vflux,wflux
C     Create matrix for the return argument.
      classid = mxClassIDFromClassName('double')
      complexflag = 0

C     Dimensions for integrated flux divergences are 3d (i,j,t)
      fdims(1)=dims(1)
      fdims(2)=dims(2)
      fdims(3)=dims(4)

      plhs(1) = mxCreateNumericArray(ndim-1, fdims, 
     &          classid, 
     &          complexflag)
      
      rint = mxGetPr(plhs(1))

CC     Call the computational subroutine.

         call calc_vertical_int(dims(1),dims(2),dims(3),dims(4),
     &                %val(kave),%val(dz),%val(mask),
     &                %val(fld),
     &                %val(rint))
     
      return
      end
      
CC------------------------COMPUTATIONAL SUBROUTINE------------------------
      
      subroutine calc_vertical_int(iMax,jMax,Nr,Nt,kLev,dz,
     &                       maskC,fld,rint)
      
      implicit none
      
      integer*4 iMax,jMax,Nr,Nt,i,j,k,t

      real*8 dz   (Nr)
      real*8 maskC(iMax,jMax,Nr)
      real*8 rint (iMax,jMax,Nr)
      real*8 kLev (iMax,jMax,Nt)
      real*8 fld  (iMax,jMax,Nr,Nt)

C     Initialize fluxes to zero
      do t=1,Nt
       do j=1,jMax
        do i=1,iMax  
         rint(i,j,t)=0.D0
        enddo
       enddo  
      enddo

C     Integrate vertically
      do t=1,Nt
       do j=1,jMax
        do i=1,iMax  
         do k=1,INT(kLev(i,j,t))    
          rint(i,j,t)=rint(i,j,t)
     &      +(fld(i,j,k,t)*dz(k)*maskC(i,j,k))
         enddo
        enddo
       enddo
      enddo
                  
      return
      end