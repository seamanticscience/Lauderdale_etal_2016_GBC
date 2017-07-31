#include "fintrf.h"
C======================================================================
C     MIT_CALC_DIVERGENCE calculates transport divergence and 
C                         integrates with depth.
C     [hflux,rflux]=mit_calc_divergence(dz,volc,kave,mask,
C                                       uflux,vflux,wflux);
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
      mwPointer dz, volc, mask
      mwPointer kave
      mwPointer uflux, vflux, wflux
      mwPointer hflux, rflux
      
C     integer*4, Dimension(:,:,:), Allocatable :: kave

C    For writing error to the screen
C     mexErrMsgIdAndTxt causes crash
      character*250 line
      integer*4 k
      integer*4 mexPrintf

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 7) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nInput     ',
C     &                           'Fourteen inputs required. ')
         write(line,*)'ERROR: Fourteen inputs required.'
         k=mexPrintf(line//achar(13)) 
         return
      elseif(nlhs .gt. 2) then
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
      ndim = mxGetNumberOfDimensions(prhs(7))
      call mxCopyPtrToInteger4(mxGetDimensions(prhs(7)),  
     &           dims, mxGetNumberOfDimensions(prhs(7)))
      
C      write(line,*)'ndims: ',ndim
C      k=mexPrintf(line//achar(13)) 
C      write(line,*)'dimensions: ',dims(1),'x',dims(2),
C     &                        'x',dims(3),'x',dims(4)
C      k=mexPrintf(line//achar(13)) 
    
C     Create Fortran arrays from the input arguments.
      dz    = mxGetPr(prhs(1))
      volc  = mxGetPr(prhs(2))
      kave  = mxGetPr(prhs(3))
      mask  = mxGetPr(prhs(4))
      uflux = mxGetPr(prhs(5))
      vflux = mxGetPr(prhs(6))
      wflux = mxGetPr(prhs(7))

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
      plhs(2) = mxCreateNumericArray(ndim-1, fdims, 
     &          classid, 
     &          complexflag)

      hflux = mxGetPr(plhs(1))
      rflux = mxGetPr(plhs(2))

CC     Call the computational subroutine.

         call calc_flux_divergence(dims(1),dims(2),dims(3),dims(4),
     &                %val(kave),%val(dz),%val(volc),%val(mask),
     &                %val(uflux),%val(vflux),%val(wflux),
     &                %val(hflux),%val(rflux))
     
      return
      end
      
CC------------------------COMPUTATIONAL SUBROUTINE------------------------
      
      subroutine calc_flux_divergence(iMax,jMax,Nr,Nt,kLev,dz,volc,
     &                       maskC,fZon,fMer,fVerT,hflux,rflux)
      
      implicit none
      
      integer*4 iMax,jMax,Nr,Nt,i,j,k,t
      integer*4 ip1,jp1,km1,kp1

      real*8 dz   (Nr)
C      real*8 recip_drF(Nr)
      real*8 volc (iMax,jMax,Nr)
      real*8 recip_volc(iMax,jMax,Nr)
      real*8 maskC(iMax,jMax,Nr)
      real*8 hdiv (iMax,jMax,Nr)
      real*8 rdiv (iMax,jMax,Nr)
      real*8 kLev (iMax,jMax,Nt)
      real*8 hflux(iMax,jMax,Nt)
      real*8 rflux(iMax,jMax,Nt)
      real*8 fZon (iMax,jMax,Nr,Nt)
      real*8 fMer (iMax,jMax,Nr,Nt)
      real*8 fVerT(iMax,jMax,Nr,Nt)

C From config_summary
C rkSign =   /* index orientation relative to vertical coordinate */
C -1.000000000000000E+00
      real*8 rkSign
      PARAMETER(rkSign=-1.D0)
      
C     Initialize fluxes to zero and calculate reciprocals
      do t=1,Nt
       do j=1,jMax
        do i=1,iMax  
         hflux(i,j,t)=0.D0
         rflux(i,j,t)=0.D0
        enddo
       enddo  
      enddo

C      do k=1,Nr
C       recip_drF(k)=1.D0/dz(k)
C      enddo
      
      do j=1,jMax
       do i=1,iMax
        do k=1,Nr 
         recip_volc(i,j,k)=1.D0/volc(i,j,k)
        enddo
       enddo
      enddo

C--   Divergence of fluxes
C     Anelastic: scale vertical fluxes by rhoFac and leave Horizontal fluxes unchanged
C     for Stevens OBC: keep only vertical diffusive contribution on boundaries
      do t=1,Nt
       do k=1,Nr
        do j=1,jMax
         do i=1,iMax
          if (i .EQ. iMax) THEN
           ip1=1
          else
           ip1=i+1
          endif
          
          jp1=MIN(j+1,jMax)
          
          km1  = MAX(1 ,k-1)
          kp1  = MIN(Nr,k+1)
        
          hdiv(i,j,k) = -maskC(i,j,k)*recip_volc(i,j,k)
     &            *((fZon(ip1, j ,k,t)-fZon(i,j,k,t))
     &             +(fMer( i ,jp1,k,t)-fMer(i,j,k,t))
     &             )
       
          rdiv(i,j,k) = -maskC(i,j,k)*recip_volc(i,j,k)
     &             *(fVerT(i,j,kp1,t)-fVerT(i,j,k,t))*rkSign
         enddo
        enddo
       enddo

C     Integrate vertically          
       do j=1,jMax
        do i=1,iMax  
         do k=1,INT(kLev(i,j,t))    
          hflux(i,j,t)=hflux(i,j,t)
     &      +(hdiv(i,j,k)*dz(k)*maskC(i,j,k))
     
          rflux(i,j,t)=rflux(i,j,t)
     &      +(rdiv(i,j,k)*dz(k)*maskC(i,j,k))
         enddo
        enddo
       enddo
      enddo
            
      return
      end