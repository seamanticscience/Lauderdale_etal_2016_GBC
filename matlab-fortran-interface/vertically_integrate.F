#include "fintrf.h"
C-----------------------------------------------------------------------------------------
C     vertically_integrate integrates with depth.
C     [flint]=vertically_integrate(dz,kave,mask,fld);
C     
C-----------------------------------------------------------------------------------------
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer*4 nlhs, nrhs
      
C     Function declarations:
      mwPointer mxCreateNumericArray
      mwPointer mxGetDimensions
      mwPointer mxGetPr
      mwSize mxGetNumberOfDimensions
      mwSize mxGetNumberOfElements
      mwIndex mxIsDouble, mxIsNumeric, mxIsNaN
      mwIndex mxClassIDFromClassName
      mwIndex, external :: mexPrintf
      mwIndex classid, complexflag
      mwIndex iter, k, i, j, t
      mwIndex dims(4)
      mwIndex fdims(3)
      mwSize ndim, nelm
      mwSize size
      character(len=250) :: line
          
C     Allocatable arrays for fortran data inputs to gateway routine
      real*8, allocatable, dimension(:) :: dz
      real*8, allocatable, dimension(:) :: mskin
      real*8, allocatable, dimension(:) :: kin
      real*8, allocatable, dimension(:) :: fin

C     Allocatable arrays for fortran data inputs to computational routine
      real*8, allocatable, dimension(:,:,:)   :: mask
      real*8, allocatable, dimension(:,:,:)   :: kave
      real*8, allocatable, dimension(:,:,:,:) :: fld
      
C     Allocatable arrays for fortran data outputs  
      real*8, allocatable, dimension(:,:,:)   :: flint

C---------------VALIDATE------------INPUTS--------------------------------
C     Check for proper number of arguments. 
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Calculating Vertical Integral')
      call mexWarnMsgTxt('Check for proper number of arguments.')
#endif
      if(nrhs .ne. 4) then
         call mexErrMsgTxt(
     &'ERROR: 4 inputs required (dz,kave,mask,fld)')
      elseif(nlhs .gt. 1) then
         call mexErrMsgTxt(
     &'ERROR: Too many output arguments (fld_int).')
      endif

C     Check that the inputs are a number.
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check that the input is a number.')
#endif
      do iter=1,nrhs
         if(mxIsNumeric(prhs(iter)) .eq. 0) then
            call mexErrMsgTxt(
     &           'ERROR: Input must be numeric.')
         endif
      end do
      
C     Check data type of input arguments.
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check data type of input argument.')
#endif
      do iter=1,nrhs
         if (mxIsDouble(prhs(iter)) .eq. 0) then
            call mexErrMsgTxt(
     &           'ERROR: Input must be double precision.')
         endif
      end do
      
C     Get the size of the input arrays using the last as template.     
      ndim = mxGetNumberOfDimensions(prhs(nrhs))
      nelm = mxGetNumberOfElements(prhs(nrhs))      

C Initialize dims with zeros      
      do iter=1,ndim
         dims(iter)=0
      enddo

C use fortran "sizeof" to determine if dims is int4 or int8 
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt(
     & 'Use FORTRAN "sizeof" to determine if dims is int32 or int64 ')
      write(line,*)'size of dims byte storage= ',sizeof(dims(1))
      k=mexPrintf(line//achar(13)) 
#endif           
      
C     Make sure the right mxCopy routine is called!
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Get dimensions of last input field')
#endif    
      if (sizeof(dims(1)) .eq. 1) then
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Calling mxCopyPtrToInteger1')
#endif  
         call mxCopyPtrToInteger1(mxGetDimensions(prhs(nrhs)),  
     &             dims, mxGetNumberOfDimensions(prhs(nrhs)))
      elseif (sizeof(dims(1)) .eq. 2) then
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Calling mxCopyPtrToInteger2')
#endif  
         call mxCopyPtrToInteger2(mxGetDimensions(prhs(nrhs)),  
     &             dims, mxGetNumberOfDimensions(prhs(nrhs)))
      elseif (sizeof(dims(1)) .eq. 4) then
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Calling mxCopyPtrToInteger4')
#endif  
         call mxCopyPtrToInteger4(mxGetDimensions(prhs(nrhs)),  
     &             dims, mxGetNumberOfDimensions(prhs(nrhs)))
      elseif (sizeof(dims(1)) .eq. 8) then
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Calling mxCopyPtrToInteger8')
#endif  
         call mxCopyPtrToInteger8(mxGetDimensions(prhs(nrhs)),  
     &             dims, mxGetNumberOfDimensions(prhs(nrhs)))
      else
C cannot determine what size of integer is being used, so error
       call mexErrMsgTxt(
     &  'cannot determine what size of integer is being used for dims')            
      endif
    
C want to make sure if a 3d variable is given, can still handle everything
C i.e. the last dimension is not zero or a wild number.
      if (ndim.lt.4) then
         dims(4)=1.d0 
         ndim   =ndim+1
      endif

#ifdef DEBUG_MESSAGES
      write(line,*)'ndims: ',ndim
      k=mexPrintf(line//achar(13)) 
      write(line,*)'dimensions:\n ',dims(1),'\nx',dims(2),
     &                        '\nx',dims(3),'\nx',dims(4)
      k=mexPrintf(line//achar(13)) 
      write(line,*)'nelements: ',nelm
      k=mexPrintf(line//achar(13))   
#endif

C     Dimensions for integrated flux divergences are 3d (i,j,t)
      fdims(1)=dims(1)
      fdims(2)=dims(2)
      fdims(3)=dims(4)

#ifdef DEBUG_MESSAGES
      write(line,*)'dimensions:\n ',fdims(1),'\nx',fdims(2),
     &              '\nx',fdims(3) 
      k=mexPrintf(line//achar(13)) 
#endif    

C Finally allocate the input array dimensions  
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate the input array dimensions') 
      call mexWarnMsgTxt('allocate dz array dimensions'   )       
#endif
      allocate (dz   (mxGetNumberOfElements(prhs(1))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kin array dimensions'  )
#endif
      allocate (kin  (mxGetNumberOfElements(prhs(2))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mskin array dimensions')   
#endif
      allocate (mskin(mxGetNumberOfElements(prhs(3))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate fldin array dimensions' )  
#endif
      allocate (fin  (mxGetNumberOfElements(prhs(4))))

C Allocate the input to the computational function      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kave array dimensions')    
#endif
      allocate (kave(dims(1),dims(2),dims(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mask array dimensions')    
#endif
      allocate (mask(dims(1),dims(2),dims(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate fld array dimensions')    
#endif
      allocate (fld(dims(1),dims(2),dims(3),dims(4)))

C Allocate the output from the computational function      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate flint array dimensions')    
#endif
      allocate (flint(fdims(1),fdims(2),fdims(3)))

C     Fill Fortran arrays from the input argument pointers.
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill arrays from the input.') 
      call mexWarnMsgTxt('Fill dz array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(1)), dz, 
     & mxGetNumberOfElements(prhs(1)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill kin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), kin, 
     & mxGetNumberOfElements(prhs(2)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill mskin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), mskin, 
     & mxGetNumberOfElements(prhs(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill fin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(4)), fin, 
     & mxGetNumberOfElements(prhs(4)))

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Create arrays for the output argument.')  
#endif
      classid = mxClassIDFromClassName('double')
      complexflag = 0  
      plhs(1) = mxCreateNumericArray(ndim-1, fdims, 
     &          classid, 
     &          complexflag)

C Check arguments for NANs in the input
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check arguments for NANs')    
      call mexWarnMsgTxt('Check dz')    
#endif
      if (mxIsNaN(dz) .eq. 1) then
C      dz is pretty important, so error if NaN
       call mexErrMsgTxt('ERROR: Input dz contains a NaN.')
      endif
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check kin')
#endif      
      if (mxIsNaN(kin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input kin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(2))
        if( mxIsNaN(kin(iter)) .eq. 1 ) then
         kin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check mskin')
#endif          
      if (mxIsNaN(mskin) .eq. 1) then
#ifdef DEBUG_MESSAGES    
         call mexWarnMsgTxt(
     &      'Input mskin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(3))
        if( mxIsNaN(mskin(iter)) .eq. 1 ) then
         mskin(iter)=0 .d 0
        endif
       enddo
      endif
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check fin')
#endif      
      if (mxIsNaN(fin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input ufin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(4))
        if( mxIsNaN(fin(iter)) .eq. 1 ) then
         fin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Reshape input arrays') 
#endif   
      kave = reshape(kin  ,(/dims(1),dims(2),dims(4)/))
      mask = reshape(mskin,(/dims(1),dims(2),dims(3)/))
      fld  = reshape(fin  ,(/dims(1),dims(2),dims(3),dims(4)/))

C Fill the output arrays with zeros
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Fill output array with zeros') 
#endif 
      do t=1,fdims(3)
       do j=1,fdims(2)
        do i=1,fdims(1)
         flint(i,j,t)=0.D0
        enddo
       enddo  
      enddo

#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Call computational subroutine')  
#endif  
         call calc_vertical_int(dims(1),dims(2),dims(3),dims(4),
     &                kave,dz,mask,fld,flint)
     
C     Load the data into the pointers, which are the output to MATLAB.
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Load the data into the output pointer')    
#endif
      size=fdims(1)*fdims(2)*fdims(3)
      call mxCopyReal8ToPtr(flint, mxGetPr(plhs(1)), size)
     
      return
      end
      
CC------------------------COMPUTATIONAL SUBROUTINE------------------------
      
      subroutine calc_vertical_int(imax,jmax,nr,nt,klev,dz,
     &                       maskc,fldin,fint)
      
      implicit none
      
      mwIndex imax,jmax,nr,nt,i,j,k,t
      real*8 dz   (nr)
      real*8 maskc(imax,jmax,nr)
      real*8 fint(imax,jmax,nr)
      real*8 klev (imax,jmax,nt)
      real*8 fldin(imax,jmax,nr,nt)

#ifdef DEBUG_MESSAGES 
      mwIndex, external :: mexPrintf
      mwIndex msg
      character(len=250) :: line
#endif

C     Initialize fluxes to zero
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('zero the outputs')    
      write(line,*)'dimensions:\ni:',imax,'\nj:',jmax,
     &                        '\nk',nr,'\nt:',nt 
      msg=mexPrintf(line//achar(13)) 
#endif
      do t=1,nt
       do j=1,jmax
        do i=1,imax  
         fint(i,j,t)=0.D0
        enddo
       enddo  
      enddo

C     Integrate vertically
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('integrate vertically')  

#endif  
      do t=1,nt
       do j=1,jmax
        do i=1,imax  
         do k=1,INT(klev(i,j,t))    
          fint(i,j,t)=fint(i,j,t)
     &      +(fldin(i,j,k,t)*dz(k)*maskc(i,j,k))
         enddo
        enddo
       enddo
      enddo
                  
      return
      end
