#include "fintrf.h"
C-----------------------------------------------------------------------------------------
C     diffusion_flux calculates explicit fld diffusive transports
C     [uflux,vflux,wflux]=diffusion_flux( dx ,dy ,dz ,
C                                           umask,vmask,cmask.
C                                           rax,ray,rac,
C                                           kux,kuz,kvy,kvz,
C                                           kwx,kwy,fld);
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
      real*8, allocatable, dimension(:) :: dyin
      real*8, allocatable, dimension(:) :: dxin
      real*8, allocatable, dimension(:) :: racin
      real*8, allocatable, dimension(:) :: rayin
      real*8, allocatable, dimension(:) :: raxin
      real*8, allocatable, dimension(:) :: mskcin
      real*8, allocatable, dimension(:) :: mskvin
      real*8, allocatable, dimension(:) :: mskuin
      real*8, allocatable, dimension(:) :: kuxin
      real*8, allocatable, dimension(:) :: kuzin
      real*8, allocatable, dimension(:) :: kvyin
      real*8, allocatable, dimension(:) :: kvzin
      real*8, allocatable, dimension(:) :: kwxin
      real*8, allocatable, dimension(:) :: kwyin
      real*8, allocatable, dimension(:) :: fin

C     Allocatable arrays for fortran data inputs to computational routine
      real*8, allocatable, dimension(:,:)     :: dy
      real*8, allocatable, dimension(:,:)     :: dx
      real*8, allocatable, dimension(:,:)     :: rac
      real*8, allocatable, dimension(:,:,:)   :: ray
      real*8, allocatable, dimension(:,:,:)   :: rax
      real*8, allocatable, dimension(:,:,:)   :: mskc
      real*8, allocatable, dimension(:,:,:)   :: mskv
      real*8, allocatable, dimension(:,:,:)   :: msku
      real*8, allocatable, dimension(:,:,:,:) :: kux
      real*8, allocatable, dimension(:,:,:,:) :: kuz
      real*8, allocatable, dimension(:,:,:,:) :: kvy
      real*8, allocatable, dimension(:,:,:,:) :: kvz
      real*8, allocatable, dimension(:,:,:,:) :: kwx
      real*8, allocatable, dimension(:,:,:,:) :: kwy
      real*8, allocatable, dimension(:,:,:,:) :: fld

C     Allocatable arrays for fortran data outputs  
      real*8, allocatable, dimension(:,:,:,:) :: uflux
      real*8, allocatable, dimension(:,:,:,:) :: vflux
      real*8, allocatable, dimension(:,:,:,:) :: wflux

C---------------VALIDATE------------INPUTS--------------------------------
C     Check for proper number of arguments. 
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Calculating Explicit Diffusive Fluxes')
      call mexWarnMsgTxt('Check for proper number of arguments.')

      write(line,*)'Number of inputs is ',nrhs
      k=mexPrintf(line//achar(13)) 
      write(line,*)'Number of outputs is ',nlhs
      k=mexPrintf(line//achar(13)) 
#endif
      if(nrhs .ne. 16) then
         call mexErrMsgTxt('ERROR: 16 inputs required')
      elseif(nlhs .gt. 3) then
         call mexErrMsgTxt('ERROR: Too many output arguments.')
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
         dims(4)=1
         ndim=ndim+1 
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
          
C Finally allocate the input array dimensions  
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate the input array dimensions') 
      call mexWarnMsgTxt('allocate dxin array dimensions'   )       
#endif
      allocate (dxin  (mxGetNumberOfElements(prhs(1))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate dyin array dimensions'   )       
#endif
      allocate (dyin  (mxGetNumberOfElements(prhs(2))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate dz array dimensions'   )       
#endif
      allocate (dz   (mxGetNumberOfElements(prhs(3))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mskuin array dimensions' )  
#endif
      allocate (mskuin(mxGetNumberOfElements(prhs(4))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mskvin array dimensions' )    
#endif
      allocate (mskvin(mxGetNumberOfElements(prhs(5))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mskcin array dimensions' )    
#endif
      allocate (mskcin(mxGetNumberOfElements(prhs(6))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate raxin array dimensions')  
#endif
      allocate (raxin(mxGetNumberOfElements(prhs(7))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate rayin array dimensions'  )
#endif
      allocate (rayin(mxGetNumberOfElements(prhs(8))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate racin array dimensions')   
#endif
      allocate (racin(mxGetNumberOfElements(prhs(9))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kuxin array dimensions')   
#endif
      allocate (kuxin(mxGetNumberOfElements(prhs(10))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kuzin array dimensions' )  
#endif
      allocate (kuzin(mxGetNumberOfElements(prhs(11))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kvyin array dimensions' )    
#endif
      allocate (kvyin(mxGetNumberOfElements(prhs(12))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kvzin array dimensions')   
#endif
      allocate (kvzin(mxGetNumberOfElements(prhs(13))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kwxin array dimensions' )  
#endif
      allocate (kwxin(mxGetNumberOfElements(prhs(14))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kwyin array dimensions' )    
#endif
      allocate (kwyin(mxGetNumberOfElements(prhs(15))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate fin array dimensions' )    
#endif
      allocate (fin(mxGetNumberOfElements(prhs(16))))
      
C Allocate the input to the computational function
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate dx array dimensions')    
#endif
      allocate (dx(dims(1),dims(2)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate dy array dimensions')    
#endif
      allocate (dy(dims(1),dims(2))) 
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate msku array dimensions')    
#endif
      allocate (msku(dims(1),dims(2),dims(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mskv array dimensions')    
#endif
      allocate (mskv(dims(1),dims(2),dims(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mskc array dimensions')    
#endif
      allocate (mskc(dims(1),dims(2),dims(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate rax array dimensions')    
#endif
      allocate (rax(dims(1),dims(2),dims(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate ray array dimensions')    
#endif
      allocate (ray(dims(1),dims(2),dims(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate rac array dimensions')    
#endif
      allocate (rac(dims(1),dims(2)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kux array dimensions')    
#endif
      allocate (kux(dims(1),dims(2),dims(3),dims(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kuz array dimensions')    
#endif
      allocate (kuz(dims(1),dims(2),dims(3),dims(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kvy array dimensions')    
#endif
      allocate (kvy(dims(1),dims(2),dims(3),dims(4)))    
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kvz array dimensions')    
#endif
      allocate (kvz(dims(1),dims(2),dims(3),dims(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kwx array dimensions')    
#endif
      allocate (kwx(dims(1),dims(2),dims(3),dims(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kwy array dimensions')    
#endif
      allocate (kwy(dims(1),dims(2),dims(3),dims(4)))          
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate fld array dimensions')    
#endif
      allocate (fld(dims(1),dims(2),dims(3),dims(4)))     
      
C Allocate the output from the computational function      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate hflux array dimensions')    
#endif
      allocate (uflux(dims(1),dims(2),dims(3),dims(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate rflux array dimensions')    
#endif
      allocate (vflux(dims(1),dims(2),dims(3),dims(4)))          
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate rflux array dimensions')    
#endif
      allocate (wflux(dims(1),dims(2),dims(3),dims(4)))        

      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill arrays from the input.') 
      call mexWarnMsgTxt('Fill dxin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(1)), dxin, 
     & mxGetNumberOfElements(prhs(1)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill dyin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), dyin, 
     & mxGetNumberOfElements(prhs(2)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill dz array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), dz, 
     & mxGetNumberOfElements(prhs(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill mskuin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(4)), mskuin, 
     & mxGetNumberOfElements(prhs(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill mskvin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(5)), mskvin, 
     & mxGetNumberOfElements(prhs(5)))          
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill mskcin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(6)), mskcin, 
     & mxGetNumberOfElements(prhs(6)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill raxin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(7)), raxin, 
     & mxGetNumberOfElements(prhs(7)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill rayin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(8)), rayin, 
     & mxGetNumberOfElements(prhs(8)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill racin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(9)), racin, 
     & mxGetNumberOfElements(prhs(9)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill kuxin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(10)), kuxin, 
     & mxGetNumberOfElements(prhs(10)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill kuzin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(11)), kuzin, 
     & mxGetNumberOfElements(prhs(11)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill kvyin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(12)), kvyin, 
     & mxGetNumberOfElements(prhs(12)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill kvzin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(13)), kvzin, 
     & mxGetNumberOfElements(prhs(13)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill kwxin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(14)), kwyin, 
     & mxGetNumberOfElements(prhs(14)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill kwyin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(15)), kwyin, 
     & mxGetNumberOfElements(prhs(15)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill fin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(16)), fin, 
     & mxGetNumberOfElements(prhs(16)))          
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Create arrays for the output arguments.')  
#endif
      classid = mxClassIDFromClassName('double')
      complexflag = 0  
      do iter=1,nlhs
         plhs(iter) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)
      enddo

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check arguments for NANs')    
      call mexWarnMsgTxt('Check dxin')
#endif
      if (mxIsNaN(dxin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input dxin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(1))
        if( mxIsNaN(dxin(iter)) .eq. 1 ) then
         dxin(iter)=0 .d 0
        endif
       enddo
      endif
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check dyin')
#endif
      if (mxIsNaN(dyin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input dyin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(2))
        if( mxIsNaN(dyin(iter)) .eq. 1 ) then
         dyin(iter)=0 .d 0
        endif
       enddo
      endif
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check dz')    
#endif
      if (mxIsNaN(dz) .eq. 1) then
C      dz is pretty important, so error if NaN
       call mexErrMsgTxt('ERROR: Input dz contains a NaN.')
      endif

#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Check mskuin')
#endif    
      if (mxIsNaN(mskuin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input mskuin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(4))
        if( mxIsNaN(mskuin(iter)) .eq. 1 ) then
         mskuin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check mskvin')
#endif
      if (mxIsNaN(mskvin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input mskvin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(5))
        if( mxIsNaN(mskvin(iter)) .eq. 1 ) then
         mskvin(iter)=0 .d 0
        endif
       enddo
      endif
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check mskcin')    
#endif
      if (mxIsNaN(mskcin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input mskcin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(6)) 
        if( mxIsNaN(mskcin(iter)) .eq. 1 ) then
         mskcin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check raxin')    
#endif
      if (mxIsNaN(raxin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input raxin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(7))
        if( mxIsNaN(raxin(iter)) .eq. 1 ) then
         raxin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check rayin')    
#endif
       if (mxIsNaN(rayin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input rayin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(8))
        if( mxIsNaN(rayin(iter)) .eq. 1 ) then
         rayin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check racin')    
#endif
      if (mxIsNaN(racin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input racin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(9))
        if( mxIsNaN(racin(iter)) .eq. 1 ) then
         racin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check kuxin')    
#endif
      if (mxIsNaN(kuxin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input kuxin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(10))
        if( mxIsNaN(kuxin(iter)) .eq. 1 ) then
         kuxin(iter)=0 .d 0
        endif
       enddo
      endif    

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check kuzin')    
#endif
      if (mxIsNaN(kuzin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input kuzin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(11))
        if( mxIsNaN(kuzin(iter)) .eq. 1 ) then
         kuzin(iter)=0 .d 0
        endif
       enddo
      endif    

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check kvyin')    
#endif
      if (mxIsNaN(kvyin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input kvyin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(12))
        if( mxIsNaN(kvyin(iter)) .eq. 1 ) then
         kvyin(iter)=0 .d 0
        endif
       enddo
      endif    

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check kvzin')    
#endif
      if (mxIsNaN(kvzin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input kvzin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(13))
        if( mxIsNaN(kvzin(iter)) .eq. 1 ) then
         kvzin(iter)=0 .d 0
        endif
       enddo
      endif    

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check kwxin')    
#endif
      if (mxIsNaN(kwxin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input kwxin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(14))
        if( mxIsNaN(kwxin(iter)) .eq. 1 ) then
         kwxin(iter)=0 .d 0
        endif
       enddo
      endif    

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check kwyin')    
#endif
      if (mxIsNaN(kwyin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input kwyin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(15))
        if( mxIsNaN(kwyin(iter)) .eq. 1 ) then
         kwyin(iter)=0 .d 0
        endif
       enddo
      endif    

#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Check fin')
#endif    
      if (mxIsNaN(fin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input fin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(16))
        if( mxIsNaN(fin(iter)) .eq. 1 ) then
         fin(iter)=0 .d 0
        endif
       enddo
      endif    
      
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Reshape input arrays') 
#endif   
      dx   = reshape(dxin  ,(/dims(1),dims(2)/))
      dy   = reshape(dyin  ,(/dims(1),dims(2)/))
      msku = reshape(mskuin,(/dims(1),dims(2),dims(3)/))
      mskv = reshape(mskvin,(/dims(1),dims(2),dims(3)/))
      mskc = reshape(mskcin,(/dims(1),dims(2),dims(3)/))
      rax  = reshape(raxin ,(/dims(1),dims(2),dims(3)/))
      ray  = reshape(rayin ,(/dims(1),dims(2),dims(3)/))
      rac  = reshape(racin ,(/dims(1),dims(2)/))
      kux  = reshape(kuxin ,(/dims(1),dims(2),dims(3),dims(4)/))
      kuz  = reshape(kuzin ,(/dims(1),dims(2),dims(3),dims(4)/))
      kvy  = reshape(kvyin ,(/dims(1),dims(2),dims(3),dims(4)/))
      kvz  = reshape(kvzin ,(/dims(1),dims(2),dims(3),dims(4)/))
      kwx  = reshape(kwxin ,(/dims(1),dims(2),dims(3),dims(4)/))
      kwy  = reshape(kwyin ,(/dims(1),dims(2),dims(3),dims(4)/))
      fld  = reshape(fin   ,(/dims(1),dims(2),dims(3),dims(4)/))

C Fill the output arrays with zeros
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Fill output arrays with zeros') 
#endif   
      do t=1,dims(4)
       do k=1,dims(3)
        do j=1,dims(2)
         do i=1,dims(1)
          uflux(i,j,k,t)=0.D0
          vflux(i,j,k,t)=0.D0
          wflux(i,j,k,t)=0.D0
         enddo
        enddo  
       enddo
      enddo
      
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Call computational subroutine')  
#endif  
      
CC----------------DIFFUSION IN X DIRECTION----------------
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt(
     & 'Calculating explicit diffusion in the X direction') 
#endif  
         call diffx(dims(1),dims(2),dims(3),dims(4),
     &           mskc,dz,dx,rax,kux,kuz,fld,uflux)

CC----------------DIFFUSION IN Y DIRECTION----------------
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt(
     & 'Calculating explicit diffusion in the Y direction') 
#endif  
         call diffy(dims(1),dims(2),dims(3),dims(4),
     &           mskc,dz,dy,ray,kvy,kvz,fld,vflux)

CC----------------DIFFUSION IN Z DIRECTION----------------
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt(
     & 'Calculating explicit diffusion in the Z direction') 
#endif  
         call diffz(dims(1),dims(2),dims(3),dims(4),
     &           msku,mskv,mskc,dx,dy,rac,kwx,
     &           kwy,fld,wflux)
     
C     Load the data into the pointers, which are the output to MATLAB.
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Load the data into the output pointers')    
#endif
      size=dims(1)*dims(2)*dims(3)*dims(4)
      call mxCopyReal8ToPtr(uflux, mxGetPr(plhs(1)),size)
      call mxCopyReal8ToPtr(vflux, mxGetPr(plhs(2)),size)
      call mxCopyReal8ToPtr(wflux, mxGetPr(plhs(3)),size)
     
      return
      end

CC------------------------COMPUTATIONAL SUBROUTINES------------------------

      subroutine diffx(imax,jmax,nr,nt,
     &     maskc,dz,dx,xa,kux,kuz,fld,flux )

      implicit none

      mwIndex imax, jmax, nr, nt
      mwIndex i, j, k, t, kp1, km1, im1
      real*8 dz      (nr)
      real*8 recip_dz(nr)
      real*8 dx      (imax,jmax)
      real*8 recip_dx(imax,jmax)
      real*8 maskc   (imax,jmax,nr)
      real*8 xa      (imax,jmax,nr)
      real*8 kux     (imax,jmax,nr,nt)
      real*8 kuz     (imax,jmax,nr,nt)
      real*8 fld  (imax,jmax,nr,nt)
      real*8 flux      (imax,jmax,nr,nt)
      real*8 dtdz
      real*8 half
      parameter(half=1.d0/2.d0)

C     initialise flux to zero
      do t=1,nt
       do k=1,nr
         recip_dz(k)=1.d0/dz(k)  
        do j=1,jmax
         do i=1,imax
          flux(i,j,k,t) = 0.d0
          
          if (maskc(i,j,1) .ne. 0.d0) then
           recip_dx(i,j)=1.d0/dx(i,j)
          else
           recip_dx(i,j)=0.d0
          endif

         enddo
        enddo
       enddo
      enddo
               
C     area integrated zonal flux
      do t=1,nt
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          if (i .eq. 1) then
           im1=imax
          else
           im1=i-1
          endif
          flux(i,j,k,t) = flux(i,j,k,t)
     &      -xa(i,j,k)
     &      *kux(i,j,k,t)
     &      *recip_dx(i,j)
     &      *(fld(i,j,k,t)-fld(im1,j,k,t))
         enddo
        enddo
       enddo
      enddo
      
C     off-diagonal components of horizontal flux
      do t=1,nt
       do k=1,nr
       km1 = max(k-1,1)
       kp1 = min(k+1,nr)
C    vertical gradients interpolated to u points
        do j=1,jmax
         do i=1,imax
           if (i .eq. 1) then
            im1=imax
           else
            im1=i-1
           endif
             
           dtdz =  half*(
     &      half*recip_dz(k)*
     &        (maskc(im1,j,k)*(fld(im1,j,km1,t)-fld(im1,j,k,t))
     &        +maskc( i ,j,k)*(fld( i ,j,km1,t)-fld( i ,j,k,t))
     &        )
     &      +half*recip_dz(kp1)*
     &        (maskc(im1,j,kp1)*(fld(im1,j,k,t)-fld(im1,j,kp1,t))
     &        +maskc( i ,j,kp1)*(fld( i ,j,k,t)-fld( i ,j,kp1,t))
     &        )
     &      )
  
           flux(i,j,k,t) = flux(i,j,k,t) - xa(i,j,k)*kuz(i,j,k,t)*dtdz
         enddo
        enddo
       enddo
      enddo
       
      return
      end

      subroutine diffy(imax, jmax, nr, nt,
     &     maskc, dz, dy, ya,kvy, kvz,fld,flux )

      implicit none

      mwIndex imax, jmax, nr, nt
      mwIndex i, j, k, t, kp1, km1, jm1
      real*8 dz      (nr)
      real*8 recip_dz(nr)
      real*8 dy      (imax,jmax)
      real*8 recip_dy(imax,jmax)
      real*8 ya      (imax,jmax,nr)
      real*8 maskc   (imax,jmax,nr)
      real*8 kvy     (imax,jmax,nr,nt)
      real*8 kvz     (imax,jmax,nr,nt)
      real*8 fld     (imax,jmax,nr,nt)
      real*8 flux    (imax,jmax,nr,nt)
      real*8 dtdz
      real*8 half
      parameter(half=1.d0/2.d0)

C     initialise flux to zero
      do t=1,nt
       do k=1,nr
         recip_dz(k)=1.d0/dz(k)  
        do j=1,jmax
         do i=1,imax
          flux(i,j,k,t) = 0.d0
       
          if (maskc(i,j,1) .ne. 0.d0) then
           recip_dy(i,j)=1.d0/dy(i,j)
          else
           recip_dy(i,j)=0.d0
          endif
         enddo
        enddo
       enddo
      enddo
      
C     area integrated zonal flux
      do t=1,nt
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          jm1=max(j-1, 1 )
          flux(i,j,k,t) = flux(i,j,k,t)
     &      -ya(i,j,k)
     &      *kvy(i,j,k,t)
     &      *recip_dy(i,j)
     &      *(fld(i,j,k,t)-fld(i,jm1,k,t))
         enddo
        enddo
       enddo
      enddo
      
      do t=1,nt
       do k=1,nr
       km1 = max(k-1,1)
       kp1 = min(k+1,nr)
C    vertical gradients interpolated to u points
        do j=1,jmax
         do i=1,imax

          jm1=max(j-1, 1 )

C    vertical gradients interpolated to v points
           dtdz =  half*(
     &     half*recip_dz(k)*
     &        (maskc(i,jm1,k)*(fld(i,jm1,km1,t)-fld(i,jm1,k,t))
     &        +maskc(i, j ,k)*(fld(i, j ,km1,t)-fld(i, j ,k,t))
     &        )
     &     +half*recip_dz(kp1)*
     &        (maskc(i,jm1,kp1)*(fld(i,jm1,k,t)-fld(i,jm1,kp1,t))
     &        +maskc(i, j ,kp1)*(fld(i, j ,k,t)-fld(i, j ,kp1,t))
     &        )
     &     )

           flux(i,j,k,t) = flux(i,j,k,t) - ya(i,j,k)*kvz(i,j,k,t)*dtdz
         enddo
        enddo
       enddo
      enddo
      
      return
      end

      subroutine diffz(imax, jmax, nr, nt,
     &     masku, maskv, maskc,dx, dy, ra, kwx, kwy,
     &     fld, flux )

      implicit none
      mwIndex mxIsNaN
      mwIndex imax, jmax, nr, nt, i, j, k, t
      mwIndex kp1, km1, jm1, im1, ip1, jp1
      real*8 dx      (imax,jmax)
      real*8 dy      (imax,jmax)
      real*8 recip_dy(imax,jmax)
      real*8 recip_dx(imax,jmax)
      real*8 ra      (imax,jmax)
      real*8 masku   (imax,jmax,nr)
      real*8 maskv   (imax,jmax,nr)
      real*8 maskc   (imax,jmax,nr)
      real*8 kwx     (imax,jmax,nr,nt)
      real*8 kwy     (imax,jmax,nr,nt)
      real*8 fld     (imax,jmax,nr,nt)
      real*8 flux    (imax,jmax,nr,nt)

      real*8 tmp1, tmp2, dtdx, dtdy
      real*8 half
      parameter(half=1.d0/2.d0)

#ifdef DEBUG_MESSAGES 
      mwIndex, external :: mexPrintf
      mwIndex msg
      character(len=250) :: line
#endif

C     initialise flux to zero
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('DIFFZ: set flux to zero') 
#endif  
      do t=1,nt
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          flux(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo

#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('DIFFZ: calculate 1/dx and 1/dy') 
#endif  
      
      do t=1,nt
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          if (maskc(i,j,1) .ne. 0.d0) then
           recip_dx(i,j)=(1.d0/dx(i,j))
           recip_dy(i,j)=(1.d0/dy(i,j))
          else
           recip_dx(i,j)=0.d0
           recip_dy(i,j)=0.d0
          endif
          if (mxIsNaN(recip_dx(i,j)) .eq. 1) then
C          catch a NAN?
           call mexErrMsgTxt('ERROR: NAN in recip_dx during DIFFZ.')
          endif
          if (mxIsNaN(recip_dy(i,j)) .eq. 1) then
C          catch a NAN?
           call mexErrMsgTxt('ERROR: NAN in recip_dy during DIFFZ.')
          endif
         enddo
        enddo
       enddo
      enddo

#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('DIFFZ: calculate gradients and flux') 
#endif  
C    horizontal gradients interpolated to w points
      do t=1,nt
       do k=2,nr
        do j=1,jmax
         do i=1,imax
       
         if (i .eq. 1) then
          im1=imax
         else
          im1=i-1
         endif
          
         if (i .eq. imax) then
          ip1=1
         else
          ip1=i+1
         endif

         km1 = max(k-1,1)
         kp1 = min(k+1,nr)
              
         if ( fld( i ,j,k,t) .ne. 0 .AND. 
     &        fld(im1,j,k,t) .ne. 0) then
            tmp1=half*(masku(ip1,j,k)
     &            *recip_dx(ip1,j)*(fld(ip1,j,k,t)-fld( i ,j,k,t))
     &           +masku(i,j,k)
     &            *recip_dx( i ,j)*(fld( i ,j,k,t)-fld(im1,j,k,t))
     &            )
             
            tmp2=half*(masku(ip1,j,km1)
     &            *recip_dx(ip1,j)*(fld(ip1,j,km1,t)-fld( i ,j,km1,t))
     &           +masku(i,j,km1)
     &            *recip_dx( i ,j)*(fld( i ,j,km1,t)-fld(im1,j,km1,t)) 
     &            )     
   
            dtdx = half*(tmp1+tmp2)
         else
            dtdx = 0.d0
         endif
                
         jm1=max(j-1, 1 )
         jp1=min(j+1,jmax)

         if ( fld(i, j ,k,t) .ne. 0 .AND. 
     &        fld(i,jp1,k,t) .ne. 0) then

            tmp1=half*(maskv(i,j,k)
     &            *recip_dy(i,j)*(fld(i,j,k,t)-fld(i,jm1,k,t))
     &           +maskv(i,j+1,k)
     &            *recip_dy(i,jp1)*(fld(i,jp1,k,t)-fld(i,j,k,t))
     &            )

            tmp2=half*(maskv(i,j,km1)
     &            *recip_dy(i,j)*(fld(i,j,km1,t)-fld(i,jm1,km1,t))
     &           +maskv(i,jp1,km1)
     &            *recip_dy(i,jp1)*(fld(i,jp1,km1,t)-fld(i,j,km1,t))
     &            )
          
            dtdy = half*(tmp1+tmp2)
         else
            dtdx = 0.d0
         endif
         
C    off-diagonal components of vertical flux
         flux(i,j,k,t) = flux(i,j,k,t)
     &       -ra(i,j)*maskc(i,j,k)
     &       *( kwx(i,j,k,t)*dtdx
     &        + kwy(i,j,k,t)*dtdy )
         enddo
        enddo
       enddo
      enddo 

      return
      end
