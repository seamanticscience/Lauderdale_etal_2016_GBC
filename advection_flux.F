#include "fintrf.h"
C-----------------------------------------------------------------------------------------
C     advection_flux calls advection schemes for fld transports
C     [uflux,vflux,wflux]=advection_flux(timestep,advscheme,
C                                  dx,dy,dz,
C                                  rax,ray,rac,
C                                  umask,vmask,cmask,
C                                  uvel,vvel,wvel,
C                                  fld);
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

C     Arguments for computational routine:
      real*8 iAdvSch, deltat

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
      real*8, allocatable, dimension(:) :: uin
      real*8, allocatable, dimension(:) :: vin
      real*8, allocatable, dimension(:) :: win
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
      real*8, allocatable, dimension(:,:,:,:) :: uvel
      real*8, allocatable, dimension(:,:,:,:) :: vvel
      real*8, allocatable, dimension(:,:,:,:) :: wvel
      real*8, allocatable, dimension(:,:,:,:) :: fld

C     Allocatable arrays for fortran data outputs  
      real*8, allocatable, dimension(:,:,:,:) :: uflux
      real*8, allocatable, dimension(:,:,:,:) :: vflux
      real*8, allocatable, dimension(:,:,:,:) :: wflux
      
C---------------VALIDATE------------INPUTS--------------------------------
C     Check for proper number of arguments. 
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Calculating Advective Fluxes')
      call mexWarnMsgTxt('Check for proper number of arguments.')

      write(line,*)'Number of inputs is ',nrhs
      k=mexPrintf(line//achar(13)) 
      write(line,*)'Number of outputs is ',nlhs
      k=mexPrintf(line//achar(13)) 
#endif
      if(nrhs .ne. 15) then
         call mexErrMsgTxt('ERROR: 15 inputs required')
      elseif(nlhs .gt. 3) then
         call mexErrMsgTxt('ERROR: Too many output arguments.')
      endif

C     Check that the inputs are a number.
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check that the input is a number.')
#endif
      do iter=1,nrhs
         if(mxIsNumeric(prhs(iter)) == 0) then
            call mexErrMsgTxt(
     &           'ERROR: Input must be numeric.')
         endif
      end do
            
C     Check data type of input arguments.
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check data type of input argument.')
#endif
      do iter=1,nrhs
         if (mxIsDouble(prhs(iter)) == 0) then
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
      if (ndim<4) then
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

C Finally allocate the input array dimensions  
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate the input array dimensions') 
      call mexWarnMsgTxt('allocate dxin array dimensions'   )       
#endif
      allocate (dxin  (mxGetNumberOfElements(prhs(3))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate dyin array dimensions'   )       
#endif
      allocate (dyin  (mxGetNumberOfElements(prhs(4))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate dz array dimensions'   )       
#endif
      allocate (dz   (mxGetNumberOfElements(prhs(5))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate raxin array dimensions')  
#endif
      allocate (raxin(mxGetNumberOfElements(prhs(6))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate rayin array dimensions'  )
#endif
      allocate (rayin(mxGetNumberOfElements(prhs(7))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate racin array dimensions')   
#endif
      allocate (racin(mxGetNumberOfElements(prhs(8))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mskuin array dimensions' )  
#endif
      allocate (mskuin(mxGetNumberOfElements(prhs(9))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mskvin array dimensions' )    
#endif
      allocate (mskvin(mxGetNumberOfElements(prhs(10))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mskcin array dimensions' )    
#endif
      allocate (mskcin(mxGetNumberOfElements(prhs(11))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate uin array dimensions')   
#endif
      allocate (uin(mxGetNumberOfElements(prhs(12))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate vin array dimensions' )  
#endif
      allocate (vin(mxGetNumberOfElements(prhs(13))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate win array dimensions' )    
#endif
      allocate (win(mxGetNumberOfElements(prhs(14))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate fin array dimensions' )    
#endif
      allocate (fin(mxGetNumberOfElements(prhs(15))))          
          
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
      call mexWarnMsgTxt('allocate uvel array dimensions')    
#endif
      allocate (uvel(dims(1),dims(2),dims(3),dims(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate vvel array dimensions')    
#endif
      allocate (vvel(dims(1),dims(2),dims(3),dims(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate wvel array dimensions')    
#endif
      allocate (wvel(dims(1),dims(2),dims(3),dims(4)))          
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate fld array dimensions')    
#endif
      allocate (fld(dims(1),dims(2),dims(3),dims(4)))                  
          
C Allocate the output from the computational function      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate uflux array dimensions')    
#endif
      allocate (uflux(dims(1),dims(2),dims(3),dims(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate vflux array dimensions')    
#endif
      allocate (vflux(dims(1),dims(2),dims(3),dims(4)))          
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate wflux array dimensions')    
#endif
      allocate (wflux(dims(1),dims(2),dims(3),dims(4)))          
                    
C     Fill Fortran arrays from the input argument pointers.
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill arrays from the input.') 
      call mexWarnMsgTxt('Get iAdvSch from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(1)), iAdvSch, 
     & mxGetNumberOfElements(prhs(1)))
#ifdef DEBUG_MESSAGES
      write(line,*)'iAdvSch: ',iAdvSch
      k=mexPrintf(line//achar(13)) 
      call mexWarnMsgTxt('Get deltat from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), deltat, 
     & mxGetNumberOfElements(prhs(2)))
#ifdef DEBUG_MESSAGES
      write(line,*)'deltat: ',deltat
      k=mexPrintf(line//achar(13)) 
      call mexWarnMsgTxt('Fill dxin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), dxin, 
     & mxGetNumberOfElements(prhs(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill dyin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(4)), dyin, 
     & mxGetNumberOfElements(prhs(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill dz array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(5)), dz, 
     & mxGetNumberOfElements(prhs(5)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill raxin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(6)), raxin, 
     & mxGetNumberOfElements(prhs(6)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill rayin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(7)), rayin, 
     & mxGetNumberOfElements(prhs(7)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill racin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(8)), racin, 
     & mxGetNumberOfElements(prhs(8)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill mskuin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(9)), mskuin, 
     & mxGetNumberOfElements(prhs(9)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill mskvin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(10)), mskvin, 
     & mxGetNumberOfElements(prhs(10)))          
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill mskcin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(11)), mskcin, 
     & mxGetNumberOfElements(prhs(11)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill uin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(12)), uin, 
     & mxGetNumberOfElements(prhs(12)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill vin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(13)), vin, 
     & mxGetNumberOfElements(prhs(13)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill win array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(14)), win, 
     & mxGetNumberOfElements(prhs(14)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill fin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(15)), fin, 
     & mxGetNumberOfElements(prhs(15)))          
          
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

C Check arguments for NANs in the input
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check arguments for NANs')    
      call mexWarnMsgTxt('Check iAdvSch')
#endif
      if( mxIsNaN(iAdvSch) == 1 ) then
C      iAdvSch is pretty important, so error if NaN
	   call mexErrMsgTxt('ERROR: Input iAdvSch contains a NaN.')
      endif
            
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check deltat')
#endif
      if( mxIsNaN(deltat) == 1 ) then
C      deltat is pretty important, so error if NaN
	   call mexErrMsgTxt('ERROR: Input deltat contains a NaN.')
      endif
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check dxin')
#endif
      if (mxIsNaN(dxin) == 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input dxin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(3))
        if( mxIsNaN(dxin(iter)) == 1 ) then
         dxin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check dyin')
#endif
      if (mxIsNaN(dyin) == 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input dyin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(4))
        if( mxIsNaN(dyin(iter)) == 1 ) then
         dyin(iter)=0 .d 0
        endif
       enddo
      endif
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check dz')    
#endif
      if (mxIsNaN(dz) == 1) then
C      dz is pretty important, so error if NaN
       call mexErrMsgTxt('ERROR: Input dz contains a NaN.')
      endif
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check raxin')    
#endif
      if (mxIsNaN(raxin) == 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input raxin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(6))
        if( mxIsNaN(raxin(iter)) == 1 ) then
         raxin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check rayin')    
#endif
       if (mxIsNaN(rayin) == 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input rayin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(7))
        if( mxIsNaN(rayin(iter)) == 1 ) then
         rayin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check racin')    
#endif
      if (mxIsNaN(racin) == 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input racin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(8))
        if( mxIsNaN(racin(iter)) == 1 ) then
         racin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Check mskuin')
#endif    
      if (mxIsNaN(mskuin) == 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input mskuin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(9))
        if( mxIsNaN(mskuin(iter)) == 1 ) then
         mskuin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check mskvin')
#endif
      if (mxIsNaN(mskvin) == 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input mskvin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(10))
        if( mxIsNaN(mskvin(iter)) == 1 ) then
         mskvin(iter)=0 .d 0
        endif
       enddo
      endif
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check mskcin')    
#endif
      if (mxIsNaN(mskcin) == 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input mskcin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(11)) 
        if( mxIsNaN(mskcin(iter)) == 1 ) then
         mskcin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check uin')    
#endif
      if (mxIsNaN(uin) == 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input uin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(12))
        if( mxIsNaN(uin(iter)) == 1 ) then
         uin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check vin')    
#endif
      if (mxIsNaN(vin) == 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input vin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(13))
        if( mxIsNaN(vin(iter)) == 1 ) then
         vin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check win')    
#endif
      if (mxIsNaN(win) == 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input win contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(14))
        if( mxIsNaN(win(iter)) == 1 ) then
         win(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Check fin')
#endif    
      if (mxIsNaN(fin) == 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input fin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(15))
        if( mxIsNaN(fin(iter)) == 1 ) then
         fin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Reshape input arrays') 
#endif   
      dx    = reshape(dxin  ,(/dims(1),dims(2)/))
      dy    = reshape(dyin  ,(/dims(1),dims(2)/))
      rax   = reshape(raxin ,(/dims(1),dims(2),dims(3)/))
      ray   = reshape(rayin ,(/dims(1),dims(2),dims(3)/))
      rac   = reshape(racin ,(/dims(1),dims(2)/))
      msku = reshape(mskuin,(/dims(1),dims(2),dims(3)/))
      mskv = reshape(mskvin,(/dims(1),dims(2),dims(3)/))
      mskc = reshape(mskcin,(/dims(1),dims(2),dims(3)/))
      uvel  = reshape(uin  ,(/dims(1),dims(2),dims(3),dims(4)/))
      vvel  = reshape(vin  ,(/dims(1),dims(2),dims(3),dims(4)/))
      wvel  = reshape(win  ,(/dims(1),dims(2),dims(3),dims(4)/))
      fld   = reshape(fin   ,(/dims(1),dims(2),dims(3),dims(4)/))

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

CC     Call the computational subroutines.

CC----------------ADVECTION IN X DIRECTION----------------
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Calculating Advection in the X direction') 
#endif  
      if (iAdvSch == 1) then
C        first order upwind      
         call advx_uwcntr_2nd(dims(1),dims(2),dims(3),dims(4),
     &           rax,mskc,msku,
     &           uvel,fld,uflux)
      elseif (iAdvSch == 2) then
C        centre second order      
         call advx_cntr_2nd(dims(1),dims(2),dims(3),dims(4),
     &           rax,mskc,msku,
     &           uvel,fld,uflux)
      elseif (iAdvSch == 3) then
C        third order upwind biased
         call advx_uw_3rd(dims(1),dims(2),dims(3),dims(4),
     &           rax,mskc,msku,deltat,
     &           dx,uvel,fld,uflux)
      elseif (iAdvSch == 4) then
C        centre fourth order
         call advx_cntr_4th(dims(1),dims(2),dims(3),dims(4),
     &           rax,mskc,msku,
     &           uvel,fld,uflux)
      elseif (iAdvSch == 7) then
C     7th Order One Step method with Monotonicity Preserving Limiter
         call advx_mpl_7th(dims(1),dims(2),dims(3),dims(4),
     &           rax,mskc,msku,deltat,
     &           dx,uvel,fld,uflux)
      elseif (iAdvSch == 20 .or. iAdvSch == 21) then
C        second order direct-space-and-time, first upwind
         call advx_dst_2nd(dims(1),dims(2),dims(3),dims(4),
     &           iAdvSch,    
     &           rax,mskc,msku,deltat,
     &           dx,uvel,fld,uflux)
      elseif (iAdvSch == 30) then
C        third order direct-space-and-time
         call advx_dst_3rd(dims(1),dims(2),dims(3),dims(4),
     &           rax,mskc,msku,deltat,
     &           dx,uvel,fld,uflux)
      elseif (iAdvSch == 33) then
C        third order direct-space-and-time (flux limited)    
         call advx_dstfl_3rd(dims(1),dims(2),dims(3),dims(4),
     &           rax,mskc,msku,deltat,
     &           dx,uvel,fld,uflux)
      elseif (iAdvSch == 77) then
C        centre second order (flux liimited)    
         call advx_cntrfl_2nd(dims(1),dims(2),dims(3),dims(4),
     &           rax,mskc,msku,deltat,
     &           dx,uvel,fld,uflux)
      else
       call mexErrMsgTxt('ERROR: iAdvSch not recognized, values are:
     &  1,2,3,4,7,20,21,30,33,77')
      endif
     
CC----------------ADVECTION IN Y DIRECTION----------------   
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Calculating Advection in the Y direction') 
#endif     
      if (iAdvSch == 1) then
C        first order upwind     
         call advy_uwcntr_2nd(dims(1),dims(2),dims(3),dims(4),
     &           ray,mskc,mskv,
     &           vvel,fld,vflux)
      elseif (iAdvSch == 2) then
C        centre second order      
         call advy_cntr_2nd(dims(1),dims(2),dims(3),dims(4),
     &           ray,mskc,mskv,
     &           vvel,fld,vflux)
      elseif (iAdvSch == 3) then
C        third order upwind biased
         call advy_uw_3rd(dims(1),dims(2),dims(3),dims(4),
     &           ray,mskc,mskv,deltat,
     &           dy,vvel,fld,vflux)
      elseif (iAdvSch == 4) then
C        centre fourth order
         call advy_cntr_4th(dims(1),dims(2),dims(3),dims(4),
     &           ray,mskc,mskv,
     &           vvel,fld,vflux)
      elseif (iAdvSch == 7) then
C     7th Order One Step method with Monotonicity Preserving Limiter
         call advy_mpl_7th(dims(1),dims(2),dims(3),dims(4),
     &           ray,mskc,mskv,deltat,
     &           dy,vvel,fld,vflux)
      elseif (iAdvSch == 20 .or. iAdvSch == 21) then
C        second order direct-space-and-time, first upwind
         call advy_dst_2nd(dims(1),dims(2),dims(3),dims(4),
     &           iAdvSch,    
     &           ray,mskc,mskv,deltat,
     &           dy,vvel,fld,vflux)
      elseif (iAdvSch == 30) then
C        third order direct-space-and-time
         call advy_dst_3rd(dims(1),dims(2),dims(3),dims(4),
     &           ray,mskc,mskv,deltat,
     &           dy,vvel,fld,vflux)
      elseif (iAdvSch == 33) then
C        third order direct-space-and-time (flux limited)    
         call advy_dstfl_3rd(dims(1),dims(2),dims(3),dims(4),
     &           ray,mskc,mskv,deltat,
     &           dy,vvel,fld,vflux)
      elseif (iAdvSch == 77) then
C        centre second order (flux liimited)    
         call advy_cntrfl_2nd(dims(1),dims(2),dims(3),dims(4),
     &           ray,mskc,mskv,deltat,
     &           dy,vvel,fld,vflux)
      else
       call mexErrMsgTxt('ERROR: iAdvSch not recognized, values are:
     &  1,2,3,4,7,20,21,30,33,77')    
      endif

CC----------------ADVECTION IN Z DIRECTION----------------
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Calculating Advection in the Z direction') 
#endif  
      if (iAdvSch == 1) then
C        first order upwind      
         call advz_uwcntr_2nd(dims(1),dims(2),dims(3),dims(4),
     &           rac,mskc,
     &           wvel,fld,wflux)
      elseif (iAdvSch == 2) then
C        centre second order      
         call advz_cntr_2nd(dims(1),dims(2),dims(3),dims(4),
     &           rac,mskc,
     &           wvel,fld,wflux)
      elseif (iAdvSch == 3) then
C        third order upwind biased
         call advz_uw_3rd(dims(1),dims(2),dims(3),dims(4),
     &           rac,mskc,deltat,dz,
     &           wvel,fld,wflux)
      elseif (iAdvSch == 4) then
C        centre fourth order
         call advz_cntr_4th(dims(1),dims(2),dims(3),dims(4),
     &           rac,mskc,
     &           wvel,fld,wflux)
      elseif (iAdvSch == 7) then
C     7th Order One Step method with Monotonicity Preserving Limiter
         call advz_mpl_7th(dims(1),dims(2),dims(3),dims(4),
     &           rac,mskc,deltat,dz,
     &           wvel,fld,wflux)
      elseif (iAdvSch == 20 .or. iAdvSch == 21) then
C        second order direct-space-and-time, 1st upwind
         call advz_dst_2nd(dims(1),dims(2),dims(3),dims(4),
     &           iAdvSch,    
     &           rac,mskc,deltat,dz,
     &           wvel,fld,wflux)
      elseif (iAdvSch == 30) then
C        third order direct-space-and-time
         call advz_dst_3rd(dims(1),dims(2),dims(3),dims(4),
     &           rac,mskc,deltat,dz,
     &           wvel,fld,wflux)
      elseif (iAdvSch == 33) then
C        third order direct-space-and-time (flux liimited)    
         call advz_dstfl_3rd(dims(1),dims(2),dims(3),dims(4),
     &           rac,mskc,deltat,dz,
     &           wvel,fld,wflux)
      elseif (iAdvSch == 77) then
C        centre second order (flux liimited)    
         call advz_cntrfl_2nd(dims(1),dims(2),dims(3),dims(4),
     &           rac,mskc,deltat,dz,
     &           wvel,fld,wflux)
      else
       call mexErrMsgTxt('ERROR: iAdvSch not recognized, values are:
     &  1,2,3,4,7,20,21,30,33,77')    
      endif
      
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

      subroutine advx_cntr_2nd( imax,jmax,nr,nt,ra,maskc,masku,
     &           uvel,tracer,ut )

C calculates the centered second-order area integrated zonal flux due to advection of a tracer using
      implicit none

      mwIndex i,j,k,t,im1,ip1
      mwIndex imax,jmax,nr,nt
      real*8 maskc(imax,jmax,nr)
      real*8 masku(imax,jmax,nr)
      real*8 ra(imax,jmax,nr)
      real*8 uvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 ut(imax,jmax,nr,nt)
      real*8 half
      parameter(half=1.d0/2.d0)
      
         do t=1,nt
          do k=1,nr
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
             
             ut(i,j,k,t) = uvel(i,j,k,t)*ra(i,j,k)
     &        *(tracer(i,j,k,t)+tracer(im1,j,k,t))*half             

C             ut(i,j,k,t) = maskc(i,j,k)*maskc(im1,j,k )*
C     &                     masku(i,j,k)*masku(im1,j,k )*         
C     &                     uvel(i,j,k,t)*ra(i,j,k)
C     &        *(tracer(i,j,k,t)+tracer(im1,j,k,t))*half
            enddo
           enddo
          enddo
         enddo

      return
      end

      subroutine advx_cntr_4th(imax,jmax,nr,nt,
     &           ra,maskc,masku,uvel,tracer,ut )

C calculates the centered fourth-order area integrated zonal flux due to advection of a tracer
      implicit none

      mwIndex i,j,k,t,im1,im2,ip1
      mwIndex imax,jmax,nr,nt
      real*8 maskc(imax,jmax,nr)
      real*8 masku(imax,jmax,nr)
      real*8 ra(imax,jmax,nr)
      real*8 uvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 ut(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,rjjm,rjjp
      real*8 utrans,half,onesixth
      parameter(half=1.d0/2.d0)
      parameter(onesixth=1.d0/6.d0)
      
      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax           
          if (i .eq. 2) then
           im1=1
           im2=imax
          elseif (i .eq. 1) then
           im1=imax
           im2=imax-1
          else
           im1=i-1
           im2=i-2
          endif

          if (i .eq. imax) then
           ip1=1
          else
           ip1=i+1
          endif     
                    
          rjp = (tracer(ip1,j,k,t)-tracer( i ,j,k,t))
C     &           *masku(ip1,j,k)*maskc(ip1,j,k)
          rj  = (tracer( i ,j,k,t)-tracer(im1,j,k,t))
C     &           *masku( i ,j,k)*maskc( i ,j,k)
          rjm = (tracer(im1,j,k,t)-tracer(im2,j,k,t))
C     &           *masku(im1,j,k)*maskc(im1,j,k)
          rjjp=(rjp-rj)
          rjjm=(rj-rjm)
        
          utrans=uvel(i,j,k,t)*ra(i,j,k)

C          ut(i,j,k,t) = maskc( i ,j,k)*maskc(im1,j,k)*
C     &                 *masku( i ,j,k)*masku(im1,j,k)*     
          ut(i,j,k,t) = 
     &     utrans*(
     &       tracer(i,j,k,t)+tracer(im1,j,k,t)-onesixth*( rjjp+rjjm )
     &             )*half
     &       +abs( utrans )*half*onesixth*( rjjp-rjjm )
     &       *( 1.d0 - masku(im1,j,k)*masku(ip1,j,k) )
         enddo
        enddo
       enddo
      enddo  

      return
      end

      subroutine advx_cntrfl_2nd(imax,jmax,nr,nt,ra, 
     &              maskc,masku,deltat,deltax,uvel,tracer,ut )
     
C calculates the centered second-order area integrated zonal flux due to 
C advection of a tracer  with a flux limiter
      implicit none

      mwIndex i,j,k,t,im1,im2,ip1
      mwIndex imax,jmax,nr,nt
      real*8 deltat
      real*8 deltax(imax,jmax)
      real*8 recipx(imax,jmax)
      real*8 ra(imax,jmax,nr)
      real*8 maskc(imax,jmax,nr)
      real*8 masku(imax,jmax,nr)
      real*8 uvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 ut(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,ucfl,cr,lcr
      real*8 utrans,thetamax,half
      parameter(half=1.d0/2.d0)
      parameter(thetamax = 1.d+20 )

      do j=1,jmax
       do i=1,imax
          recipx(i,j)=1/deltax(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax           
          if (i .eq. 2) then
           im1=1
           im2=imax
          elseif (i .eq. 1) then
           im1=imax
           im2=imax-1
          else
           im1=i-1
           im2=i-2
          endif
          
          rjp=(tracer(im1,j,k,t)-tracer( i ,j,k,t))
C     &     *masku(im1,j,k)*maskc(im1,j,k)
          rj =(tracer( i ,j,k,t)-tracer(im1,j,k,t))
C     &     *masku( i ,j,k)*maskc( i ,j,k)
          rjm=(tracer(im1,j,k,t)-tracer(im2,j,k,t))
C     &     *masku(im1,j,k)*maskc(im1,j,k)

          ucfl = abs( uvel(i,j,k,t)*deltat*recipx(i,j) )
          utrans=uvel(i,j,k,t)*ra(i,j,k)

        if (rj.ne.0.) then
         if (utrans.gt.0) then
           cr=rjm/rj
         else
           cr=rjp/rj
         endif
        else
         if (utrans.gt.0) then
           cr=rjm*thetamax
         else
           cr=rjp*thetamax
         endif
        endif
        
C The flux limiter can be -
C upwind:
C     limiter(cr)=0.
C lax-wendroff:
C     limiter(cr)=1.
C min-mod:
C     limiter(cr)=max(0.d0,min(1.d0,cr))
C suberbee:
        lcr=max(0.d0,max(min(1.d0,2.d0*cr),
     &                         min(2.d0,cr)))
        
C        ut(i,j,k,t) = masku( i ,j,k)*masku(im1,j,k)*
C     &                maskc( i ,j,k)*maskc(im1,j,k)* 
        ut(i,j,k,t) =   
     &   utrans*(tracer(i,j,k,t)+tracer(im1,j,k,t))*half
     &   -abs(utrans)*((1.d0-lcr)+ucfl*lcr)
     &                    *rj*half
         enddo
        enddo
       enddo
      enddo  

      return
      end

      subroutine advx_dst_2nd(
     &           imax,jmax,nr,nt,advectionscheme,ra,
     &           maskc,masku,deltat,deltax,uvel,tracer,ut )

C  calculates the second-order direct space and time (dst-2/lax-wendroff, or simple 
C  first order upwind scheme) area integrated zonal flux due to advection of a tracer 
      implicit none

      mwIndex i,j,k,t,im1
      mwIndex imax,jmax,nr,nt
      real*8 advectionscheme, deltat
      real*8 ucfl, xlimit, uabs, utrans
      real*8 deltax(imax,jmax)
      real*8 recipx(imax,jmax)
      real*8 ra    (imax,jmax,nr)
      real*8 maskc (imax,jmax,nr)
      real*8 masku (imax,jmax,nr)
      real*8 uvel  (imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 ut    (imax,jmax,nr,nt)
      
      if ( advectionscheme.eq.20 ) then
       xlimit = 1.d0
      else
       xlimit = 0.d0
      endif 

      do j=1,jmax
       do i=1,imax
          recipx(i,j)=1/deltax(i,j)
       enddo
      enddo
         
	  do t=1,nt		  
	   do k=1,nr		
	    do j=1,jmax
         do i=1,imax
		  
		  if (i .eq. 1) then
		   im1=imax
		  else
		   im1=i-1
		  endif
	  
		  ucfl = abs( uvel(i,j,k,t)*deltat*recipx(i,j) )

		  utrans=uvel(i,j,k,t)*ra(i,j,k)
		
		  uabs = abs(utrans)*( 1.d0 - xlimit*(1.d0 - ucfl) )
		  
		  ut(i,j,k,t) = 
     &		  ( utrans+uabs )* 0.5d0 * tracer(im1,j,k,t)
     &		+ ( utrans-uabs )* 0.5d0 * tracer(i  ,j,k,t)
	     enddo
	    enddo
	   enddo
	  enddo
      
      return
      end

      subroutine advx_dst_3rd( imax,jmax,nr,nt,ra, 
     &              maskc,masku,deltat,deltax,uvel,tracer,ut )
   
C  calculates the 3rd-order direct space and time (dst-3) area integrated zonal flux 
C  due to advection of a tracer.
      implicit none

      mwIndex i,j,k,t,im1,im2
      mwIndex imax,jmax,nr,nt
      real*8 deltat
      real*8 deltax(imax,jmax)
      real*8 recipx(imax,jmax)
      real*8 ra(imax,jmax,nr)
      real*8 maskc(imax,jmax,nr)
      real*8 masku(imax,jmax,nr)
      real*8 uvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 ut(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,d0,d1,onesixth
      real*8 utrans,ucfl,half
      parameter(half=1.d0/2.d0)
      parameter(onesixth=1.d0/6.d0)
      
      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          ut(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        
      
      do j=1,jmax
       do i=1,imax
          recipx(i,j)=1/deltax(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax           
          if (i .eq. 2) then
           im1=1
           im2=imax
          elseif (i .eq. 1) then
           im1=imax
           im2=imax-1
          else
           im1=i-1
           im2=i-2
          endif
             
          rjp=(tracer(im1,j,k,t)-tracer( i ,j,k,t))
     &     *maskc( i ,j,k)*maskc(im1,j,k)
          rj =(tracer( i ,j,k,t)-tracer(im1,j,k,t))
     &     *masku( i ,j,k)*maskc( i ,j,k)
          rjm=(tracer(im1,j,k,t)-tracer(im2,j,k,t))
     &     *maskc(im1,j,k)*maskc(im2,j,k)

          ucfl = abs( uvel(i,j,k,t)*deltat*recipx(i,j) )
          utrans=uvel(i,j,k,t)*ra(i,j,k)
          
          d0=(2.d0-ucfl)*(1.d0-ucfl)*onesixth
          d1=(1.d0-ucfl*ucfl)*onesixth

          ut(i,j,k,t)=maskc( i ,j,k)*maskc(im1,j,k)
C     &              *masku( i ,j,k)*masku(im1,j,k)
C          ut(i,j,k,t)= (
     &   *(half*(utrans+abs(utrans))
     &      *( tracer(im1,j,k,t) + (d0*rj+d1*rjm) )
     &    +half*(utrans-abs(utrans))
     &      *( tracer( i ,j,k,t) - (d0*rj+d1*rjp) ))

         enddo
        enddo
       enddo
      enddo  

      return
      end

      subroutine advx_dstfl_3rd( imax,jmax,nr,nt,ra, 
     &              maskc,masku,deltat,deltax,uvel,tracer,ut )
     
C  calculates the 3rd-order direct space and time (dst-3) area integrated zonal flux 
C  due to advection of a tracer with flux limiting.
      implicit none

      mwIndex i,j,k,t,im1,im2,ip1
      mwIndex imax,jmax,nr,nt
      real*8 deltat
      real*8 deltax(imax,jmax)
      real*8 recipx(imax,jmax)
      real*8 ra(imax,jmax,nr)
      real*8 maskc(imax,jmax,nr)
      real*8 masku(imax,jmax,nr)
      real*8 uvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 ut(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,d0,d1,psip,psim
      real*8 thetamin,thetamax,thetap,thetam
      real*8 utrans,ucfl,onesixth,half
      parameter(half=1.d0/2.d0)
      parameter(thetamax = 1.d+20 )
      parameter(thetamin = 1.d-20 )
      parameter(onesixth=1.d0/6.d0)

      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          ut(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do j=1,jmax
       do i=1,imax
          recipx(i,j)=1/deltax(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax           
          if (i .eq. 2) then
           im1=1
           im2=imax
          elseif (i .eq. 1) then
           im1=imax
           im2=imax-1
          else
           im1=i-1
           im2=i-2
          endif
          
          rjp=(tracer(im1,j,k,t)-tracer( i ,j,k,t))
     &     *maskc( i ,j,k)*maskc(im1,j,k)
          rj =(tracer( i ,j,k,t)-tracer(im1,j,k,t))
     &     *maskc( i ,j,k)*maskc(im1,j,k)
          rjm=(tracer(im1,j,k,t)-tracer(im2,j,k,t))
     &     *maskc(im1,j,k)*maskc(im2,j,k)

          ucfl = abs( uvel(i,j,k,t)*deltat*recipx(i,j) )
          utrans=uvel(i,j,k,t)*ra(i,j,k)
          
          d0=(2.d0-ucfl)*(1.d0-ucfl)*onesixth
          d1=(1.d0-ucfl*ucfl)*onesixth
          
C      prevent |thetap,m| to reach too big value
           if ( abs(rj)*thetamax .le. abs(rjm) ) then
             thetap=sign(thetamax,rjm*rj)
           else
             thetap=rjm/rj
           endif
           if ( abs(rj)*thetamax .le. abs(rjp) ) then
             thetam=sign(thetamax,rjp*rj)
           else
             thetam=rjp/rj
           endif

           psip=d0+d1*thetap
           psip=max(0.d0,min(min(1.d0,psip),
     &          thetap*(1.d0 -ucfl)/(ucfl+thetamin) ))
           psim=d0+d1*thetam
           psim=max(0.d0,min(min(1.d0,psim),
     &          thetam*(1.d0 -ucfl)/(ucfl+thetamin) ))

           ut(i,j,k,t)=maskc( i ,j,k)*maskc(im1,j,k)
C     &                *masku( i ,j,k)*masku(im1,j,k)
C           ut(i,j,k,t)=(
     &   *(half*(utrans+abs(utrans))
     &      *( tracer(i,j, k ,t) + psim*rj )
     &  +half*(utrans-abs(utrans))
     &      *( tracer(im1,j,k,t) - psip*rj ))
         enddo
        enddo
       enddo
      enddo  

      return
      end

      subroutine advx_mpl_7th(imax,jmax,nr,nt,ra, 
     &           maskc,masku,deltat,deltax,uvel,tracer,ut)

C calculates the 7th order dst area integrated zonal flux due to 
C advection of a tracer with a monotone preserving limiter
      implicit none

      mwIndex i,j,k,t,imax,jmax,nr,nt
      mwIndex im1,im2,im3,im4,ip1,ip2,ip3
      real*8 deltat
      real*8 deltax(imax,jmax)
      real*8 recipx(imax,jmax)
      real*8 ra(imax,jmax,nr)
      real*8 maskc(imax,jmax,nr)
      real*8 masku(imax,jmax,nr)
      real*8 uvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 ut(imax,jmax,nr,nt)
      real*8 cfl,psi,utrans
      real*8 uloc,fac,delip,deli,phi,eps,rp1h,rp1h_cfl
      real*8 recip_delip, recip_deli
      real*8 qippp,qipp,qip,qi,qim,qimm,qimmm
      real*8 mskipp,mskip,mski,mskim,mskimm,mskimmm
      real*8 d2,d2p1,d2m1,a,b,c,d
      real*8 dp1h,dm1h, phimd,philc,phimin,phimax
      real*8 delm,delp,delmm,delpp,delmmm,delppp
      real*8 del2mm,del2m,del2,del2p,del2pp
      real*8 del3mm,del3m,del3p,del3pp
      real*8 del4m,del4,del4p
      real*8 del5m,del5p
      real*8 del6
      parameter(eps = 1.d-20)
      
      
      do j=1,jmax
       do i=1,imax
          recipx(i,j)=1/deltax(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax        
      
          if (i .eq. 4) then
           im1=i-1
           im2=1-2
           im3=1-3
           im4=imax
          elseif (i .eq. 3) then
           im1=i-1
           im2=i-2
           im3=imax
           im4=imax-1          
          elseif (i .eq. 2) then
           im1=1
           im2=imax
           im3=imax-1
           im4=imax-2
          elseif (i .eq. 1) then
           im1=imax
           im2=imax-1 
           im3=imax-2
           im4=imax-3
          else
           im1=i-1
           im2=i-2
           im3=i-3
           im4=i-4
          endif
          
          if (i .eq. imax ) then
           ip1=1
           ip2=2
           ip3=3
          elseif (i .eq. imax-1 ) then
           ip1=i+1
           ip2=1
           ip3=2
          elseif (i .eq. imax-2 ) then
           ip1=i+1
           ip2=i+2
           ip3=1
          else
           ip1=i+1
           ip2=i+2
           ip3=i+3   
          endif
          
          cfl = abs(uvel(i,j,k,t)*deltat*recipx(i,j) )
          utrans=uvel(i,j,k,t)*ra(i,j,k)
        
          if (utrans.gt.0.d0) then
           qippp = tracer(ip2,j,k,t)
           qipp  = tracer(ip1,j,k,t)
           qip   = tracer(i  ,j,k,t)
           qi    = tracer(im1,j,k,t)
           qim   = tracer(im2,j,k,t)
           qimm  = tracer(im3,j,k,t)
           qimmm = tracer(im4,j,k,t)

           mskipp  = masku(ip2,j,k)
           mskip   = masku(ip1,j,k)
           mski    = masku(i  ,j,k)
           mskim   = masku(im1,j,k)
           mskimm  = masku(im2,j,k)
           mskimmm = masku(im3,j,k)
          elseif (utrans.lt.0.d0) then
           qippp = tracer(im3,j,k,t)
           qipp  = tracer(im2,j,k,t)
           qip   = tracer(im1,j,k,t)
           qi    = tracer(i  ,j,k,t)
           qim   = tracer(ip1,j,k,t)
           qimm  = tracer(ip2,j,k,t)
           qimmm = tracer(ip3,j,k,t)

           mskipp  = masku(im2,j,k)
           mskip   = masku(im1,j,k)
           mski    = masku(i  ,j,k)
           mskim   = masku(ip1,j,k)
           mskimm  = masku(ip2,j,k)
           mskimmm = masku(ip3,j,k)
          else
           qippp = 0.d0
           qipp  = 0.d0
           qip   = 0.d0
           qi    = 0.d0
           qim   = 0.d0
           qimm  = 0.d0
           qimmm = 0.d0

           mskipp  = 0.d0
           mskip   = 0.d0
           mski    = 0.d0
           mskim   = 0.d0
           mskimm  = 0.d0
           mskimmm = 0.d0
          endif

          if (utrans.ne.0.d0) then
C          2nd order correction [i im1]
           fac = 1.d0
           delp = (qip-qi)*mski
           phi = fac * delp
C          3rd order correction [i im1 im2]
           fac = fac * ( cfl + 1.d0 )/3.d0
           delm = (qi-qim)*mskim
           del2 = delp - delm
           phi = phi - fac * del2
C          4th order correction [ip1 i im1 im2]
           fac = fac * ( cfl - 2.d0 )/4.d0
           delpp = (qipp-qip)*mskip*mski
           del2p = delpp - delp
           del3p = del2p - del2
           phi = phi + fac * del3p
C          5th order correction [ip1 i im1 im2 im3]
           fac = fac * ( cfl - 3.d0 )/5.d0
           delmm = (qim-qimm)*mskimm*mskim
           del2m = delm - delmm
           del3m = del2 - del2m
           del4 = del3p - del3m
           phi = phi + fac * del4
C          6th order correction [ip2 ip1 i im1 im2 im3]
           fac = fac * ( cfl + 2.d0 )/6.d0
           delppp = (qippp-qipp)*mskipp*mskip*mski
           del2pp = delpp - delp
           del3pp = del2pp - del2p
           del4p = del3pp - del3p
           del5p = del4p - del4
           phi = phi + fac * del5p
C          7th order correction [ip2 ip1 i im1 im2 im3 im4]
           fac = fac * ( cfl + 2.d0 )/7.d0
           delmmm = (qimm-qimmm)*mskimmm*mskimm*mskim
           del2mm = delmm - delmmm
           del3mm = del2m - del2mm
           del4m = del3m - del3mm
           del5m = del4 - del4m
           del6 = del5p - del5m
           phi = phi - fac * del6

           delip = ( qip - qi ) * mski
C        phi = sign(1.d0,phi)*sign(1.d0,delip)
C    &        *abs(phi+eps)/abs(delip+eps)
C     simplify and avoid division by zero
           recip_delip = sign(1.d0,delip)/max(abs(delip),eps)
           phi = phi*recip_delip

           deli = ( qi - qim ) * mskim
C        rp1h =sign(1.d0,deli)*sign(1.d0,delip)
C    &        *abs(deli+eps)/abs(delip+eps)
C     simplify and avoid division by zero
           recip_deli = sign(1.d0,deli)/max(abs(deli),eps)
           rp1h = deli*recip_delip
           rp1h_cfl = rp1h/(cfl+eps)

C        tvd limiter
C        phi = max(0.d0, min( 2./(1-cfl), phi, 2.*rp1h_cfl ) )

C        mp limiter
           d2   = del2 !( ( qip + qim ) - 2.*qi  ) * mski * mskim
           d2p1 = del2p !( ( qipp + qi ) - 2.*qip ) * mskip * mski
           d2m1 = del2m !( ( qi + qimm ) - 2.*qim ) * mskim * mskimm
           a = 4.d0*d2 - d2p1
           b = 4.d0*d2p1 - d2
           c = d2
           d = d2p1
           dp1h = max(min(a,b,c,d),0.d0)+min(max(a,b,c,d),0.d0)
           a = 4.d0*d2m1 - d2
           b = 4.d0*d2 - d2m1
           c = d2m1
           d = d2
           dm1h = max(min(a,b,c,d),0.d0)+min(max(a,b,c,d),0.d0)
C        qmd = 0.5*( ( qi + qip ) - dp1h )
C        qmd = 0.5d0*( ( 2.d0*qi + delip ) - dp1h )
C        qul = qi + (1.d0-cfl)/(cfl+eps)*deli
C        qlc = qi + 0.5d0*( 1.d0+dm1h/(deli+eps) )*(qul-qi)
C        phimd = 2.d0/(1.d0-cfl)*(qmd-qi+eps)/(delip+eps)
C        philc = 2.d0*rp1h_cfl*(qlc-qi+eps)/(qul-qi+eps)
C     simplify and avoid division by zero
           phimd = 1.d0/(1.d0-cfl)*(delip-dp1h)*recip_delip
           philc = rp1h_cfl*( 1.d0+dm1h*recip_deli )
C  
           phimin = max( min(0.d0,phimd),
     &                 min(0.d0,2.d0*rp1h_cfl,philc) )
           phimax = min( max(2.d0/(1.d0-cfl),phimd),
     &                 max(0.d0,2.d0*rp1h_cfl,philc) )
           phi = max(phimin,min(phi,phimax))

           psi = phi * 0.5d0 * (1.d0 - cfl)
           ut(i,j,k,t) = utrans*( qi + psi*delip )
          else
           ut(i,j,k,t) = 0.d0
          endif
         enddo 
        enddo
       enddo
      enddo

      return
      end

      subroutine advx_uw_3rd(imax,jmax,nr,nt,ra, 
     &              maskc,masku,deltat,deltax,uvel,tracer,ut )

C calculates the third-order upwind biased area integrated zonal flux due to 
C advection of a tracer
      implicit none

      mwIndex i,j,k,t,im1,im2,ip1
      mwIndex imax,jmax,nr,nt
      real*8 deltat
      real*8 deltax(imax,jmax)
      real*8 ra(imax,jmax,nr)
      real*8 masku(imax,jmax,nr)
      real*8 maskc(imax,jmax,nr)
      real*8 uvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 ut(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,rjjp,rjjm
      real*8 utrans,ucfl,onesixth,half
      parameter(half=1.d0/2.d0)
      parameter(onesixth=1.d0/6.d0)

      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          ut(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax           
          if (i .eq. 2) then
           im1=1
           im2=imax
          elseif (i .eq. 1) then
           im1=imax
           im2=imax-1
          else
           im1=i-1
           im2=i-2
          endif
          
          if (i .eq. imax) then
           ip1=1
          else
           ip1=i+1
          endif          
          
          rjp=(tracer(ip1,j,k,t)-tracer( i ,j,k,t))
     &     *maskc(ip1,j,k)*maskc( i ,j,k)
          rj =(tracer( i ,j,k,t)-tracer(im1,j,k,t))
     &     *maskc( i ,j,k)*maskc(im1,j,k)
          rjm=(tracer(im1,j,k,t)-tracer(im2,j,k,t))
     &     *maskc(im1,j,k)*maskc(im2,j,k)

          rjjp=rjp-rj
          rjjm=rj-rjm
        
          utrans=uvel(i,j,k,t)*ra(i,j,k)
        
          ut(i,j,k,t) = maskc( i ,j,k)*
C     &                  maskc(im1,j,k)*
C     &                  masku( i ,j,k)*masku(im1,j,k)*     
C          ut(i,j,k,t) =
     &   utrans*(
     &     tracer(i,j,k,t)+tracer(im1,j,k,t)-onesixth*( rjjp+rjjm )
     &               )*half
     &     +abs( utrans )*half*onesixth*( rjjp-rjjm )
         enddo
        enddo
       enddo
      enddo  

      return
      end

      subroutine advx_uwcntr_2nd( imax,jmax,nr,nt,ra,maskc,
     &           masku,uvel,tracer,ut )

C calculates the center second upwind biased area integrated zonal flux due to 
C advection of a tracer
      implicit none

      mwIndex i,j,k,t,im1
      mwIndex imax,jmax,nr,nt
      real*8 maskc(imax,jmax,nr)
      real*8 masku(imax,jmax,nr)
      real*8 ra(imax,jmax,nr)
      real*8 uvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 ut(imax,jmax,nr,nt)
      
         do t=1,nt
          do k=1,nr
           do j=1,jmax
            do i=1,imax
             if (i .eq. 1) then
                im1=imax
             else
                im1=i-1
             endif
             
             if (uvel(i,j,k,t) .ge. 0.d0) then
C                ut(i,j,k,t) = masku(i,j,k)*masku(im1,j,k )*
C     &                        maskc(i,j,k)*maskc(im1,j,k )* 
                ut(i,j,k,t) = 
     &                        uvel(i,j,k,t)*tracer(im1,j,k,t)*ra(i,j,k)
             else
C                ut(i,j,k,t) = masku(i,j,k )*maskc(i,j,k)*
                ut(i,j,k,t) = 
     &                        uvel(i,j,k,t)*tracer( i ,j,k,t)*ra(i,j,k)
             endif

            enddo
           enddo
          enddo
         enddo

      return
      end

      subroutine advy_cntr_2nd( imax,jmax,nr,nt,ra,maskc,maskv,
     &           vvel,tracer,vt )

C calculates the centered second-order area integrated meridional flux due to advection of a tracer
      implicit none

      mwIndex i,j,k,t,jm1,jp1
      mwIndex imax,jmax,nr,nt
      real*8 maskc(imax,jmax,nr)
      real*8 maskv(imax,jmax,nr)
      real*8 ra(imax,jmax,nr)
      real*8 vvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 vt(imax,jmax,nr,nt)
      real*8 half
      parameter(half=1.d0/2.d0)      
      
         do t=1,nt
          do k=1,nr
           do j=1,jmax
            do i=1,imax
             jm1=max(j-1, 1  )
             jp1=min(j+1,jmax)
             
             vt(i,j,k,t) = vvel(i,j,k,t)*ra(i,j,k)
     &             *(tracer(i,j,k,t)+tracer(i,jm1,k,t))*half

C             vt(i,j,k,t) = maskc(i,j,k)*maskc(i,jm1,k)*
C     &                     maskv(i,j,k)*maskv(i,jm1,k)*
C     &                     vvel(i,j,k,t)*ra(i,j,k)
C     &             *(tracer(i,j,k,t)+tracer(i,jm1,k,t))*half
            enddo
           enddo  
          enddo
         enddo

      return
      end

      subroutine advy_cntr_4th(imax,jmax,nr,nt,
     &           ra,maskc,maskv,vvel,tracer,vt)

C calculates the centered fourth-order area integrated meridional flux due to advection of a tracer
      implicit none
      
      mwIndex i,j,k,t,jm1,jm2,jp1
      mwIndex imax,jmax,nr,nt
      real*8 maskv(imax,jmax,nr)
      real*8 maskc(imax,jmax,nr)
      real*8 ra(imax,jmax,nr)
      real*8 vvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 vt(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,rjjm,rjjp
      real*8 vtrans,half,onesixth
      parameter(half=1.d0/2.d0)
      parameter(onesixth=1.d0/6.d0)
      
      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax           
          jp1=min(j+1,jmax)
          jm1=max(j-1, 1  )
          jm2=max(j-2, 1  )
       
          rjp = (tracer(i,jp1,k,t)-tracer(i, j ,k,t))
C     &           *maskv(i,jp1,k)*maskc(i,jp1,k)
          rj  = (tracer(i, j ,k,t)-tracer(i,jm1,k,t))
C     &           *maskv(i, j ,k)*maskc(i, j ,k)
          rjm = (tracer(i,jm1,k,t)-tracer(i,jm2,k,t))
C     &           *maskv(i,jm1,k)*maskc(i,jm1,k)
          rjjp=(rjp-rj)
          rjjm=(rj-rjm)
          
          vtrans=vvel(i,j,k,t)*ra(i,j,k)

C          vt(i,j,k,t) = maskc(i, j ,k)*maskc(i,jm1,k)*
C     &                 *maskv(i, j ,k)*maskv(i,jm1,k)*     
          vt(i,j,k,t) = 
     &     vtrans*(
     &       tracer(i,j,k,t)+tracer(i,jm1,k,t)-onesixth*( rjjp+rjjm )
     &            )*half
     &       +abs( vtrans )*half*onesixth*( rjjp-rjjm )
     &       *( 1.d0 - maskv(i,jm1,k)*maskv(i,jp1,k) )
         enddo
        enddo
       enddo
      enddo  

      return
      end

      subroutine advy_cntrfl_2nd(imax,jmax,nr,nt,ra, 
     &              maskc,maskv,deltat,deltay,vvel,tracer,vt )

C calculates the centered second-order area integrated meridional flux due to 
C advection of a tracer with a flux limiter
      implicit none
      
      mwIndex i,j,k,t,jm1,jm2,jp1
      mwIndex imax,jmax,nr,nt
      real*8 deltat
      real*8 deltay(imax,jmax)
      real*8 recipy(imax,jmax)
      real*8 ra(imax,jmax,nr)
      real*8 maskc(imax,jmax,nr)
      real*8 maskv(imax,jmax,nr)
      real*8 vvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 vt(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,vcfl,cr,lcr
      real*8 vtrans,thetamax,half
      parameter(half=1.d0/2.d0)
      parameter(thetamax = 1.d+20 )

      do j=1,jmax
       do i=1,imax
          recipy(i,j)=1/deltay(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax           
          jp1=min(j+1,jmax)
          jm1=max(j-1, 1  )
          jm2=max(j-2, 1  )
          
          rjp=(tracer(i,jp1,k,t)-tracer(i, j ,k,t))
C     &     *maskv(i,jp1,k)*maskc(i,jp1,k)
          rj =(tracer(i, j ,k,t)-tracer(i,jm1,k,t))
C     &     *maskv(i, j ,k)*maskc(i, j ,k)
          rjm=(tracer(i,jm1,k,t)-tracer(i,jm2,k,t))
C     &     *maskv(i,jm1,k)*maskc(i,jm1,k)

          vcfl = abs( vvel(i,j,k,t)*deltat*recipy(i,j) )
          vtrans=vvel(i,j,k,t)*ra(i,j,k)

          if (rj.ne.0.) then
           if (vtrans.gt.0) then
            cr=rjm/rj
          else
            cr=rjp/rj
          endif
         else
          if (vtrans.gt.0) then
            cr=rjm*thetamax
          else
            cr=rjp*thetamax
          endif
         endif

C the flux limiter can be -
C upwind:
C     limiter(cr)=0.
C lax-wendroff:
C     limiter(cr)=1.
C min-mod:
C     limiter(cr)=max(0.d0,min(1.d0,cr))
C suberbee:
        lcr=max(0.d0,max(min(1.d0,2.d0*cr),
     &                         min(2.d0,cr)))
             
C         vt(i,j,k,t) =maskv(i, j ,k)*maskv(i,jm1,k)*
C     &                maskc(i, j ,k)*maskc(i,jm1,k)* 
         vt(i,j,k,t) =   
     &   vtrans*(tracer(i,j,k,t)+tracer(i,jm1,k,t))*half
     &   -abs(vtrans)*((1.d0-lcr)+vcfl*lcr)
     &                    *rj*half
         enddo
        enddo
       enddo
      enddo  

      return
      end

      subroutine advy_dst_2nd(
     &           imax,jmax,nr,nt,advectionscheme,ra,
     &           maskc,maskv,deltat,deltay,vvel,tracer,vt)

C  calculates the second-order direct space and time (dst-2/lax-wendroff, or simple 
C  first order upwind scheme) area integrated meridional flux due to advection of a tracer 
      implicit none

      mwIndex i,j,k,t,jm1
      mwIndex imax,jmax,nr,nt
      real*8 advectionscheme, deltat
      real*8 vcfl, ylimit, vabs, vtrans
      real*8 deltay(imax,jmax)
      real*8 recipy(imax,jmax)
      real*8 ra    (imax,jmax,nr)
      real*8 maskc (imax,jmax,nr)
      real*8 maskv (imax,jmax,nr)
      real*8 vvel  (imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 vt    (imax,jmax,nr,nt)

      if ( advectionscheme.eq.20 ) then
       ylimit = 1.d0
      else
       ylimit = 0.d0
      endif
    
      do j=1,jmax
       do i=1,imax
          recipy(i,j)=1/deltay(i,j)
       enddo
      enddo
      
      do t=1,nt
       do k=1,nr
        do j=1,jmax
         
         jm1=max(j-1, 1  )
         
         do i=1,imax
          
          vcfl = abs( vvel(i,j,k,t)*deltat*recipy(i,j) )
          vtrans=vvel(i,j,k,t)*ra(i,j,k)

          vabs = abs(vtrans)*( 1.d0 - ylimit*(1.d0 - vcfl) )
          vt(i,j,k,t) = ( vtrans+vabs )* 0.5d0 * tracer(i,jm1,k,t)
     &                + ( vtrans-vabs )* 0.5d0 * tracer(i,j  ,k,t)

         enddo
        enddo
       enddo
      enddo
      
      return
      end

      subroutine advy_dst_3rd( imax,jmax,nr,nt,ra, 
     &              maskc,maskv,deltat,deltay,vvel,tracer,vt )

C  calculates the 3rd-order direct space and time (dst-3) area integrated meridional flux 
C  due to advection of a tracer.
      implicit none

      mwIndex i,j,k,t,jm1,jm2,jp1
      mwIndex imax,jmax,nr,nt
      real*8 deltat
      real*8 deltay(imax,jmax)
      real*8 recipy(imax,jmax)
      real*8 ra(imax,jmax,nr)
      real*8 maskc(imax,jmax,nr)
      real*8 maskv(imax,jmax,nr)
      real*8 vvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 vt(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,d0,d1,onesixth
      real*8 vtrans,vcfl,half
      parameter(half=1.d0/2.d0)
      parameter(onesixth=1.d0/6.d0)

      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          vt(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do j=1,jmax
       do i=1,imax
          recipy(i,j)=1/deltay(i,j)
       enddo
      enddo
            
      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax           
          jp1=min(j+1,jmax)
          jm1=max(j-1, 1  )
          jm2=max(j-2, 1  )
          
          rjp=(tracer(i,jp1,k,t)-tracer(i, j ,k,t))
     &     *maskc(i, j ,k)*maskc(i,jp1,k)
          rj =(tracer(i, j ,k,t)-tracer(i,jm1,k,t))
     &     *maskc(i, j ,k)*maskc(i,jm1,k)
          rjm=(tracer(i,jm1,k,t)-tracer(i,jm2,k,t))
     &     *maskc(i,jm1,k)*maskc(i,jm2,k)

          vcfl = abs( vvel(i,j,k,t)*deltat*recipy(i,j) )
          vtrans=vvel(i,j,k,t)*ra(i,j,k)
          
          d0=(2.d0-vcfl)*(1.d0-vcfl)*onesixth
          d1=(1.d0-vcfl*vcfl)*onesixth

          vt(i,j,k,t)=maskc(i, j ,k)*maskc(i,jm1,k)
C     &               *maskv(i, j ,k)*maskv(i,jm1,k)   
C          vt(i,j,k,t)= (  
     &   *( half*(vtrans+abs(vtrans))
     &      *( tracer(i,jm1,k,t) + (d0*rj+d1*rjm) )
     &  +half*(vtrans-abs(vtrans))
     &      *( tracer(i, j ,k,t) - (d0*rj+d1*rjp) ))
         enddo
        enddo
       enddo
      enddo  
      
      return
      end

      subroutine advy_dstfl_3rd( imax,jmax,nr,nt,ra, 
     &              maskc,maskv,deltat,deltay,vvel,tracer,vt )

C  calculates the 3rd-order direct space and time (dst-3) area integrated zonal flux 
C  due to advection of a tracer with flux limiting.
      implicit none

      mwIndex i,j,k,t,jm1,jm2,jp1
      mwIndex imax,jmax,nr,nt
      real*8 deltat
      real*8 deltay(imax,jmax)
      real*8 recipy(imax,jmax)
      real*8 ra(imax,jmax,nr)
      real*8 maskc(imax,jmax,nr)
      real*8 maskv(imax,jmax,nr)
      real*8 vvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 vt(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,d0,d1,psip,psim
      real*8 thetamin,thetamax,thetap,thetam
      real*8 vtrans,vcfl,onesixth,half
      parameter(half=1.d0/2.d0)
      parameter(thetamax = 1.d+20 )
      parameter(thetamin = 1.d-20 )
      parameter(onesixth=1.d0/6.d0)

      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          vt(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do j=1,jmax
       do i=1,imax
          recipy(i,j)=1/deltay(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax           
          jp1=min(j+1,jmax)
          jm1=max(j-1, 1  )
          jm2=max(j-2, 1  )
          
          rjp=(tracer(i,jp1,k,t)-tracer(i, j ,k,t))
     &     *maskc(i, j ,k)*maskc(i,jp1,k)
          rj =(tracer(i, j ,k,t)-tracer(i,jm1,k,t))
     &     *maskc(i, j ,k)*maskc(i,jm1,k)
          rjm=(tracer(i,jm1,k,t)-tracer(i,jm2,k,t))
     &     *maskc(i,jm1,k)*maskc(i,jm2,k)

          vcfl = abs( vvel(i,j,k,t)*deltat*recipy(i,j) )
          vtrans=vvel(i,j,k,t)*ra(i,j,k)
          
          d0=(2.d0-vcfl)*(1.d0-vcfl)*onesixth
          d1=(1.d0-vcfl*vcfl)*onesixth

C      prevent |thetap,m| to reach too big value
           if ( abs(rj)*thetamax .le. abs(rjm) ) then
             thetap=sign(thetamax,rjm*rj)
           else
             thetap=rjm/rj
           endif
           if ( abs(rj)*thetamax .le. abs(rjp) ) then
             thetam=sign(thetamax,rjp*rj)
           else
             thetam=rjp/rj
           endif

           psip=d0+d1*thetap
           psip=max(0.d0,min(min(1.d0,psip),
     &          thetap*(1.d0 -vcfl)/(vcfl+thetamin) ))
           psim=d0+d1*thetam
           psim=max(0.d0,min(min(1.d0,psim),
     &          thetam*(1.d0 -vcfl)/(vcfl+thetamin) ))

           vt(i,j,k,t)=maskc(i, j ,k)*maskc(i,jm1,k)
C     &                *maskv(i, j ,k)*maskv(i,jm1,k)
C           vt(i,j,k,t)=(      
     &   *( half*(vtrans+abs(vtrans))
     &      *( tracer(i, j ,k,t) + psim*rj )
     &  +half*(vtrans-abs(vtrans))
     &      *( tracer(i,jm1,k,t) - psip*rj ))
         enddo
        enddo
       enddo
      enddo  
      
      return
      end

      subroutine advy_mpl_7th(imax,jmax,nr,nt,ra, 
     &              maskc,maskv,deltat,deltay,vvel,tracer,vt )

C calculates the 7th order dst area integrated meridional flux due to 
C advection of a tracer with a monotone preserving limiter
      implicit none

      mwIndex imax,jmax,nr,nt,i,j,k,t
      mwIndex jm1,jm2,jm3,jm4,jp1,jp2,jp3
      real*8 deltat
      real*8 deltay(imax,jmax)
      real*8 recipy(imax,jmax)
      real*8 ra    (imax,jmax,nr)
      real*8 maskc (imax,jmax,nr)
      real*8 maskv (imax,jmax,nr)
      real*8 vvel  (imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 vt    (imax,jmax,nr,nt)
      real*8 cfl,psi,vtrans
      real*8 vloc,fac,delip,deli,phi,eps,rp1h,rp1h_cfl
      real*8 recip_delip, recip_deli
      real*8 qippp,qipp,qip,qi,qim,qimm,qimmm
      real*8 mskipp,mskip,mski,mskim,mskimm,mskimmm
      real*8 d2,d2p1,d2m1,a,b,c,d
      real*8 dp1h,dm1h, phimd,philc,phimin,phimax
      real*8 delm,delp,delmm,delpp,delmmm,delppp
      real*8 del2mm,del2m,del2,del2p,del2pp
      real*8 del3mm,del3m,del3p,del3pp
      real*8 del4m,del4,del4p
      real*8 del5m,del5p
      real*8 del6
      parameter(eps = 1.d-20)

      do j=1,jmax
       do i=1,imax
          recipy(i,j)=1/deltay(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,jmax
          
          jp1=min(j+1,jmax)
          jp2=min(j+2,jmax)
          jp3=min(j+3,jmax)
          jm1=max(j-1, 1  )
          jm2=max(j-2, 1  )
          jm3=max(j-3, 1  )
          jm4=max(j-4, 1  )
          
         do i=1,imax           

          cfl = abs( vvel(i,j,k,t)*deltat*recipy(i,j) )
          vtrans=vvel(i,j,k,t)*ra(i,j,k)

			if (vtrans.gt.0.d0) then
			 qippp = tracer(i,j+2,k,t)
			 qipp  = tracer(i,j+1,k,t)
			 qip   = tracer(i,j  ,k,t)
			 qi    = tracer(i,j-1,k,t)
			 qim   = tracer(i,j-2,k,t)
			 qimm  = tracer(i,j-3,k,t)
			 qimmm = tracer(i,j-4,k,t)

			 mskipp  = maskv(i,j+2,k)
			 mskip   = maskv(i,j+1,k)
			 mski    = maskv(i,j  ,k)
			 mskim   = maskv(i,j-1,k)
			 mskimm  = maskv(i,j-2,k)
			 mskimmm = maskv(i,j-3,k)
			elseif (vtrans.lt.0.d0) then
			 qippp = tracer(i,j-3,k,t)
			 qipp  = tracer(i,j-2,k,t)
			 qip   = tracer(i,j-1,k,t)
			 qi    = tracer(i,j  ,k,t)
			 qim   = tracer(i,j+1,k,t)
			 qimm  = tracer(i,j+2,k,t)
			 qimmm = tracer(i,j+3,k,t)

			 mskipp  = maskv(i,j-2,k)
			 mskip   = maskv(i,j-1,k)
			 mski    = maskv(i,j  ,k)
			 mskim   = maskv(i,j+1,k)
			 mskimm  = maskv(i,j+2,k)
			 mskimmm = maskv(i,j+3,k)
			else
			 qippp = 0.d0
			 qipp  = 0.d0
			 qip   = 0.d0
			 qi    = 0.d0
			 qim   = 0.d0
			 qimm  = 0.d0
			 qimmm = 0.d0

			 mskipp  = 0.d0
			 mskip   = 0.d0
			 mski    = 0.d0
			 mskim   = 0.d0
			 mskimm  = 0.d0
			 mskimmm = 0.d0
			endif

        if (vtrans.ne.0.d0) then
C        2nd order correction [i i-1]
         fac = 1.d0
         delp = (qip-qi)*mski
         phi = fac * delp
C        3rd order correction [i i-1 i-2]
         fac = fac * ( cfl + 1.d0 )/3.d0
         delm = (qi-qim)*mskim
         del2 = delp - delm
         phi = phi - fac * del2
C        4th order correction [i+1 i i-1 i-2]
         fac = fac * ( cfl - 2.d0 )/4.d0
         delpp = (qipp-qip)*mskip*mski
         del2p = delpp - delp
         del3p = del2p - del2
         phi = phi + fac * del3p
C        5th order correction [i+1 i i-1 i-2 i-3]
         fac = fac * ( cfl - 3.d0 )/5.d0
         delmm = (qim-qimm)*mskimm*mskim
         del2m = delm - delmm
         del3m = del2 - del2m
         del4 = del3p - del3m
         phi = phi + fac * del4
C        6th order correction [i+2 i+1 i i-1 i-2 i-3]
         fac = fac * ( cfl + 2.d0 )/6.d0
         delppp = (qippp-qipp)*mskipp*mskip*mski
         del2pp = delpp - delp
         del3pp = del2pp - del2p
         del4p = del3pp - del3p
         del5p = del4p - del4
         phi = phi + fac * del5p
C        7th order correction [i+2 i+1 i i-1 i-2 i-3 i-4]
         fac = fac * ( cfl + 2.d0 )/7.d0
         delmmm = (qimm-qimmm)*mskimmm*mskimm*mskim
         del2mm = delmm - delmmm
         del3mm = del2m - del2mm
         del4m = del3m - del3mm
         del5m = del4 - del4m
         del6 = del5p - del5m
         phi = phi - fac * del6

         delip = ( qip - qi ) * mski
C        phi = sign(1.d0,phi)*sign(1.d0,delip)
C    &        *abs(phi+eps)/abs(delip+eps)
C     simplify and avoid division by zero
         recip_delip = sign(1.d0,delip)/max(abs(delip),eps)
         phi = phi*recip_delip

         deli = ( qi - qim ) * mskim
C        rp1h =sign(1.d0,deli)*sign(1.d0,delip)
C    &        *abs(deli+eps)/abs(delip+eps)
C     simplify and avoid division by zero
         recip_deli = sign(1.d0,deli)/max(abs(deli),eps)
         rp1h = deli*recip_delip
         rp1h_cfl = rp1h/(cfl+eps)

C        tvd limiter
C        phi = max(0.d0, min( 2./(1-cfl), phi, 2.*rp1h_cfl ) )

C        mp limiter
         d2   = del2 !( ( qip + qim ) - 2.*qi  ) * mski * mskim
         d2p1 = del2p !( ( qipp + qi ) - 2.*qip ) * mskip * mski
         d2m1 = del2m !( ( qi + qimm ) - 2.*qim ) * mskim * mskimm
         a = 4.d0*d2 - d2p1
         b = 4.d0*d2p1 - d2
         c = d2
         d = d2p1
         dp1h = max(min(a,b,c,d),0.d0)+min(max(a,b,c,d),0.d0)
         a = 4.d0*d2m1 - d2
         b = 4.d0*d2 - d2m1
         c = d2m1
         d = d2
         dm1h = max(min(a,b,c,d),0.d0)+min(max(a,b,c,d),0.d0)
C        qmd = 0.5*( ( qi + qip ) - dp1h )
C        qmd = 0.5d0*( ( 2.d0*qi + delip ) - dp1h )
C        qul = qi + (1.d0-cfl)/(cfl+eps)*deli
C        qlc = qi + 0.5d0*( 1.d0+dm1h/(deli+eps) )*(qul-qi)
C        phimd = 2.d0/(1.d0-cfl)*(qmd-qi+eps)/(delip+eps)
C        philc = 2.d0*rp1h_cfl*(qlc-qi+eps)/(qul-qi+eps)
C     simplify and avoid division by zero
         phimd = 1.d0/(1.d0-cfl)*(delip-dp1h)*recip_delip
         philc = rp1h_cfl*( 1.d0+dm1h*recip_deli )
C  
         phimin = max( min(0.d0,phimd),
     &                 min(0.d0,2.d0*rp1h_cfl,philc) )
         phimax = min( max(2.d0/(1.d0-cfl),phimd),
     &                 max(0.d0,2.d0*rp1h_cfl,philc) )
         phi = max(phimin,min(phi,phimax))

         psi = phi * 0.5d0 * (1.d0 - cfl)
         vt(i,j,k,t) = vtrans*( qi + psi*delip )
        else
         vt(i,j,k,t) = 0.d0
        endif
       enddo
       enddo
       enddo
      enddo

      return
      end

      subroutine advy_uw_3rd( imax,jmax,nr,nt,ra, 
     &              maskc,maskv,deltat,deltay,vvel,tracer,vt )

C calculates the third-order upwind biased area integrated meridional flux due to 
C advection of a tracer
      implicit none

      mwIndex i,j,k,t,jm1,jm2,jp1
      mwIndex imax,jmax,nr,nt
      real*8 deltat
      real*8 deltay(imax,jmax)
      real*8 ra(imax,jmax,nr)
      real*8 maskv(imax,jmax,nr)
      real*8 maskc(imax,jmax,nr)
      real*8 vvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 vt(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,rjjm,rjjp
      real*8 vtrans,vcfl,onesixth,half
      parameter(half=1.d0/2.d0)
      parameter(onesixth=1.d0/6.d0)
      
      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          vt(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax           
          jp1=min(j+1,jmax)
          jm1=max(j-1, 1  )
          jm2=max(j-2, 1  )
          
          rjp=(tracer(i,jp1,k,t)-tracer(i, j ,k,t))
     &     *maskc(i,jp1,k)*maskc(i, j ,k)
          rj =(tracer(i, j ,k,t)-tracer(i,jm1,k,t))
     &     *maskc(i, j ,k)*maskc(i,jm1,k)
          rjm=(tracer(i,jm1,k,t)-tracer(i,jm2,k,t))
     &     *maskc(i,jm1,k)*maskc(i,jm2,k)
          rjjp=rjp-rj
          rjjm=rj-rjm

          vtrans=vvel(i,j,k,t)*ra(i,j,k)
          
          vt(i,j,k,t) = maskc(i, j ,k)*
C     &                  maskc(i,jm1,k)*
C     &                  maskv(i, j ,k)*maskv(i,jm1,k)*     
C          vt(i,j,k,t) = 
     &   vtrans*(
     &     tracer(i,j,k,t)+tracer(i,jm1,k,t)-onesixth*( rjjp+rjjm )
     &               )*half
     &  +abs( vtrans )*half*onesixth*( rjjp-rjjm )
         enddo
        enddo
       enddo
      enddo  

      return
      end

      subroutine advy_uwcntr_2nd( imax,jmax,nr,nt,ra,maskc,
     &           maskv,vvel,tracer,vt )

C calculates the center second upwind biased area integrated meridional flux due to 
C advection of a tracer
      implicit none

      mwIndex i,j,k,t,jm1
      mwIndex imax,jmax,nr,nt
      real*8 maskc(imax,jmax,nr)
      real*8 maskv(imax,jmax,nr)
      real*8 ra(imax,jmax,nr)
      real*8 vvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 vt(imax,jmax,nr,nt)
      
         do t=1,nt
          do k=1,nr
           do j=1,jmax
            do i=1,imax
             jm1=max(j-1, 1  )
 
             if (vvel(i,j,k,t) .ge. 0.0) then
C               vt(i,j,k,t) = maskv(i,j,k) *maskv(i,jm1,k )*
C     &                       maskc(i,j,k) *maskc(i,jm1,k )*      
               vt(i,j,k,t) =     
     &                       vvel(i,j,k,t)*tracer(i,jm1,k,t)*ra(i,j,k)
             else
C               vt(i,j,k,t) = maskv(i,j,k )*maskc(i,j,k)*
               vt(i,j,k,t) = 
     &                       vvel(i,j,k,t)*tracer(i, j ,k,t)*ra(i,j,k)
             endif
            
            enddo
           enddo  
          enddo
         enddo

      return
      end

      subroutine advz_cntr_2nd( imax,jmax,nr,nt,ra,maskc,
     &           wvel,tracer,wt )

C calculates the centered second-order area integrated vertical flux due to advection of a tracer
      implicit none

      mwIndex i,j,k,t,km1,kp1
      mwIndex imax,jmax,nr,nt
      real*8 ra(imax,jmax)
      real*8 maskc(imax,jmax,nr)
      real*8 wvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 wt(imax,jmax,nr,nt)
      real*8 half
      parameter(half=1.d0/2.d0)
               
      do t=1,nt    
       do k=1,nr
        km1=max(1 ,k-1)
        kp1=min(nr,k+1)
        
        if ( k.eq.1 .or. k.gt.nr) then
          do j=1,jmax
           do i=1,imax
           wt(i,j,k,t) = 0.0000
          enddo
         enddo
        else
          do j=1,jmax
           do i=1,imax         
C             wt(i,j,k,t) = maskc(i,j,km1)*maskc(i,j,k)*
             wt(i,j,k,t) = 
     &       wvel(i,j,k,t)*ra(i,j)
     &       *(tracer(i,j,k,t)+tracer(i,j,km1,t))*half         
          enddo
         enddo
        endif
       enddo
      enddo
      
      return
      end

      subroutine advz_cntr_4th(imax,jmax,nr,nt,
     &           ra,maskc,wvel,tracer,wt )

C calculates the centered fourth-order area integrated vertical flux due to advection of a tracer
      implicit none

      mwIndex i,j,k,t,km1,km2,kp1
      mwIndex imax,jmax,nr,nt
      real*8 maskc(imax,jmax,nr)
      real*8 ra(imax,jmax)
      real*8 wvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 wt(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,rjjm,rjjp
      real*8 rtrans,half,onesixth,maskpm,maskbound
      parameter(half=1.d0/2.d0)
      parameter(onesixth=1.d0/6.d0)

      do t=1,nt    
       do k=1,nr
        km2=max(1,k-2)
        km1=max(1,k-1)
        kp1=min(nr,k+1)
        maskpm = 1.d0
        if (k.le.2 .or. k.ge.nr) maskpm = 0.d0
   
        if ( k.eq.1 .or. k.gt.nr) then
         do j=1,jmax
          do i=1,imax
           wt(i,j,k,t) = 0.d0
          enddo
         enddo
        else
         do j=1,jmax
          do i=1,imax
           maskbound = maskpm*maskc(i,j,km2)*maskc(i,j,kp1)
           rjp=(tracer(i,j,kp1,t)-tracer(i,j, k ,t))
C     &          *maskc(i,j,kp1)
           rj =(tracer(i,j, k ,t)-tracer(i,j,km1,t))
           rjm=(tracer(i,j,km1,t)-tracer(i,j,km2,t))
C     &          *maskc(i,j,km1)
           rjjp=(rjp-rj)
           rjjm=(rj-rjm)
           
           rtrans=wvel(i,j,k,t)*ra(i,j)

C           wt(i,j,k,t) = maskc(i,j,km1)*maskc(i,j,k)
C     &       *(
           wt(i,j,k,t) = 

     &       rtrans*(
     &         (tracer(i,j,k,t)+tracer(i,j,km1,t))*half
     &         -onesixth*(rjjm+rjjp)*half 
     &              )
     &         +abs(rtrans)*
     &         onesixth*(rjjm-rjjp)*half*(1.d0 - maskbound)
C     &        )                      
          enddo
         enddo
        endif
       enddo
      enddo
      
      return
      end

      subroutine advz_cntrfl_2nd( imax,jmax,nr,nt,ra, 
     &              maskc,deltat,deltaz,wvel,tracer,wt )

C calculates the centered second-order area integrated meridional flux due to 
C advection of a tracer with a flux limiter
      implicit none

      mwIndex i,j,k,t,km1,km2,kp1
      mwIndex imax,jmax,nr,nt
      real*8 deltat
      real*8 deltaz(nr)
      real*8 recipz(nr)
      real*8 ra(imax,jmax)
      real*8 maskc(imax,jmax,nr)
      real*8 maskup(imax,jmax,nr)
      real*8 wvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 wt(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,wcfl,cr,lcr
      real*8 rtrans,thetamax,half
      parameter(half=1.d0/2.d0)
      parameter(thetamax = 1.d+20 )

      km2=max(1,k-2)
      km1=max(1,k-1)
      kp1=min(nr,k+1)

      do k=1,nr
       recipz(k)=1/deltaz(k)
      enddo
      
      do t=1,nt    
       do k=1,nr
        km2=max(1,k-2)
        km1=max(1,k-1)
        kp1=min(nr,k+1)

         do j=1,jmax
          do i=1,imax
           rjp=(tracer(i,j,k,t)-tracer(i,j,kp1,t))
C     &         *maskc(i,j,kp1)
           rj =(tracer(i,j,km1,t)-tracer(i,j,k,t))
C     &         *maskc(i,j,k)*maskc(i,j,km1)
           rjm=(tracer(i,j,km2,t)-tracer(i,j,km1,t))
C     &         *maskc(i,j,km1)

           rtrans= wvel(i,j,k,t)*ra(i,j)
		   wcfl = abs(wvel(i,j,k,t)*deltat*recipz(k))

         if (rj.ne.0.) then
          if (rtrans.lt.0.) then
            cr=rjm/rj
          else
            cr=rjp/rj
          endif
         else
          if (rtrans.lt.0.) then
            cr=rjm*thetamax
          else
            cr=rjp*thetamax
          endif
         endif

C the flux limiter can be -
C upwind:
C     limiter(cr)=0.
C lax-wendroff:
C     limiter(cr)=1.
C min-mod:
C     limiter(cr)=max(0.d0,min(1.d0,cr))
C suberbee:
         lcr=max(0.d0,max(min(1.d0,2.d0*cr),
     &                         min(2.d0,cr)))
              
C         wt(i,j,k,t) = maskc(i,j,km1)*maskc(i,j,k)*(
         wt(i,j,k,t) = (
     &     rtrans*
     &        (tracer(i,j,k,t)+tracer(i,j,km1,t))*half
     &    +abs(rtrans)*((1.d0-lcr)+wcfl*lcr)
     &                     *rj*half )
         enddo
        enddo
       enddo
      enddo  

      return
      end

      subroutine advz_dst_2nd(
     &           imax,jmax,nr,nt,advectionscheme,ra,
     &           maskc,deltat,deltaz,wvel,tracer,wt )

C  calculates the second-order direct space and time (dst-2/lax-wendroff, or simple 
C  first order upwind scheme) area integrated vertical flux due to advection of a tracer 
      implicit none

      mwIndex i,j,k,t,km1
      mwIndex imax,jmax,nr,nt
      real*8 advectionscheme, deltat
      real*8 wcfl, rlimit, wabs, rtrans
      real*8 deltaz(nr)
      real*8 recipz(nr)
      real*8 ra(imax,jmax)
      real*8 maskc(imax,jmax,nr)
      real*8 maskup(imax,jmax,nr)
      real*8 wvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 wt(imax,jmax,nr,nt)
      real*8 rksign
      parameter(rksign=-1.d0)
       
      if ( advectionscheme.eq.20 ) then
       rlimit = 1.d0
      else
       rlimit = 0.d0
      endif  

      do k=1,nr
       recipz(k)=1/deltaz(k)
      enddo
      
      do t=1,nt    
       do k=1,nr
        
        km1=max(1,k-1)
        
        if ( k.le.1 .or. k.gt.nr) then
         do j=1,jmax
          do i=1,imax
           wt(i,j,k,t) = 0.d0
          enddo
         enddo
        else
         do j=1,jmax
          do i=1,imax
          wcfl = abs( wvel(i,j,k,t)*deltat*recipz(k) )
          rtrans= wvel(i,j,k,t)*ra(i,j)

          wabs = abs(rtrans)*rksign
     &       *( 1.d0 - rlimit*(1.d0 - wcfl) )
          wt(i,j,k,t) = maskc(i,j,km1)*(
     &             ( rtrans+wabs )* 0.5d0 * tracer(i,j,km1,t)
     &           + ( rtrans-wabs )* 0.5d0 * tracer(i,j,k  ,t)
     &                                  )
          enddo
         enddo
        endif
       enddo
      enddo
      return
      end

      subroutine advz_dst_3rd( imax,jmax,nr,nt,ra, 
     &              maskc,deltat,deltaz,wvel,tracer,wt )

C  calculates the 3rd-order direct space and time (dst-3) area integrated vertical flux 
C  due to advection of a tracer
      implicit none

      mwIndex i,j,k,t,km1,km2,kp1
      mwIndex imax,jmax,nr,nt
      real*8 deltat
      real*8 deltaz(nr)
      real*8 recipz(nr)
      real*8 ra(imax,jmax)
      real*8 maskc(imax,jmax,nr)
      real*8 maskup(imax,jmax,nr)
      real*8 wvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 wt(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,d0,d1,onesixth
      real*8 rtrans,wcfl,half
      parameter(half=1.d0/2.d0)
      parameter(onesixth=1.d0/6.d0)

      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          wt(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do k=1,nr
       recipz(k)=1/deltaz(k)
      enddo
      
      do t=1,nt    
       do k=1,nr
        km2=max(1,k-2)
        km1=max(1,k-1)
        kp1=min(nr,k+1)

         do j=1,jmax
          do i=1,imax
           rjp=(tracer(i,j,k,t)-tracer(i,j,kp1,t))
     &         *maskc(i,j,kp1)*maskc(i,j,k)
           rj =(tracer(i,j,km1,t)-tracer(i,j,k,t))
     &         *maskc(i,j,k)*maskc(i,j,km1)
           rjm=(tracer(i,j,km2,t)-tracer(i,j,km1,t))
     &         *maskc(i,j,km1)*maskc(i,j,km2)

           rtrans= wvel(i,j,k,t)*ra(i,j)
           wcfl=abs( wvel(i,j,k,t)*deltat*recipz(k) )

           d0=(2.d0-wcfl)*(1.d0-wcfl)*onesixth
           d1=(1.d0-wcfl*wcfl)*onesixth

           wt(i,j,k,t)=maskc(i,j,k)*maskc(i,j,km1)           
C           wt(i,j,k,t)=(
     &    *(half*(rtrans+abs(rtrans))
     &       *( tracer(i,j, k ,t) + (d0*rj+d1*rjp) )
     &   +half*(rtrans-abs(rtrans))
     &       *( tracer(i,j,km1,t) - (d0*rj+d1*rjm) ))
         enddo
        enddo
       enddo
      enddo
      
      return
      end

      subroutine advz_dstfl_3rd( imax,jmax,nr,nt,ra, 
     &              maskc,deltat,deltaz,wvel,tracer,wt )

C  calculates the 3rd-order direct space and time (dst-3) area integrated zonal flux 
C  due to advection of a tracer with flux limiting.
      implicit none

      mwIndex i,j,k,t,km1,km2,kp1
      mwIndex imax,jmax,nr,nt
      real*8 deltat
      real*8 deltaz(nr)
      real*8 recipz(nr)
      real*8 ra(imax,jmax)
      real*8 maskc(imax,jmax,nr)
      real*8 maskup(imax,jmax,nr)
      real*8 wvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 wt(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,d0,d1,psip,psim
      real*8 thetamin,thetamax,thetap,thetam
      real*8 rtrans,wcfl,onesixth,half
      parameter(half=1.d0/2.d0)
      parameter(thetamax = 1.d+20 )
      parameter(thetamin = 1.d-20 )
      parameter(onesixth=1.d0/6.d0)

      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          wt(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do k=1,nr
       recipz(k)=1/deltaz(k)
      enddo
      
      do t=1,nt    
       do k=1,nr
        km2=max(1,k-2)
        km1=max(1,k-1)
        kp1=min(nr,k+1)

         do j=1,jmax
          do i=1,imax
           rjp=(tracer(i,j, k ,t)-tracer(i,j,kp1,t))
     &         *maskc(i,j,kp1)*maskc(i,j, k )
           rj =(tracer(i,j,km1,t)-tracer(i,j, k ,t))
     &         *maskc(i,j, k )*maskc(i,j,km1)
           rjm=(tracer(i,j,km2,t)-tracer(i,j,km1,t))
     &         *maskc(i,j,km1)*maskc(i,j,km2)

		   wcfl = abs( wvel(i,j,k,t)*deltat*recipz(k) )
           rtrans= wvel(i,j,k,t)*ra(i,j)
		   
		   d0=(2.d0 -wcfl)*(1.d0 -wcfl)*onesixth
		   d1=(1.d0 -wcfl*wcfl)*onesixth

C      prevent |thetap,m| to reach too big value
           if ( abs(rj)*thetamax .le. abs(rjm) ) then
             thetap=sign(thetamax,rjm*rj)
           else
             thetap=rjm/rj
           endif
           if ( abs(rj)*thetamax .le. abs(rjp) ) then
             thetam=sign(thetamax,rjp*rj)
           else
             thetam=rjp/rj
           endif

           psip=d0+d1*thetap
           psip=max(0.d0,min(min(1.d0,psip),
     &          thetap*(1.d0 -wcfl)/(wcfl+thetamin) ))
           psim=d0+d1*thetam
           psim=max(0.d0,min(min(1.d0,psim),
     &          thetam*(1.d0 -wcfl)/(wcfl+thetamin) ))

           wt(i,j,k,t)=maskc(i,j,k)*maskc(i,j,km1)      
C           wt(i,j,k,t)=(
     &   *( half*(rtrans+abs(rtrans))
     &      *( tracer(i,j, k ,t) + psim*rj )
     &  +half*(rtrans-abs(rtrans))
     &      *( tracer(i,j,km1,t) - psip*rj ))
         enddo
        enddo
       enddo
      enddo  

      return
      end

      subroutine advz_mpl_7th(imax,jmax,nr,nt,ra, 
     &              maskc,deltat,deltaz,wvel,tracer,wt )
     
C calculates the 7th order dst area integrated zonal flux due to 
C advection of a tracer with a monotone preserving limiter
      implicit none

      mwIndex imax,jmax,nr,nt,i,j,k,t
      mwIndex kp3,kp2,kp1,km1,km2,km3,km4
      real*8 deltat
      real*8 deltaz(nr)
      real*8 recipz(nr)
      real*8 ra(imax,jmax)
      real*8 maskc(imax,jmax,nr)
      real*8 maskup(imax,jmax,nr)
      real*8 wvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 wt(imax,jmax,nr,nt)
      real*8 cfl,psi,wtrans
      real*8 wloc,fac,delip,deli,phi,eps,rp1h,rp1h_cfl
      real*8 recip_delip, recip_deli
      real*8 qippp,qipp,qip,qi,qim,qimm,qimmm
      real*8 mskipp,mskip,mski,mskim,mskimm,mskimmm
      real*8 d2,d2p1,d2m1,a,b,c,d
      real*8 dp1h,dm1h, phimd,philc,phimin,phimax
      real*8 delm,delp,delmm,delpp,delmmm,delppp
      real*8 del2mm,del2m,del2,del2p,del2pp
      real*8 del3mm,del3m,del3p,del3pp
      real*8 del4m,del4,del4p
      real*8 del5m,del5p
      real*8 del6
      parameter(eps = 1.d-20)
      
      do k=1,nr
       recipz(k)=1/deltaz(k)
      enddo
      
      do t=1,nt    
       do k=1,nr
       
        km4=max(1,k-4)
        km3=max(1,k-3)
        km2=max(1,k-2)
        km1=max(1,k-1)
        kp1=min(nr,k+1)
        kp2=min(nr,k+2)
        kp3=min(nr,k+3)

         do j=1,jmax
          do i=1,imax
           wtrans= wvel(i,j,k,t)*ra(i,j)
           cfl=abs( wvel(i,j,k,t)*deltat*recipz(k) )

			if (wtrans.lt.0.d0) then
			 qippp = tracer(i,j,kp2,t)
			 qipp  = tracer(i,j,kp1,t)
			 qip   = tracer(i,j,k  ,t)
			 qi    = tracer(i,j,km1,t)
			 qim   = tracer(i,j,km2,t)
			 qimm  = tracer(i,j,km3,t)
			 qimmm = tracer(i,j,km4,t)

			 mskipp  = maskc(i,j,kp2) * float(kp2-kp1)
			 mskip   = maskc(i,j,kp1) * float(kp1-k)
			 mski    = maskc(i,j,k  ) * float(k-km1)
			 mskim   = maskc(i,j,km1) * float(km1-km2)
			 mskimm  = maskc(i,j,km2) * float(km2-km3)
			 mskimmm = maskc(i,j,km3) * float(km3-km4)
			elseif (wtrans.gt.0.d0) then
			 qippp = tracer(i,j,km3,t)
			 qipp  = tracer(i,j,km2,t)
			 qip   = tracer(i,j,km1,t)
			 qi    = tracer(i,j,k  ,t)
			 qim   = tracer(i,j,kp1,t)
			 qimm  = tracer(i,j,kp2,t)
			 qimmm = tracer(i,j,kp3,t)

			 mskipp  = maskc(i,j,km2) * float(km2-km3)
			 mskip   = maskc(i,j,km1) * float(km1-km2)
			 mski    = maskc(i,j,k  ) * float(k-km1)
			 mskim   = maskc(i,j,kp1) * float(kp1-k)
			 mskimm  = maskc(i,j,kp2) * float(kp2-kp1)
			 mskimmm = maskc(i,j,kp3) * float(kp3-kp2)
			else
			 qippp = 0.d0
			 qipp  = 0.d0
			 qip   = 0.d0
			 qi    = 0.d0
			 qim   = 0.d0
			 qimm  = 0.d0
			 qimmm = 0.d0

			 mskipp  = 0.d0
			 mskip   = 0.d0
			 mski    = 0.d0
			 mskim   = 0.d0
			 mskimm  = 0.d0
			 mskimmm = 0.d0
			endif

        if (wtrans.ne.0.d0) then
C        2nd order correction [i i-1]
         fac = 1.d0
         delp = (qip-qi)*mski
         phi = fac * delp
C        3rd order correction [i i-1 i-2]
         fac = fac * ( cfl + 1.d0 )/3.d0
         delm = (qi-qim)*mskim
         del2 = delp - delm
         phi = phi - fac * del2
C        4th order correction [i+1 i i-1 i-2]
         fac = fac * ( cfl - 2.d0 )/4.d0
         delpp = (qipp-qip)*mskip*mski
         del2p = delpp - delp
         del3p = del2p - del2
         phi = phi + fac * del3p
C        5th order correction [i+1 i i-1 i-2 i-3]
         fac = fac * ( cfl - 3.d0 )/5.d0
         delmm = (qim-qimm)*mskimm*mskim
         del2m = delm - delmm
         del3m = del2 - del2m
         del4 = del3p - del3m
         phi = phi + fac * del4
C        6th order correction [i+2 i+1 i i-1 i-2 i-3]
         fac = fac * ( cfl + 2.d0 )/6.d0
         delppp = (qippp-qipp)*mskipp*mskip*mski
         del2pp = delpp - delp
         del3pp = del2pp - del2p
         del4p = del3pp - del3p
         del5p = del4p - del4
         phi = phi + fac * del5p
C        7th order correction [i+2 i+1 i i-1 i-2 i-3 i-4]
         fac = fac * ( cfl + 2.d0 )/7.d0
         delmmm = (qimm-qimmm)*mskimmm*mskimm*mskim
         del2mm = delmm - delmmm
         del3mm = del2m - del2mm
         del4m = del3m - del3mm
         del5m = del4 - del4m
         del6 = del5p - del5m
         phi = phi - fac * del6

         delip = ( qip - qi ) * mski
C        phi = sign(1.d0,phi)*sign(1.d0,delip)
C    &        *abs(phi+eps)/abs(delip+eps)
C     simplify and avoid division by zero
         recip_delip = sign(1.d0,delip)/max(abs(delip),eps)
         phi = phi*recip_delip

         deli = ( qi - qim ) * mskim
C        rp1h =sign(1.d0,deli)*sign(1.d0,delip)
C    &        *abs(deli+eps)/abs(delip+eps)
C     simplify and avoid division by zero
         recip_deli = sign(1.d0,deli)/max(abs(deli),eps)
         rp1h = deli*recip_delip
         rp1h_cfl = rp1h/(cfl+eps)

C        tvd limiter
C        phi = max(0.d0, min( 2./(1-cfl), phi, 2.*rp1h_cfl ) )

C        mp limiter
         d2   = del2 !( ( qip + qim ) - 2.*qi  ) * mski * mskim
         d2p1 = del2p !( ( qipp + qi ) - 2.*qip ) * mskip * mski
         d2m1 = del2m !( ( qi + qimm ) - 2.*qim ) * mskim * mskimm
         a = 4.d0*d2 - d2p1
         b = 4.d0*d2p1 - d2
         c = d2
         d = d2p1
         dp1h = max(min(a,b,c,d),0.d0)+min(max(a,b,c,d),0.d0)
         a = 4.d0*d2m1 - d2
         b = 4.d0*d2 - d2m1
         c = d2m1
         d = d2
         dm1h = max(min(a,b,c,d),0.d0)+min(max(a,b,c,d),0.d0)
C        qmd = 0.5*( ( qi + qip ) - dp1h )
C        qmd = 0.5d0*( ( 2.d0*qi + delip ) - dp1h )
C        qul = qi + (1.d0-cfl)/(cfl+eps)*deli
C        qlc = qi + 0.5d0*( 1.d0+dm1h/(deli+eps) )*(qul-qi)
C        phimd = 2.d0/(1.d0-cfl)*(qmd-qi+eps)/(delip+eps)
C        philc = 2.d0*rp1h_cfl*(qlc-qi+eps)/(qul-qi+eps)
C     simplify and avoid division by zero
         phimd = 1.d0/(1.d0-cfl)*(delip-dp1h)*recip_delip
         philc = rp1h_cfl*( 1.d0+dm1h*recip_deli )
C  
         phimin = max(min(0.d0,phimd),
     &        min(0.d0,2.d0*rp1h_cfl,philc))
         phimax = min(max(2.d0/(1.d0-cfl),phimd),
     &        max(0.d0,2.d0*rp1h_cfl,philc))
         phi = max(phimin,min(phi,phimax))

         psi = phi * 0.5d0 * (1.d0 - cfl)
         wt(i,j,k,t) = wtrans*( qi + psi*delip )
        else
         wt(i,j,k,t) = 0.d0
        endif

       enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine advz_uw_3rd( imax,jmax,nr,nt,ra, 
     &              maskc,deltat,deltaz,wvel,tracer,wt )
     
C calculates the third-order upwind biased area integrated vertical flux due to 
C advection of a tracer
      implicit none

      mwIndex i,j,k,t,km1,km2,kp1
      mwIndex imax,jmax,nr,nt
      real*8 deltat
      real*8 deltaz(nr)
      real*8 ra(imax,jmax)
      real*8 maskc(imax,jmax,nr)
      real*8 maskup(imax,jmax,nr)
      real*8 wvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 wt(imax,jmax,nr,nt)
      real*8 rjm,rj,rjp,rjjp,rjjm
      real*8 rtrans,onesixth,half
      parameter(half=1.d0/2.d0)            
      parameter(onesixth=1.d0/6.d0)

      do t=1,nt    
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          wt(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo             
      
      do t=1,nt    
       do k=1,nr
        km2=max(1,k-2)
        km1=max(1,k-1)
        kp1=min(nr,k+1)

         do j=1,jmax
          do i=1,imax
           rjp=(tracer(i,j,kp1,t)-tracer(i,j, k ,t))
     &         *maskc(i,j,kp1)*maskc(i,j,k)
           rj =(tracer(i,j, k ,t)-tracer(i,j,km1,t))
     &         *maskc(i,j,k)*maskc(i,j,km1)
           rjm=(tracer(i,j,km1,t)-tracer(i,j,km2,t))
     &         *maskc(i,j,km1)*maskc(i,j,km2)
           
           rjjp = rjp-rj
           rjjm = rj-rjm

           rtrans= wvel(i,j,k,t)*ra(i,j)
         
           wt(i,j,k,t) = maskc(i,j,k)
C     &                  *maskc(i,j,km1)
C            wt(i,j,k,t) =
     &       *(rtrans*( (tracer(i,j,k,t)+tracer(i,j,km1,t))*half
     &                  -onesixth*(rjjm+rjjp)*half  )
     &      +abs(rtrans)*
     &                   onesixth*(rjjm-rjjp)*half
     &                                  )
         enddo
        enddo
       enddo
      enddo  

      return
      end

      subroutine advz_uwcntr_2nd( imax,jmax,nr,nt,ra,maskc,
     &           wvel,tracer,wt )

C calculates the center second upwind biased area integrated vertical flux due to 
C advection of a tracer
      implicit none

      mwIndex i,j,k,t,km1,kp1
      mwIndex imax,jmax,nr,nt

      real*8 ra(imax,jmax)
      real*8 maskc(imax,jmax,nr)
      real*8 wvel(imax,jmax,nr,nt)
      real*8 tracer(imax,jmax,nr,nt)
      real*8 wt(imax,jmax,nr,nt)
               
      do t=1,nt    
       do k=1,nr
        km1=max(1,k-1)
        kp1=min(k+1,nr)
         do j=1,jmax
          do i=1,imax
            if ( k .eq. 1 ) then
              wt(i,j,1,t) = 0.d0
            endif
            
            if ( wvel(i,j,kp1,t) .ge. 0.d0 ) then
C              wt(i,j,kp1,t) = maskc(i,j,km1)*maskc(i,j,kp1)*
              wt(i,j,kp1,t) = 
     &                        wvel(i,j,kp1,t)*tracer(i,j,kp1,t)*ra(i,j)
            else
C              wt(i,j,kp1,t) = maskc(i,j,km1)*maskc(i,j,kp1)*
              wt(i,j,kp1,t) = 
     &                        wvel(i,j,kp1,t)*tracer(i,j, k ,t)*ra(i,j)
            endif
          enddo
         enddo
       enddo
      enddo
      
      return
      end
