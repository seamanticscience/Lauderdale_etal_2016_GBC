#include "fintrf.h"
C-----------------------------------------------------------------------------------------
C     implicit_diffusion_flux calls implicit routine to calculate diffusive transport
C     [wflux,itertrac]=implicit_diffusion_flux(timestep,niters,rac,
C                               dzf,dzc,mask,diffkz,fld);
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
      real*8 deltat, niters

C     Allocatable arrays for fortran data inputs to gateway routine
      real*8, allocatable, dimension(:) :: dzf
      real*8, allocatable, dimension(:) :: dzc
      real*8, allocatable, dimension(:) :: racin
      real*8, allocatable, dimension(:) :: mskcin
      real*8, allocatable, dimension(:) :: diffkzin
      real*8, allocatable, dimension(:) :: fin

C     Allocatable arrays for fortran data inputs to computational routine
      real*8, allocatable, dimension(:,:)     :: rac
      real*8, allocatable, dimension(:,:,:)   :: mskc
      real*8, allocatable, dimension(:,:,:,:) :: diffkz
      real*8, allocatable, dimension(:,:,:,:) :: fld

C     Allocatable arrays for fortran data outputs  
      real*8, allocatable, dimension(:,:,:,:) :: wflux

C---------------VALIDATE------------INPUTS--------------------------------
C     Check for proper number of arguments. 
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Calculating Implicit Diffusive Fluxes')
      call mexWarnMsgTxt('Check for proper number of arguments.')

      write(line,*)'Number of inputs is ',nrhs
      k=mexPrintf(line//achar(13)) 
      write(line,*)'Number of outputs is ',nlhs
      k=mexPrintf(line//achar(13)) 
#endif
      if(nrhs .ne. 8) then
         call mexErrMsgTxt('ERROR: 15 inputs required')
      elseif(nlhs .gt. 1) then
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
      call mexWarnMsgTxt('allocate racin array dimensions'   )       
#endif
      allocate (racin  (mxGetNumberOfElements(prhs(3))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate dzf array dimensions'   )       
#endif
      allocate (dzf  (mxGetNumberOfElements(prhs(4))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate dzc array dimensions'   )       
#endif
      allocate (dzc   (mxGetNumberOfElements(prhs(5))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mskcin array dimensions')  
#endif
      allocate (mskcin(mxGetNumberOfElements(prhs(6))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate diffkzin array dimensions'  )
#endif
      allocate (diffkzin(mxGetNumberOfElements(prhs(7))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate fin array dimensions')   
#endif
      allocate (fin(mxGetNumberOfElements(prhs(8))))

C Allocate the input to the computational function
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate rac array dimensions')    
#endif
      allocate (rac(dims(1),dims(2)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mskc array dimensions')    
#endif
      allocate (mskc(dims(1),dims(2),dims(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate diffkz array dimensions')    
#endif
      allocate (diffkz(dims(1),dims(2),dims(3),dims(4)))          
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate fld array dimensions')    
#endif
      allocate (fld(dims(1),dims(2),dims(3),dims(4))) 

C Allocate the output from the computational function      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate wflux array dimensions')    
#endif
      allocate (wflux(dims(1),dims(2),dims(3),dims(4)))          

C     Fill Fortran arrays from the input argument pointers.
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill arrays from the input.') 
      call mexWarnMsgTxt('Get deltat from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(1)), deltat, 
     & mxGetNumberOfElements(prhs(1)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Get niters from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), niters, 
     & mxGetNumberOfElements(prhs(2)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill racin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), racin, 
     & mxGetNumberOfElements(prhs(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill dzf array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(4)), dzf, 
     & mxGetNumberOfElements(prhs(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill dzc array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(5)), dzc, 
     & mxGetNumberOfElements(prhs(5)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill mskcin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(6)), mskcin, 
     & mxGetNumberOfElements(prhs(6)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill diffkzin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(7)), diffkzin, 
     & mxGetNumberOfElements(prhs(7)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill fin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(8)), fin, 
     & mxGetNumberOfElements(prhs(8)))
     
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Create arrays for the output arguments.')  
#endif
      classid = mxClassIDFromClassName('double')
      complexflag = 0  
      plhs(1) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)
      
C Check arguments for NANs in the input
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check arguments for NANs')    
      call mexWarnMsgTxt('Check deltat')
#endif
      if( mxIsNaN(deltat) .eq. 1 ) then
C      deltat is pretty important, so error if NaN
	   call mexErrMsgTxt('ERROR: Input deltat contains a NaN.')
      endif
            
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check niters')
#endif
      if( mxIsNaN(niters) .eq. 1 ) then
C        call mexErrMsgTxt('ERROR: Input niters contains a NaN.')
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Input niters contains a NaN.')
#endif
         niters=1 .d 0
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
       do iter=1,mxGetNumberOfElements(prhs(3))
        if( mxIsNaN(racin(iter)) .eq. 1 ) then
         racin(iter)=0 .d 0
        endif
       enddo
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check dzf')    
#endif
      if (mxIsNaN(dzf) .eq. 1) then
C      dzf is pretty important, so error if NaN
       call mexErrMsgTxt('ERROR: Input dzf contains a NaN.')
      endif
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check dzc')    
#endif
      if (mxIsNaN(dzc) .eq. 1) then
C      dzc is pretty important, so error if NaN
       call mexErrMsgTxt('ERROR: Input dzc contains a NaN.')
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
      call mexWarnMsgTxt('Check diffkzin')    
#endif
      if (mxIsNaN(diffkzin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input diffkzin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(7)) 
        if( mxIsNaN(diffkzin(iter)) .eq. 1 ) then
         diffkzin(iter)=0 .d 0
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
       do iter=1,mxGetNumberOfElements(prhs(8))
        if( mxIsNaN(fin(iter)) .eq. 1 ) then
         fin(iter)=0 .d 0
        endif
       enddo
      endif
      
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Reshape input arrays') 
#endif   
      rac   = reshape(racin   ,(/dims(1),dims(2)/))
      mskc  = reshape(mskcin  ,(/dims(1),dims(2),dims(3)/))
      diffkz= reshape(diffkzin,(/dims(1),dims(2),dims(3),dims(4)/))
      fld   = reshape(fin     ,(/dims(1),dims(2),dims(3),dims(4)/))

C Fill the output arrays with zeros
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Fill output arrays with zeros') 
#endif   
      do t=1,dims(4)
       do k=1,dims(3)
        do j=1,dims(2)
         do i=1,dims(1)
          wflux(i,j,k,t)=0.D0
         enddo
        enddo  
       enddo
      enddo
      
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Call computational subroutine')  
#endif        
C     Call the computational subroutines.
         call implicit_mixing(dims(1),dims(2),dims(3),dims(4),
     &           deltat, niters, mskc, rac, dzf, dzc, diffkz,fld,
     &           wflux)

C     Load the data into the pointers, which are the output to MATLAB.
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Load the data into the output pointers')    
#endif
      size=dims(1)*dims(2)*dims(3)*dims(4)
      call mxCopyReal8ToPtr(wflux, mxGetPr(plhs(1)),size)
     
      return
      end
      
CC------------------------COMPUTATIONAL SUBROUTINES------------------------
      subroutine implicit_mixing(imax, jmax, nr, nt,
     &     deltat,niters,maskc,ra,dzf,dzc,kappaz,
     &     fld,flux)

c     solve implicitly vertical advection and diffusion
      implicit none

      mwIndex imax,jmax,nr,nt
      mwIndex i,j,k,t,iter
      real*8 niters
      real*8 deltat
      real*8 recip_drf(nr)
      real*8 recip_drc(nr)
      real*8 dzf      (nr)
      real*8 dzc      (nr)
      real*8 ra       (imax,jmax)
      real*8 maskc    (imax,jmax,nr)
      real*8 tmpfld   (imax,jmax,nr)
      real*8 iterfld  (imax,jmax,nr)
      real*8 fld      (imax,jmax,nr,nt)
      real*8 kappaz   (imax,jmax,nr,nt)
      real*8 flux     (imax,jmax,nr,nt)
      real*8 b5d      (imax,jmax,nr)
      real*8 c5d      (imax,jmax,nr)
      real*8 d5d      (imax,jmax,nr)    
      real*8 rksign
      parameter(rksign=-1.d0)

      do t=1,nt
C     initialise
       do k=1,nr        
        recip_drf(k)=1.d0/dzf(k)
        recip_drc(k)=1.d0/dzc(k)
        do j=1,jmax
         do i=1,imax
          b5d(i,j,k) = 0.d0
          c5d(i,j,k) = 1.d0
          d5d(i,j,k) = 0.d0
         enddo
        enddo
       enddo

C     set the tri-diagonal matrix to solve the implicit diffusion problem
C      1rst lower diagonal :
        do k=2,nr
         do j=1,jmax
          do i=1,imax
            b5d(i,j,k) =-deltat
     &                  *maskc(i,j,k-1)*maskc (i,j,k)  
     &                  *recip_drf(k)
     &                  *kappaz(i,j,k,t)*recip_drc(k)
          enddo
         enddo
        enddo
C       1rst upper diagonal :
        do k=1,nr-1
         do j=1,jmax
          do i=1,imax
            d5d(i,j,k) = -deltat
     &                 *maskc(i,j,k+1)*maskc (i,j,k)    
     &                 *recip_drf(k)
     &                 *kappaz(i,j,k+1,t)*recip_drc(k+1)
          enddo
         enddo
        enddo
C       main diagonal :
        do k=1,nr
         do j=1,jmax
          do i=1,imax
            c5d(i,j,k) = 1.0000 - b5d(i,j,k) - d5d(i,j,k)
          enddo
         enddo
        enddo

C       get local copy of fld for the tridiag solver
        do k=1,nr
         do j=1,jmax
          do i=1,imax
           tmpfld(i,j,k)=fld(i,j,k,t)
          enddo
         enddo
        enddo
      
C iterate a few times
      do iter=1,int(niters)
C       solve tri-diagonal system :
         call solve_matrix( imax, jmax, nr, b5d, c5d, d5d, tmpfld)

        do k=1,nr
         do j=1,jmax
          do i=1,imax
           iterfld(i,j,k)=tmpfld(i,j,k)
          enddo
         enddo
        enddo
C end of the iteration loop
       enddo
       
C      calculate diffusive flux
       do k=nr,1,-1
        if ( k.ge.2 ) then
         do j=1,jmax
          do i=1,imax
           flux(i,j,k,t) =
     &         -ra(i,j)*kappaz(i,j,k,t)*recip_drc(k)*rksign
     &         * (tmpfld(i,j,k) - tmpfld(i,j,k-1))
     &         * maskc(i,j,k)
     &         * maskc(i,j,k-1)
          enddo
         enddo
        else
         do j=1,jmax
          do i=1,imax
           flux(i,j,k,t) = 0.d0
          enddo
         enddo
        endif
       enddo

C     end of time loop
      enddo
      
      return
      end
      
      subroutine solve_matrix(imax, jmax, nr, a3d, b3d, c3d, y3d)

C      used to solve implicitly vertical advection & diffusion
C     a3d :: matrix lower diagnonal
C     b3d :: matrix main  diagnonal
C     c3d :: matrix upper diagnonal
C     y3d :: input = y vector ; output = x = solution of a*x=y
      implicit none

      mwIndex imax,jmax,nr,i,j,k
      real*8 a3d(imax,jmax,nr)
      real*8 b3d(imax,jmax,nr)
      real*8 c3d(imax,jmax,nr)
      real*8 y3d(imax,jmax,nr)
      real*8 tmpvar, recvar 
      real*8 c3d_prime (imax,jmax,nr)
      real*8 y3d_prime (imax,jmax,nr)      
      real*8 y3d_update(imax,jmax,nr)

      do k=1,nr
       do j=1,jmax
        do i=1,imax
         c3d_prime(i,j,k) = 0.d0
         y3d_prime(i,j,k) = 0.d0
         y3d_update(i,j,k) = 0.d0
        enddo
       enddo
      enddo

C     forward sweep
      do k=1,nr
       do j=1,jmax
        do i=1,imax
         if ( k .eq. 1 ) then
           if ( b3d(i,j,1) .ne. 0.d0 ) then
             recvar = 1.d0 / b3d(i,j,1)
             c3d_prime(i,j,1) = c3d(i,j,1)*recvar
             y3d_prime(i,j,1) = y3d(i,j,1)*recvar
           else
C             c3d_prime(i,j,1) = 0.d0
C             y3d_prime(i,j,1) = 0.d0
             call mexErrMsgTxt(
     &         'ERROR: forward sweep of implicit mixing matrix solver.')
           endif
         else
           tmpvar = b3d(i,j,k) - a3d(i,j,k)*c3d_prime(i,j,k-1)
           if ( tmpvar .ne. 0.d0 ) then
             recvar = 1.d0 / tmpvar
             c3d_prime(i,j,k) = c3d(i,j,k)*recvar
             y3d_prime(i,j,k) = ( y3d(i,j,k)
     &                          - a3d(i,j,k)*y3d_prime(i,j,k-1)
     &                          )*recvar
           else
C             c3d_prime(i,j,k) = 0.d0
C             y3d_prime(i,j,k) = 0.d0
             call mexErrMsgTxt(
     &         'ERROR: forward sweep of implicit mixing matrix solver.')
           endif
         endif
        enddo
       enddo
      enddo

C     backward sweep
      do k=nr,1,-1
       do j=1,jmax
        do i=1,imax
         if ( k .eq. nr ) then
          y3d_update(i,j,k) = y3d_prime(i,j,k)
         else
          y3d_update(i,j,k) = y3d_prime(i,j,k)
     &                      - c3d_prime(i,j,k)*y3d_update(i,j,k+1)
         endif
        enddo
       enddo
      enddo

C     update array
      do k=1,nr
       do j=1,jmax
        do i=1,imax
         y3d(i,j,k) = y3d_update(i,j,k)
        enddo
       enddo
      enddo
      return
      end

