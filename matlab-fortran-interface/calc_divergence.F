#include "fintrf.h"
C-----------------------------------------------------------------------------------------
C     MIT_CALC_DIVERGENCE calculates transport divergence and 
C                         integrates with depth.
C     [hflux,rflux]=mit_calc_divergence(dz,volc,kave,mask,
C                                       uflux,vflux,wflux);
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
      real*8, allocatable, dimension(:) :: kin
      real*8, allocatable, dimension(:) :: dz
      real*8, allocatable, dimension(:) :: volin
      real*8, allocatable, dimension(:) :: mskin
      real*8, allocatable, dimension(:) :: ufin
      real*8, allocatable, dimension(:) :: vfin
      real*8, allocatable, dimension(:) :: wfin

C     Allocatable arrays for fortran data inputs to computational routine
      real*8, allocatable, dimension(:,:,:)   :: kave
      real*8, allocatable, dimension(:,:,:)   :: volc
      real*8, allocatable, dimension(:,:,:)   :: mask
      real*8, allocatable, dimension(:,:,:,:) :: uflux
      real*8, allocatable, dimension(:,:,:,:) :: vflux
      real*8, allocatable, dimension(:,:,:,:) :: wflux      

C     Allocatable arrays for fortran data outputs  
      real*8, allocatable, dimension(:,:,:) :: hflux
      real*8, allocatable, dimension(:,:,:) :: rflux
      
C---------------VALIDATE------------INPUTS--------------------------------
C     Check for proper number of arguments. 
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Calculating Flux Divergence')
      call mexWarnMsgTxt('Check for proper number of arguments.')

      write(line,*)'Number of inputs is ',nrhs
      k=mexPrintf(line//achar(13)) 
      write(line,*)'Number of outputs is ',nlhs
      k=mexPrintf(line//achar(13)) 
#endif
      if(nrhs .ne. 7) then
         call mexErrMsgTxt(
     &'ERROR: 7 inputs required (dz,volc,kave,mask,uflux,vflux,wflux)')
      elseif(nlhs .gt. 2) then
         call mexErrMsgTxt(
     &'ERROR: Too many output arguments (hflux, rflux).')
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
      if (ndim .lt. 4) then
         dims(4)=1.d0 
         ndim   =ndim+1
      endif

#ifdef DEBUG_MESSAGES
      write(line,*)'ndims: ',ndim
      k=mexPrintf(line//achar(13)) 
      write(line,*)'input dims:\n ',dims(1),'\nx',dims(2),
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
      write(line,*)'output dims:\n ',fdims(1),'\nx',fdims(2),
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
      call mexWarnMsgTxt('allocate volin array dimensions')  
#endif
      allocate (volin(mxGetNumberOfElements(prhs(2))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate kin array dimensions'  )
#endif
      allocate (kin  (mxGetNumberOfElements(prhs(3))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate mskin array dimensions')   
#endif
      allocate (mskin(mxGetNumberOfElements(prhs(4))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate ufin array dimensions' )  
#endif
      allocate (ufin (mxGetNumberOfElements(prhs(5))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate vfin array dimensions' )    
#endif
      allocate (vfin (mxGetNumberOfElements(prhs(6))))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate wfin array dimensions' )    
#endif
      allocate (wfin (mxGetNumberOfElements(prhs(7))))
 
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
      call mexWarnMsgTxt('allocate volc array dimensions')    
#endif
      allocate (volc(dims(1),dims(2),dims(3)))
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
      
C Allocate the output from the computational function      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate hflux array dimensions')    
#endif
      allocate (hflux(fdims(1),fdims(2),fdims(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('allocate rflux array dimensions')    
#endif
      allocate (rflux(fdims(1),fdims(2),fdims(3)))

C     Fill Fortran arrays from the input argument pointers.
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill arrays from the input.') 
      call mexWarnMsgTxt('Fill dz array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(1)), dz, 
     & mxGetNumberOfElements(prhs(1)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill volin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), volin, 
     & mxGetNumberOfElements(prhs(2)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill kin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), kin, 
     & mxGetNumberOfElements(prhs(3)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill mskin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(4)), mskin, 
     & mxGetNumberOfElements(prhs(4)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill ufin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(5)), ufin, 
     & mxGetNumberOfElements(prhs(5)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill vfin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(6)), vfin, 
     & mxGetNumberOfElements(prhs(6)))
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Fill wfin array from input')       
#endif
      call mxCopyPtrToReal8(mxGetPr(prhs(7)), wfin, 
     & mxGetNumberOfElements(prhs(7)))

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Create arrays for the output arguments.')  
#endif
      classid = mxClassIDFromClassName('double')
      complexflag = 0  
      do iter=1,nlhs
         plhs(iter) = mxCreateNumericArray(ndim-1, fdims, 
     &          classid, 
     &          complexflag)
      enddo

C Check arguments for NANs in the input
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check arguments for NANs')    
      call mexWarnMsgTxt('Check dz')    
#endif
      if (mxIsNaN(dz) .eq. 1) then
C       dz is pretty important, so error if NaN
        call mexErrMsgTxt('ERROR: Input dz contains a NaN.')
      endif

#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check volin')
#endif      
      if (mxIsNaN(volin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input volin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(2))
        if( mxIsNaN(volin(iter)) .eq. 1 ) then
         volin(iter)=0 .d 0
        endif
       enddo
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
       do iter=1,mxGetNumberOfElements(prhs(3))
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
       do iter=1,mxGetNumberOfElements(prhs(4))
        if( mxIsNaN(mskin(iter)) .eq. 1 ) then
         mskin(iter)=0 .d 0
        endif
       enddo
      endif
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check ufin')
#endif      
      if (mxIsNaN(ufin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input ufin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(5))
        if( mxIsNaN(ufin(iter)) .eq. 1 ) then
         ufin(iter)=0 .d 0
        endif
       enddo
      endif
      
#ifdef DEBUG_MESSAGES
      call mexWarnMsgTxt('Check vfin')
#endif          
      if (mxIsNaN(vfin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input vfin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero
       do iter=1,mxGetNumberOfElements(prhs(6))
        if( mxIsNaN(vfin(iter)) .eq. 1 ) then
         vfin(iter)=0 .d 0
        endif
       enddo
      endif
      
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Check wfin')
#endif            
      if (mxIsNaN(wfin) .eq. 1) then
#ifdef DEBUG_MESSAGES
         call mexWarnMsgTxt(
     &      'Input wfin contains a NaN.')
#endif
C Find out where the NaN is and set it to zero  
       do iter=1,mxGetNumberOfElements(prhs(7))
        if( mxIsNaN(wfin(iter)) .eq. 1 ) then
         wfin(iter)=0 .d 0
        endif
       enddo
      endif
        
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Reshape input arrays') 
#endif   
      kave  = reshape(kin  ,(/dims(1),dims(2),dims(4)/))
      volc  = reshape(volin,(/dims(1),dims(2),dims(3)/))
      mask  = reshape(mskin,(/dims(1),dims(2),dims(3)/))
      uflux = reshape(ufin ,(/dims(1),dims(2),dims(3),dims(4)/))
      vflux = reshape(vfin ,(/dims(1),dims(2),dims(3),dims(4)/))
      wflux = reshape(wfin ,(/dims(1),dims(2),dims(3),dims(4)/))

C Fill the output arrays with zeros
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Fill output arrays with zeros') 
#endif   
      do t=1,fdims(3)
       do j=1,fdims(2)
        do i=1,fdims(1)
         hflux(i,j,t)=0.D0
         rflux(i,j,t)=0.D0
        enddo
       enddo  
      enddo

#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Call computational subroutine')  
#endif  
      call calc_flux_divergence(dims(1),dims(2),dims(3),dims(4),
     &                kave,dz,volc,mask,
     &                uflux,vflux,wflux,
     &                hflux,rflux)   
     
C     Load the data into the pointers, which are the output to MATLAB.
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('Load the data into the output pointers')    
#endif
      size=fdims(1)*fdims(2)*fdims(3)
      call mxCopyReal8ToPtr(hflux, mxGetPr(plhs(1)),size)
      call mxCopyReal8ToPtr(rflux, mxGetPr(plhs(2)),size)  
      
      return
      end

CC------------------------COMPUTATIONAL SUBROUTINE------------------------
      
      subroutine calc_flux_divergence(imax,jmax,nr,nt,klev,dz,volc,
     &                       maskc,fzon,fmer,fvert,hflux,rflux)
      
      implicit none
      
      mwIndex imax,jmax,nr,nt,i,j,k,t
      mwIndex ip1,jp1,km1,kp1
      real*8 dz   (nr)
      real*8 volc (imax,jmax,nr)
      real*8 recip_volc(imax,jmax,nr)
      real*8 maskc(imax,jmax,nr)
      real*8 klev (imax,jmax,nt)
      real*8 hdiv (imax,jmax,nr)
      real*8 rdiv (imax,jmax,nr)
      real*8 hflux(imax,jmax,nt)
      real*8 rflux(imax,jmax,nt)
      real*8 fzon (imax,jmax,nr,nt)
      real*8 fmer (imax,jmax,nr,nt)
      real*8 fvert(imax,jmax,nr,nt)

#ifdef DEBUG_MESSAGES 
      mwIndex, external :: mexPrintf
      mwIndex msg
      character(len=250) :: line
#endif

C From config_summary
C rkSign =   /* index orientation relative to vertical coordinate */
C -1.000000000000000E+00
      real*8 rkSign
      PARAMETER(rkSign=-1.D0)
      
C     Initialize fluxes to zero and calculate reciprocals
#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('zero the outputs')    
      write(line,*)'dimensions:\ni:',imax,'\nj:',jmax,
     &                        '\nk',nr,'\nt:',nt 
      msg=mexPrintf(line//achar(13)) 
#endif
      do t=1,nt
       do j=1,jmax
        do i=1,imax  
         hflux(i,j,t)=0.D0
         rflux(i,j,t)=0.D0
        enddo
       enddo  
      enddo

#ifdef DEBUG_MESSAGES 
      call mexWarnMsgTxt('calculate 1/volume')    
#endif
      do j=1,jmax
       do i=1,imax
        do k=1,nr 
         recip_volc(i,j,k)=1.D0/volc(i,j,k)
        enddo
       enddo
      enddo

C--   Divergence of fluxes
#ifdef DEBUG_MESSAGES 
	  call mexWarnMsgTxt('calculate divergence')    
#endif
      do t=1,nt
       do k=1,nr
        do j=1,jmax
         do i=1,imax
          if (i .eq. imax) THEN
           ip1=1
          else
           ip1=i+1
          endif
          
          jp1=MIN(j+1,jmax)
          
          km1  = MAX(1 ,k-1)
          kp1  = MIN(nr,k+1)
        
          hdiv(i,j,k) = -maskc(i,j,k)*recip_volc(i,j,k)
     &            *((fzon(ip1, j ,k,t)-fzon(i,j,k,t))
     &             +(fmer( i ,jp1,k,t)-fmer(i,j,k,t))
     &             )
       
          rdiv(i,j,k) = -maskc(i,j,k)*recip_volc(i,j,k)
     &             *(fvert(i,j,kp1,t)-fvert(i,j,k,t))*rkSign
         enddo
        enddo
       enddo

C     Integrate vertically 
#ifdef DEBUG_MESSAGES 
      if (t .eq. 1) then
         call mexWarnMsgTxt('integrate vertically')  
      endif
#endif  
       do j=1,jmax
        do i=1,imax  
         do k=1,INT(klev(i,j,t))    
          hflux(i,j,t)=hflux(i,j,t)
     &      +(hdiv(i,j,k)*dz(k)*maskc(i,j,k))
     
          rflux(i,j,t)=rflux(i,j,t)
     &      +(rdiv(i,j,k)*dz(k)*maskc(i,j,k))
         enddo
        enddo
       enddo
      enddo
            
      return
      end
