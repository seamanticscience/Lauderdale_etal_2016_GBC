C
C      if (advectionscheme == 1) then
CC        first order upwind      
C         call advx_uwcntr_2nd(dims(1),dims(2),dims(3),dims(4),
C     &           rax,mskc,msku,
C     &           uvel,tracer,uflux)
C      elseif (advectionscheme == 2) then
CC        centre second order      
C         call advx_cntr_2nd(dims(1),dims(2),dims(3),dims(4),
C     &           rax,mskc,msku,
C     &           uvel,tracer,uflux)
C      elseif (advectionscheme == 3) then
CC        third order upwind biased
C         call advx_uw_3rd(dims(1),dims(2),dims(3),dims(4),
C     &           rax,mskc,msku,deltat,
C     &           dx,uvel,tracer,uflux)
C      elseif (advectionscheme == 4) then
CC        centre fourth order
C         call advx_cntr_4th(dims(1),dims(2),dims(3),dims(4),
C     &           rax,mskc,msku,
C     &           uvel,tracer,uflux)
C      elseif (advectionscheme == 7) then
CC     7th Order One Step method with Monotonicity Preserving Limiter
C         call advx_mpl_7th(dims(1),dims(2),dims(3),dims(4),
C     &           rax,mskc,msku,deltat,
C     &           dx,uvel,tracer,uflux)
C      elseif (advectionscheme == 20 .or. advectionscheme == 21) then
CC        second order direct-space-and-time, first upwind
C         call advx_dst_2nd(dims(1),dims(2),dims(3),dims(4),
C     &           advectionscheme,    
C     &           rax,mskc,msku,deltat,
C     &           dx,uvel,tracer,uflux)
C      elseif (advectionscheme == 30) then
CC        third order direct-space-and-time
C         call advx_dst_3rd(dims(1),dims(2),dims(3),dims(4),
C     &           rax,mskc,msku,deltat,
C     &           dx,uvel,tracer,uflux)
C      elseif (advectionscheme == 33) then
CC        third order direct-space-and-time (flux limited)    
C         call advx_dstfl_3rd(dims(1),dims(2),dims(3),dims(4),
C     &           rax,mskc,msku,deltat,
C     &           dx,uvel,tracer,uflux)
C      elseif (advectionscheme == 77) then
CC        centre second order (flux liimited)    
C         call advx_cntrfl_2nd(dims(1),dims(2),dims(3),dims(4),
C     &           rax,mskc,msku,deltat,
C     &           dx,uvel,tracer,uflux)
C      else
C       call mexErrMsgTxt('ERROR: advectionscheme not recognized, values are:
C     &  1,2,3,4,7,20,21,30,33,77')
C      endif

       subroutine advx_cntr_2nd(
     &           nx,ny,nr,nt,advectionscheme,rax,
     &           maskc,masku,deltat,dx,uvel,tracer,uflux)

C calculates the centered second-order area integrated zonal flux due to advection of a tracer using
      implicit none
      include "offline_consts.h"
      include "offline_zonal_advection.h"
C Local variable definitions
      integer i,j,k,t,im1,ip1
      real*8 deltat
      
         do t=1,nt
          do k=1,nr
           do j=1,ny
            do i=1,nx
             if (i .eq. 1) then
                im1=nx
             else
                im1=i-1
             endif
             
             if (i .eq. nx) then
                ip1=1
             else
                ip1=i+1
             endif
             
             uflux(i,j,k,t) = uvel(i,j,k,t)*rax(i,j,k)
     &        *(tracer(i,j,k,t)+tracer(im1,j,k,t))*half             

C             uflux(i,j,k,t) = maskc(i,j,k)*maskc(im1,j,k )*
C     &                     masku(i,j,k)*masku(im1,j,k )*         
C     &                     uvel(i,j,k,t)*rax(i,j,k)
C     &        *(tracer(i,j,k,t)+tracer(im1,j,k,t))*half
            enddo
           enddo
          enddo
         enddo

      return
      end

      subroutine advx_cntr_4th(
     &           nx,ny,nr,nt,advectionscheme,rax,
     &           maskc,masku,deltat,dx,uvel,tracer,uflux)

C calculates the centered fourth-order area integrated zonal flux due to advection of a tracer
      implicit none
      include "offline_consts.h"
      include "offline_zonal_advection.h"
C Local variable definitions
      integer i,j,k,t,im1,im2,ip1
      real*8 deltat
      real*8 rjm,rj,rjp,rjjm,rjjp
      real*8 utrans
      
      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx           
          if (i .eq. 2) then
           im1=1
           im2=nx
          elseif (i .eq. 1) then
           im1=nx
           im2=nx-1
          else
           im1=i-1
           im2=i-2
          endif

          if (i .eq. nx) then
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
        
          utrans=uvel(i,j,k,t)*rax(i,j,k)

C          uflux(i,j,k,t) = maskc( i ,j,k)*maskc(im1,j,k)*
C     &                 *masku( i ,j,k)*masku(im1,j,k)*     
          uflux(i,j,k,t) = 
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

      subroutine advx_cntrfl_2nd(
     &           nx,ny,nr,nt,advectionscheme,rax,
     &           maskc,masku,deltat,dx,uvel,tracer,uflux)
     
C calculates the centered second-order area integrated zonal flux due to 
C advection of a tracer  with a flux limiter
      implicit none
      include "offline_consts.h"
      include "offline_zonal_advection.h"
C Local variable definitions
      integer i,j,k,t,im1,im2
      real*8 deltat
      real*8 recipx(nx,ny)
      real*8 rjm,rj,rjp,ucfl,cr,lcr
      real*8 utrans

      do j=1,ny
       do i=1,nx
          recipx(i,j)=1/dx(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx           
          if (i .eq. 2) then
           im1=1
           im2=nx
          elseif (i .eq. 1) then
           im1=nx
           im2=nx-1
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
          utrans=uvel(i,j,k,t)*rax(i,j,k)

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
        
C        uflux(i,j,k,t) = masku( i ,j,k)*masku(im1,j,k)*
C     &                maskc( i ,j,k)*maskc(im1,j,k)* 
        uflux(i,j,k,t) =   
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
     &           nx,ny,nr,nt,advectionscheme,rax,
     &           maskc,masku,deltat,dx,uvel,tracer,uflux)

C  calculates the second-order direct space and time (dst-2/lax-wendroff, or simple 
C  first order upwind scheme) area integrated zonal flux due to advection of a tracer 
      implicit none
      include "offline_consts.h"
      include "offline_zonal_advection.h"
C Local variable definitions
      integer i,j,k,t,im1
      real*8 deltat
      real*8 ucfl, xlimit, uabs, utrans
      real*8 recipx(nx,ny)
      
      if ( advectionscheme.eq.20 ) then
       xlimit = 1.d0
      else
       xlimit = 0.d0
      endif 

      do j=1,ny
       do i=1,nx
          recipx(i,j)=1/dx(i,j)
       enddo
      enddo
         
      do t=1,nt          
       do k=1,nr        
        do j=1,ny
         do i=1,nx
          
          if (i .eq. 1) then
           im1=nx
          else
           im1=i-1
          endif
      
          ucfl = abs( uvel(i,j,k,t)*deltat*recipx(i,j) )

          utrans=uvel(i,j,k,t)*rax(i,j,k)
        
          uabs = abs(utrans)*( 1.d0 - xlimit*(1.d0 - ucfl) )
          
          uflux(i,j,k,t) = 
     &          ( utrans+uabs )* 0.5d0 * tracer(im1,j,k,t)
     &        + ( utrans-uabs )* 0.5d0 * tracer(i  ,j,k,t)
         enddo
        enddo
       enddo
      enddo
      
      return
      end

      subroutine advx_dst_3rd(
     &           nx,ny,nr,nt,advectionscheme,rax,
     &           maskc,masku,deltat,dx,uvel,tracer,uflux)
   
C  calculates the 3rd-order direct space and time (dst-3) area integrated zonal flux 
C  due to advection of a tracer.
      implicit none
      include "offline_consts.h"
      include "offline_zonal_advection.h"
C Local variable definitions
      integer i,j,k,t,im1,im2
      real*8 deltat
      real*8 recipx(nx,ny)
      real*8 rjm,rj,rjp,d0,d1
      real*8 utrans,ucfl
      
      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx
          uflux(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        
      
      do j=1,ny
       do i=1,nx
          recipx(i,j)=1/dx(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx           
          if (i .eq. 2) then
           im1=1
           im2=nx
          elseif (i .eq. 1) then
           im1=nx
           im2=nx-1
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
          utrans=uvel(i,j,k,t)*rax(i,j,k)
          
          d0=(2.d0-ucfl)*(1.d0-ucfl)*onesixth
          d1=(1.d0-ucfl*ucfl)*onesixth

          uflux(i,j,k,t)=maskc( i ,j,k)*maskc(im1,j,k)
C     &              *masku( i ,j,k)*masku(im1,j,k)
C          uflux(i,j,k,t)= (
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

      subroutine advx_dstfl_3rd(
     &           nx,ny,nr,nt,advectionscheme,rax,
     &           maskc,masku,deltat,dx,uvel,tracer,uflux)
     
C  calculates the 3rd-order direct space and time (dst-3) area integrated zonal flux 
C  due to advection of a tracer with flux limiting.
      implicit none
      include "offline_consts.h"
      include "offline_zonal_advection.h"
C Local variable definitions
      integer i,j,k,t,im1,im2
      real*8 deltat
      real*8 recipx(nx,ny)
      real*8 rjm,rj,rjp,d0,d1,psip,psim
      real*8 thetap,thetam
      real*8 utrans,ucfl

      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx
          uflux(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do j=1,ny
       do i=1,nx
          recipx(i,j)=1/dx(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx           
          if (i .eq. 2) then
           im1=1
           im2=nx
          elseif (i .eq. 1) then
           im1=nx
           im2=nx-1
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
          utrans=uvel(i,j,k,t)*rax(i,j,k)
          
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

           uflux(i,j,k,t)=maskc( i ,j,k)*maskc(im1,j,k)
C     &                *masku( i ,j,k)*masku(im1,j,k)
C           uflux(i,j,k,t)=(
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

      subroutine advx_mpl_7th(
     &           nx,ny,nr,nt,advectionscheme,rax,
     &           maskc,masku,deltat,dx,uvel,tracer,uflux)

C calculates the 7th order dst area integrated zonal flux due to 
C advection of a tracer with a monotone preserving limiter
      implicit none
      include "offline_consts.h"
      include "offline_zonal_advection.h"
C Local variable definitions
      integer i,j,k,t
      integer im1,im2,im3,im4,ip1,ip2,ip3
      real*8 deltat
      real*8 recipx(nx,ny)
      real*8 cfl,psi,utrans
      real*8 fac,delip,deli,phi,rp1h,rp1h_cfl
      real*8 recip_delip, recip_deli
      real*8 qippp,qipp,qip,qi,qim,qimm,qimmm
      real*8 mskipp,mskip,mski,mskim,mskimm,mskimmm
      real*8 d2,d2p1,d2m1,a,b,c,d
      real*8 dp1h,dm1h, phimd,philc,phimin,phnx
      real*8 delm,delp,delmm,delpp,delmmm,delppp
      real*8 del2mm,del2m,del2,del2p,del2pp
      real*8 del3mm,del3m,del3p,del3pp
      real*8 del4m,del4,del4p
      real*8 del5m,del5p
      real*8 del6
      
      do j=1,ny
       do i=1,nx
          recipx(i,j)=1/dx(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx        
      
          if (i .eq. 4) then
           im1=i-1
           im2=1-2
           im3=1-3
           im4=nx
          elseif (i .eq. 3) then
           im1=i-1
           im2=i-2
           im3=nx
           im4=nx-1          
          elseif (i .eq. 2) then
           im1=1
           im2=nx
           im3=nx-1
           im4=nx-2
          elseif (i .eq. 1) then
           im1=nx
           im2=nx-1 
           im3=nx-2
           im4=nx-3
          else
           im1=i-1
           im2=i-2
           im3=i-3
           im4=i-4
          endif
          
          if (i .eq. nx ) then
           ip1=1
           ip2=2
           ip3=3
          elseif (i .eq. nx-1 ) then
           ip1=i+1
           ip2=1
           ip3=2
          elseif (i .eq. nx-2 ) then
           ip1=i+1
           ip2=i+2
           ip3=1
          else
           ip1=i+1
           ip2=i+2
           ip3=i+3   
          endif
          
          cfl = abs(uvel(i,j,k,t)*deltat*recipx(i,j) )
          utrans=uvel(i,j,k,t)*rax(i,j,k)
        
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
           phnx = min( max(2.d0/(1.d0-cfl),phimd),
     &                 max(0.d0,2.d0*rp1h_cfl,philc) )
           phi = max(phimin,min(phi,phnx))

           psi = phi * 0.5d0 * (1.d0 - cfl)
           uflux(i,j,k,t) = utrans*( qi + psi*delip )
          else
           uflux(i,j,k,t) = 0.d0
          endif
         enddo 
        enddo
       enddo
      enddo

      return
      end

      subroutine advx_uw_3rd(
     &           nx,ny,nr,nt,advectionscheme,rax,
     &           maskc,masku,deltat,dx,uvel,tracer,uflux)

C calculates the third-order upwind biased area integrated zonal flux due to 
C advection of a tracer
      implicit none
      include "offline_consts.h"
      include "offline_zonal_advection.h"
C Local variable definitions
      integer i,j,k,t,im1,im2,ip1
      real*8 deltat
      real*8 rjm,rj,rjp,rjjp,rjjm
      real*8 utrans

      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx
          uflux(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx           
          if (i .eq. 2) then
           im1=1
           im2=nx
          elseif (i .eq. 1) then
           im1=nx
           im2=nx-1
          else
           im1=i-1
           im2=i-2
          endif
          
          if (i .eq. nx) then
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
        
          utrans=uvel(i,j,k,t)*rax(i,j,k)
        
          uflux(i,j,k,t) = maskc( i ,j,k)*
C     &                  maskc(im1,j,k)*
C     &                  masku( i ,j,k)*masku(im1,j,k)*     
C          uflux(i,j,k,t) =
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

      subroutine advx_uwcntr_2nd(
     &           nx,ny,nr,nt,advectionscheme,rax,
     &           maskc,masku,deltat,dx,uvel,tracer,uflux)

C calculates the center second upwind biased area integrated zonal flux due to 
C advection of a tracer
      implicit none
      include "offline_consts.h"
      include "offline_zonal_advection.h"
C Local variable definitions
      integer i,j,k,t,im1
      real*8 deltat
      
         do t=1,nt
          do k=1,nr
           do j=1,ny
            do i=1,nx
             if (i .eq. 1) then
                im1=nx
             else
                im1=i-1
             endif
             
             if (uvel(i,j,k,t) .ge. 0.d0) then
C                uflux(i,j,k,t) = masku(i,j,k)*masku(im1,j,k )*
C     &                        maskc(i,j,k)*maskc(im1,j,k )* 
                uflux(i,j,k,t) = 
     &                        uvel(i,j,k,t)*tracer(im1,j,k,t)*rax(i,j,k)
             else
C                uflux(i,j,k,t) = masku(i,j,k )*maskc(i,j,k)*
                uflux(i,j,k,t) = 
     &                        uvel(i,j,k,t)*tracer( i ,j,k,t)*rax(i,j,k)
             endif

            enddo
           enddo
          enddo
         enddo

      return
      end

CCC----------------ADVECTION IN Y DIRECTION----------------   
C      if (advectionscheme == 1) then
CC        first order upwind     
C         call advy_uwcntr_2nd(dims(1),dims(2),dims(3),dims(4),
C     &           ray,mskc,mskv,
C     &           vvel,tracer,vflux)
C      elseif (advectionscheme == 2) then
CC        centre second order      
C         call advy_cntr_2nd(dims(1),dims(2),dims(3),dims(4),
C     &           ray,mskc,mskv,
C     &           vvel,tracer,vflux)
C      elseif (advectionscheme == 3) then
CC        third order upwind biased
C         call advy_uw_3rd(dims(1),dims(2),dims(3),dims(4),
C     &           ray,mskc,mskv,deltat,
C     &           dy,vvel,tracer,vflux)
C      elseif (advectionscheme == 4) then
CC        centre fourth order
C         call advy_cntr_4th(dims(1),dims(2),dims(3),dims(4),
C     &           ray,mskc,mskv,
C     &           vvel,tracer,vflux)
C      elseif (advectionscheme == 7) then
CC     7th Order One Step method with Monotonicity Preserving Limiter
C         call advy_mpl_7th(dims(1),dims(2),dims(3),dims(4),
C     &           ray,mskc,mskv,deltat,
C     &           dy,vvel,tracer,vflux)
C      elseif (advectionscheme == 20 .or. advectionscheme == 21) then
CC        second order direct-space-and-time, first upwind
C         call advy_dst_2nd(dims(1),dims(2),dims(3),dims(4),
C     &           advectionscheme,    
C     &           ray,mskc,mskv,deltat,
C     &           dy,vvel,tracer,vflux)
C      elseif (advectionscheme == 30) then
CC        third order direct-space-and-time
C         call advy_dst_3rd(dims(1),dims(2),dims(3),dims(4),
C     &           ray,mskc,mskv,deltat,
C     &           dy,vvel,tracer,vflux)
C      elseif (advectionscheme == 33) then
CC        third order direct-space-and-time (flux limited)    
C         call advy_dstfl_3rd(dims(1),dims(2),dims(3),dims(4),
C     &           ray,mskc,mskv,deltat,
C     &           dy,vvel,tracer,vflux)
C      elseif (advectionscheme == 77) then
CC        centre second order (flux liimited)    
C         call advy_cntrfl_2nd(dims(1),dims(2),dims(3),dims(4),
C     &           ray,mskc,mskv,deltat,
C     &           dy,vvel,tracer,vflux)
C      else
C       call mexErrMsgTxt('ERROR: advectionscheme not recognized, values are:
C     &  1,2,3,4,7,20,21,30,33,77')    
C      endif

      subroutine advy_cntr_2nd(
     &           nx,ny,nr,nt,advectionscheme,ray,
     &           maskc,maskv,deltat,dy,vvel,tracer,vflux)

C calculates the centered second-order area integrated meridional flux due to advection of a tracer
      implicit none
      include "offline_consts.h"
      include "offline_merid_advection.h"
C Local variable definitions
      integer i,j,k,t,jm1,jp1
      real*8 deltat
      
         do t=1,nt
          do k=1,nr
           do j=1,ny
            do i=1,nx
             jm1=max(j-1, 1  )
             jp1=min(j+1,ny)
             
             vflux(i,j,k,t) = vvel(i,j,k,t)*ray(i,j,k)
     &             *(tracer(i,j,k,t)+tracer(i,jm1,k,t))*half

C             vflux(i,j,k,t) = maskc(i,j,k)*maskc(i,jm1,k)*
C     &                     maskv(i,j,k)*maskv(i,jm1,k)*
C     &                     vvel(i,j,k,t)*ray(i,j,k)
C     &             *(tracer(i,j,k,t)+tracer(i,jm1,k,t))*half
            enddo
           enddo  
          enddo
         enddo

      return
      end

      subroutine advy_cntr_4th(
     &           nx,ny,nr,nt,advectionscheme,ray,
     &           maskc,maskv,deltat,dy,vvel,tracer,vflux)

C calculates the centered fourth-order area integrated meridional flux due to advection of a tracer
      implicit none
      include "offline_consts.h"
      include "offline_merid_advection.h"
C Local variable definitions     
      integer i,j,k,t,jm1,jm2,jp1
      real*8 deltat
      real*8 rjm,rj,rjp,rjjm,rjjp
      real*8 vtrans
      
      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx           
          jp1=min(j+1,ny)
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
          
          vtrans=vvel(i,j,k,t)*ray(i,j,k)

C          vflux(i,j,k,t) = maskc(i, j ,k)*maskc(i,jm1,k)*
C     &                 *maskv(i, j ,k)*maskv(i,jm1,k)*     
          vflux(i,j,k,t) = 
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

      subroutine advy_cntrfl_2nd(
     &           nx,ny,nr,nt,advectionscheme,ray,
     &           maskc,maskv,deltat,dy,vvel,tracer,vflux)

C calculates the centered second-order area integrated meridional flux due to 
C advection of a tracer with a flux limiter
      implicit none
      include "offline_consts.h"
      include "offline_merid_advection.h"
C Local variable definitions    
      integer i,j,k,t,jm1,jm2,jp1
      real*8 deltat
      real*8 recipy(nx,ny)
      real*8 rjm,rj,rjp,vcfl,cr,lcr
      real*8 vtrans

      do j=1,ny
       do i=1,nx
          recipy(i,j)=1/dy(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx           
          jp1=min(j+1,ny)
          jm1=max(j-1, 1  )
          jm2=max(j-2, 1  )
          
          rjp=(tracer(i,jp1,k,t)-tracer(i, j ,k,t))
C     &     *maskv(i,jp1,k)*maskc(i,jp1,k)
          rj =(tracer(i, j ,k,t)-tracer(i,jm1,k,t))
C     &     *maskv(i, j ,k)*maskc(i, j ,k)
          rjm=(tracer(i,jm1,k,t)-tracer(i,jm2,k,t))
C     &     *maskv(i,jm1,k)*maskc(i,jm1,k)

          vcfl = abs( vvel(i,j,k,t)*deltat*recipy(i,j) )
          vtrans=vvel(i,j,k,t)*ray(i,j,k)

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
             
C         vflux(i,j,k,t) =maskv(i, j ,k)*maskv(i,jm1,k)*
C     &                maskc(i, j ,k)*maskc(i,jm1,k)* 
         vflux(i,j,k,t) =   
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
     &           nx,ny,nr,nt,advectionscheme,ray,
     &           maskc,maskv,deltat,dy,vvel,tracer,vflux)

C  calculates the second-order direct space and time (dst-2/lax-wendroff, or simple 
C  first order upwind scheme) area integrated meridional flux due to advection of a tracer 
      implicit none
      include "offline_consts.h"
      include "offline_merid_advection.h"
C Local variable definitions
      integer i,j,k,t,jm1
      real*8 deltat
      real*8 vcfl, ylimit, vabs, vtrans
      real*8 recipy(nx,ny)

      if ( advectionscheme.eq.20 ) then
       ylimit = 1.d0
      else
       ylimit = 0.d0
      endif
    
      do j=1,ny
       do i=1,nx
          recipy(i,j)=1/dy(i,j)
       enddo
      enddo
      
      do t=1,nt
       do k=1,nr
        do j=1,ny
         
         jm1=max(j-1, 1  )
         
         do i=1,nx
          
          vcfl = abs( vvel(i,j,k,t)*deltat*recipy(i,j) )
          vtrans=vvel(i,j,k,t)*ray(i,j,k)

          vabs = abs(vtrans)*( 1.d0 - ylimit*(1.d0 - vcfl) )
          vflux(i,j,k,t) = ( vtrans+vabs )* 0.5d0 * tracer(i,jm1,k,t)
     &                + ( vtrans-vabs )* 0.5d0 * tracer(i,j  ,k,t)

         enddo
        enddo
       enddo
      enddo
      
      return
      end

      subroutine advy_dst_3rd(
     &           nx,ny,nr,nt,advectionscheme,ray,
     &           maskc,maskv,deltat,dy,vvel,tracer,vflux)

C  calculates the 3rd-order direct space and time (dst-3) area integrated meridional flux 
C  due to advection of a tracer.
      implicit none
      include "offline_consts.h"
      include "offline_merid_advection.h"
C Local variable definitions
      integer i,j,k,t,jm1,jm2,jp1
      real*8 deltat
      real*8 recipy(nx,ny)
      real*8 rjm,rj,rjp,d0,d1
      real*8 vtrans,vcfl

      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx
          vflux(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do j=1,ny
       do i=1,nx
          recipy(i,j)=1/dy(i,j)
       enddo
      enddo
            
      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx           
          jp1=min(j+1,ny)
          jm1=max(j-1, 1  )
          jm2=max(j-2, 1  )
          
          rjp=(tracer(i,jp1,k,t)-tracer(i, j ,k,t))
     &     *maskc(i, j ,k)*maskc(i,jp1,k)
          rj =(tracer(i, j ,k,t)-tracer(i,jm1,k,t))
     &     *maskc(i, j ,k)*maskc(i,jm1,k)
          rjm=(tracer(i,jm1,k,t)-tracer(i,jm2,k,t))
     &     *maskc(i,jm1,k)*maskc(i,jm2,k)

          vcfl = abs( vvel(i,j,k,t)*deltat*recipy(i,j) )
          vtrans=vvel(i,j,k,t)*ray(i,j,k)
          
          d0=(2.d0-vcfl)*(1.d0-vcfl)*onesixth
          d1=(1.d0-vcfl*vcfl)*onesixth

          vflux(i,j,k,t)=maskc(i, j ,k)*maskc(i,jm1,k)
C     &               *maskv(i, j ,k)*maskv(i,jm1,k)   
C          vflux(i,j,k,t)= (  
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

      subroutine advy_dstfl_3rd(
     &           nx,ny,nr,nt,advectionscheme,ray,
     &           maskc,maskv,deltat,dy,vvel,tracer,vflux)

C  calculates the 3rd-order direct space and time (dst-3) area integrated zonal flux 
C  due to advection of a tracer with flux limiting.
      implicit none
      include "offline_consts.h"
      include "offline_merid_advection.h"
C Local variable definitions
      integer i,j,k,t,jm1,jm2,jp1
      real*8 deltat
      real*8 recipy(nx,ny)
      real*8 rjm,rj,rjp,d0,d1,psip,psim
      real*8 thetap,thetam
      real*8 vtrans,vcfl

      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx
          vflux(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do j=1,ny
       do i=1,nx
          recipy(i,j)=1/dy(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx           
          jp1=min(j+1,ny)
          jm1=max(j-1, 1  )
          jm2=max(j-2, 1  )
          
          rjp=(tracer(i,jp1,k,t)-tracer(i, j ,k,t))
     &     *maskc(i, j ,k)*maskc(i,jp1,k)
          rj =(tracer(i, j ,k,t)-tracer(i,jm1,k,t))
     &     *maskc(i, j ,k)*maskc(i,jm1,k)
          rjm=(tracer(i,jm1,k,t)-tracer(i,jm2,k,t))
     &     *maskc(i,jm1,k)*maskc(i,jm2,k)

          vcfl = abs( vvel(i,j,k,t)*deltat*recipy(i,j) )
          vtrans=vvel(i,j,k,t)*ray(i,j,k)
          
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

           vflux(i,j,k,t)=maskc(i, j ,k)*maskc(i,jm1,k)
C     &                *maskv(i, j ,k)*maskv(i,jm1,k)
C           vflux(i,j,k,t)=(      
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

      subroutine advy_mpl_7th(
     &           nx,ny,nr,nt,advectionscheme,ray,
     &           maskc,maskv,deltat,dy,vvel,tracer,vflux)

C calculates the 7th order dst area integrated meridional flux due to 
C advection of a tracer with a monotone preserving limiter
      implicit none
      include "offline_consts.h"
      include "offline_merid_advection.h"
C Local variable definitions
      integer i,j,k,t
      integer jm1,jm2,jm3,jm4,jp1,jp2,jp3
      real*8 deltat
      real*8 recipy(nx,ny)
      real*8 cfl,psi,vtrans
      real*8 fac,delip,deli,phi,rp1h,rp1h_cfl
      real*8 recip_delip, recip_deli
      real*8 qippp,qipp,qip,qi,qim,qimm,qimmm
      real*8 mskipp,mskip,mski,mskim,mskimm,mskimmm
      real*8 d2,d2p1,d2m1,a,b,c,d
      real*8 dp1h,dm1h, phimd,philc,phimin,phnx
      real*8 delm,delp,delmm,delpp,delmmm,delppp
      real*8 del2mm,del2m,del2,del2p,del2pp
      real*8 del3mm,del3m,del3p,del3pp
      real*8 del4m,del4,del4p
      real*8 del5m,del5p
      real*8 del6

      do j=1,ny
       do i=1,nx
          recipy(i,j)=1/dy(i,j)
       enddo
      enddo
      
      do t=1,nt    
       do k=1,nr
        do j=1,ny
          
          jp1=min(j+1,ny)
          jp2=min(j+2,ny)
          jp3=min(j+3,ny)
          jm1=max(j-1, 1  )
          jm2=max(j-2, 1  )
          jm3=max(j-3, 1  )
          jm4=max(j-4, 1  )
          
         do i=1,nx           

          cfl = abs( vvel(i,j,k,t)*deltat*recipy(i,j) )
          vtrans=vvel(i,j,k,t)*ray(i,j,k)

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
         phnx = min( max(2.d0/(1.d0-cfl),phimd),
     &                 max(0.d0,2.d0*rp1h_cfl,philc) )
         phi = max(phimin,min(phi,phnx))

         psi = phi * 0.5d0 * (1.d0 - cfl)
         vflux(i,j,k,t) = vtrans*( qi + psi*delip )
        else
         vflux(i,j,k,t) = 0.d0
        endif
       enddo
       enddo
       enddo
      enddo

      return
      end

      subroutine advy_uw_3rd(
     &           nx,ny,nr,nt,advectionscheme,ray,
     &           maskc,maskv,deltat,dy,vvel,tracer,vflux)

C calculates the third-order upwind biased area integrated meridional flux due to 
C advection of a tracer
      implicit none
      include "offline_consts.h"
      include "offline_merid_advection.h"
C Local variable definitions
      integer i,j,k,t,jm1,jm2,jp1
      real*8 deltat
      real*8 rjm,rj,rjp,rjjm,rjjp
      real*8 vtrans
      
      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx
          vflux(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx           
          jp1=min(j+1,ny)
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

          vtrans=vvel(i,j,k,t)*ray(i,j,k)
          
          vflux(i,j,k,t) = maskc(i, j ,k)*
C     &                  maskc(i,jm1,k)*
C     &                  maskv(i, j ,k)*maskv(i,jm1,k)*     
C          vflux(i,j,k,t) = 
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

      subroutine advy_uwcntr_2nd(
     &           nx,ny,nr,nt,advectionscheme,ray,
     &           maskc,maskv,deltat,dy,vvel,tracer,vflux)

C calculates the center second upwind biased area integrated meridional flux due to 
C advection of a tracer
      implicit none
      include "offline_consts.h"
      include "offline_merid_advection.h"
C Local variable definitions
      integer i,j,k,t,jm1
      real*8 deltat
      
         do t=1,nt
          do k=1,nr
           do j=1,ny
            do i=1,nx
             jm1=max(j-1, 1  )
 
             if (vvel(i,j,k,t) .ge. 0.0) then
C               vflux(i,j,k,t) = maskv(i,j,k) *maskv(i,jm1,k )*
C     &                       maskc(i,j,k) *maskc(i,jm1,k )*      
               vflux(i,j,k,t) =     
     &                       vvel(i,j,k,t)*tracer(i,jm1,k,t)*ray(i,j,k)
             else
C               vflux(i,j,k,t) = maskv(i,j,k )*maskc(i,j,k)*
               vflux(i,j,k,t) = 
     &                       vvel(i,j,k,t)*tracer(i, j ,k,t)*ray(i,j,k)
             endif
            
            enddo
           enddo  
          enddo
         enddo

      return
      end

CCC----------------ADVECTION IN Z DIRECTION----------------
C      if (advectionscheme == 1) then
CC        first order upwind      
C         call advz_uwcntr_2nd(dims(1),dims(2),dims(3),dims(4),
C     &           rac,mskc,
C     &           wvel,tracer,wflux)
C      elseif (advectionscheme == 2) then
CC        centre second order      
C         call advz_cntr_2nd(dims(1),dims(2),dims(3),dims(4),
C     &           rac,mskc,
C     &           wvel,tracer,wflux)
C      elseif (advectionscheme == 3) then
CC        third order upwind biased
C         call advz_uw_3rd(dims(1),dims(2),dims(3),dims(4),
C     &           rac,mskc,deltat,dz,
C     &           wvel,tracer,wflux)
C      elseif (advectionscheme == 4) then
CC        centre fourth order
C         call advz_cntr_4th(dims(1),dims(2),dims(3),dims(4),
C     &           rac,mskc,
C     &           wvel,tracer,wflux)
C      elseif (advectionscheme == 7) then
CC     7th Order One Step method with Monotonicity Preserving Limiter
C         call advz_mpl_7th(dims(1),dims(2),dims(3),dims(4),
C     &           rac,mskc,deltat,dz,
C     &           wvel,tracer,wflux)
C      elseif (advectionscheme == 20 .or. advectionscheme == 21) then
CC        second order direct-space-and-time, 1st upwind
C         call advz_dst_2nd(dims(1),dims(2),dims(3),dims(4),
C     &           advectionscheme,    
C     &           rac,mskc,deltat,dz,
C     &           wvel,tracer,wflux)
C      elseif (advectionscheme == 30) then
CC        third order direct-space-and-time
C         call advz_dst_3rd(dims(1),dims(2),dims(3),dims(4),
C     &           rac,mskc,deltat,dz,
C     &           wvel,tracer,wflux)
C      elseif (advectionscheme == 33) then
CC        third order direct-space-and-time (flux liimited)    
C         call advz_dstfl_3rd(dims(1),dims(2),dims(3),dims(4),
C     &           rac,mskc,deltat,dz,
C     &           wvel,tracer,wflux)
C      elseif (advectionscheme == 77) then
CC        centre second order (flux liimited)    
C         call advz_cntrfl_2nd(dims(1),dims(2),dims(3),dims(4),
C     &           rac,mskc,deltat,dz,
C     &           wvel,tracer,wflux)
C      else
C       call mexErrMsgTxt('ERROR: advectionscheme not recognized, values are:
C     &  1,2,3,4,7,20,21,30,33,77')    
C      endif
C      
C      return

      subroutine advz_cntr_2nd(
     &           nx,ny,nr,nt,advectionscheme,rac,
     &           maskc,deltat,dz,wvel,tracer,wflux)


C calculates the centered second-order area integrated vertical flux due to advection of a tracer
      implicit none
      include "offline_consts.h"
      include "offline_vert_advection.h"
C Local variable definitions
      integer i,j,k,t,km1,kp1
      real*8 deltat
               
      do t=1,nt    
       do k=1,nr
        km1=max(1 ,k-1)
        kp1=min(nr,k+1)
        
        if ( k.eq.1 .or. k.gt.nr) then
          do j=1,ny
           do i=1,nx
           wflux(i,j,k,t) = 0.0000
          enddo
         enddo
        else
          do j=1,ny
           do i=1,nx         
C             wflux(i,j,k,t) = maskc(i,j,km1)*maskc(i,j,k)*
             wflux(i,j,k,t) = 
     &       wvel(i,j,k,t)*rac(i,j)
     &       *(tracer(i,j,k,t)+tracer(i,j,km1,t))*half         
          enddo
         enddo
        endif
       enddo
      enddo
      
      return
      end

      subroutine advz_cntr_4th(
     &           nx,ny,nr,nt,advectionscheme,rac,
     &           maskc,deltat,dz,wvel,tracer,wflux)


C calculates the centered fourth-order area integrated vertical flux due to advection of a tracer
      implicit none
      include "offline_consts.h"
      include "offline_vert_advection.h"
C Local variable definitions

      integer i,j,k,t,km1,km2,kp1
      real*8 deltat
      real*8 rjm,rj,rjp,rjjm,rjjp
      real*8 rtrans,maskpm,maskbound

      do t=1,nt    
       do k=1,nr
        km2=max(1,k-2)
        km1=max(1,k-1)
        kp1=min(nr,k+1)
        maskpm = 1.d0
        if (k.le.2 .or. k.ge.nr) maskpm = 0.d0
   
        if ( k.eq.1 .or. k.gt.nr) then
         do j=1,ny
          do i=1,nx
           wflux(i,j,k,t) = 0.d0
          enddo
         enddo
        else
         do j=1,ny
          do i=1,nx
           maskbound = maskpm*maskc(i,j,km2)*maskc(i,j,kp1)
           rjp=(tracer(i,j,kp1,t)-tracer(i,j, k ,t))
C     &          *maskc(i,j,kp1)
           rj =(tracer(i,j, k ,t)-tracer(i,j,km1,t))
           rjm=(tracer(i,j,km1,t)-tracer(i,j,km2,t))
C     &          *maskc(i,j,km1)
           rjjp=(rjp-rj)
           rjjm=(rj-rjm)
           
           rtrans=wvel(i,j,k,t)*rac(i,j)

C           wflux(i,j,k,t) = maskc(i,j,km1)*maskc(i,j,k)
C     &       *(
           wflux(i,j,k,t) = 
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

      subroutine advz_cntrfl_2nd(
     &           nx,ny,nr,nt,advectionscheme,rac,
     &           maskc,deltat,dz,wvel,tracer,wflux)


C calculates the centered second-order area integrated meridional flux due to 
C advection of a tracer with a flux limiter
      implicit none
      include "offline_consts.h"
      include "offline_vert_advection.h"
C Local variable definitions
      integer i,j,k,t,km1,km2,kp1
      real*8 deltat
      real*8 recipz(nr)
      real*8 rjm,rj,rjp,wcfl,cr,lcr
      real*8 rtrans

      km2=max(1,k-2)
      km1=max(1,k-1)
      kp1=min(nr,k+1)

      do k=1,nr
       recipz(k)=1/dz(k)
      enddo
      
      do t=1,nt    
       do k=1,nr
        km2=max(1,k-2)
        km1=max(1,k-1)
        kp1=min(nr,k+1)

         do j=1,ny
          do i=1,nx
           rjp=(tracer(i,j,k,t)-tracer(i,j,kp1,t))
C     &         *maskc(i,j,kp1)
           rj =(tracer(i,j,km1,t)-tracer(i,j,k,t))
C     &         *maskc(i,j,k)*maskc(i,j,km1)
           rjm=(tracer(i,j,km2,t)-tracer(i,j,km1,t))
C     &         *maskc(i,j,km1)

           rtrans= wvel(i,j,k,t)*rac(i,j)
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
              
C         wflux(i,j,k,t) = maskc(i,j,km1)*maskc(i,j,k)*(
         wflux(i,j,k,t) = (
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
     &           nx,ny,nr,nt,advectionscheme,rac,
     &           maskc,deltat,dz,wvel,tracer,wflux)

C  calculates the second-order direct space and time (dst-2/lax-wendroff, or simple 
C  first order upwind scheme) area integrated vertical flux due to advection of a tracer 
      implicit none
      include "offline_consts.h"
      include "offline_vert_advection.h"
C Local variable definitions
      integer i,j,k,t,km1
      real*8 deltat
      real*8 wcfl, rlimit, wabs, rtrans
      real*8 recipz(nr)
       
      if ( advectionscheme.eq.20 ) then
       rlimit = 1.d0
      else
       rlimit = 0.d0
      endif  

      do k=1,nr
       recipz(k)=1/dz(k)
      enddo
      
      do t=1,nt    
       do k=1,nr
        
        km1=max(1,k-1)
        
        if ( k.le.1 .or. k.gt.nr) then
         do j=1,ny
          do i=1,nx
           wflux(i,j,k,t) = 0.d0
          enddo
         enddo
        else
         do j=1,ny
          do i=1,nx
          wcfl = abs( wvel(i,j,k,t)*deltat*recipz(k) )
          rtrans= wvel(i,j,k,t)*rac(i,j)

          wabs = abs(rtrans)*rksign
     &       *( 1.d0 - rlimit*(1.d0 - wcfl) )
          wflux(i,j,k,t) = maskc(i,j,km1)*(
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

      subroutine advz_dst_3rd(
     &           nx,ny,nr,nt,advectionscheme,rac,
     &           maskc,deltat,dz,wvel,tracer,wflux)


C  calculates the 3rd-order direct space and time (dst-3) area integrated vertical flux 
C  due to advection of a tracer
      implicit none
      include "offline_consts.h"
      include "offline_vert_advection.h"
C Local variable definitions
      integer i,j,k,t,km1,km2,kp1
      real*8 deltat
      real*8 recipz(nr)
      real*8 rjm,rj,rjp,d0,d1
      real*8 rtrans,wcfl

      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx
          wflux(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do k=1,nr
       recipz(k)=1/dz(k)
      enddo
      
      do t=1,nt    
       do k=1,nr
        km2=max(1,k-2)
        km1=max(1,k-1)
        kp1=min(nr,k+1)

         do j=1,ny
          do i=1,nx
           rjp=(tracer(i,j,k,t)-tracer(i,j,kp1,t))
     &         *maskc(i,j,kp1)*maskc(i,j,k)
           rj =(tracer(i,j,km1,t)-tracer(i,j,k,t))
     &         *maskc(i,j,k)*maskc(i,j,km1)
           rjm=(tracer(i,j,km2,t)-tracer(i,j,km1,t))
     &         *maskc(i,j,km1)*maskc(i,j,km2)

           rtrans= wvel(i,j,k,t)*rac(i,j)
           wcfl=abs( wvel(i,j,k,t)*deltat*recipz(k) )

           d0=(2.d0-wcfl)*(1.d0-wcfl)*onesixth
           d1=(1.d0-wcfl*wcfl)*onesixth

           wflux(i,j,k,t)=maskc(i,j,k)*maskc(i,j,km1)           
C           wflux(i,j,k,t)=(
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

      subroutine advz_dstfl_3rd(
     &           nx,ny,nr,nt,advectionscheme,rac,
     &           maskc,deltat,dz,wvel,tracer,wflux)


C  calculates the 3rd-order direct space and time (dst-3) area integrated zonal flux 
C  due to advection of a tracer with flux limiting.
      implicit none
      include "offline_consts.h"
      include "offline_vert_advection.h"
C Local variable definitions
      integer i,j,k,t,km1,km2,kp1
      real*8 deltat
      real*8 recipz(nr)
      real*8 rjm,rj,rjp,d0,d1,psip,psim
      real*8 thetap,thetam
      real*8 rtrans,wcfl

      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx
          wflux(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo        

      do k=1,nr
       recipz(k)=1/dz(k)
      enddo
      
      do t=1,nt    
       do k=1,nr
        km2=max(1,k-2)
        km1=max(1,k-1)
        kp1=min(nr,k+1)

         do j=1,ny
          do i=1,nx
           rjp=(tracer(i,j, k ,t)-tracer(i,j,kp1,t))
     &         *maskc(i,j,kp1)*maskc(i,j, k )
           rj =(tracer(i,j,km1,t)-tracer(i,j, k ,t))
     &         *maskc(i,j, k )*maskc(i,j,km1)
           rjm=(tracer(i,j,km2,t)-tracer(i,j,km1,t))
     &         *maskc(i,j,km1)*maskc(i,j,km2)

           wcfl = abs( wvel(i,j,k,t)*deltat*recipz(k) )
           rtrans= wvel(i,j,k,t)*rac(i,j)
           
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

           wflux(i,j,k,t)=maskc(i,j,k)*maskc(i,j,km1)      
C           wflux(i,j,k,t)=(
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

      subroutine advz_mpl_7th(
     &           nx,ny,nr,nt,advectionscheme,rac,
     &           maskc,deltat,dz,wvel,tracer,wflux)

     
C calculates the 7th order dst area integrated zonal flux due to 
C advection of a tracer with a monotone preserving limiter
      implicit none
      include "offline_consts.h"
      include "offline_vert_advection.h"
C Local variable definitions
C      integer nx,ny,nr,nt
      integer i,j,k,t
      integer kp3,kp2,kp1,km1,km2,km3,km4
      real*8 deltat
      real*8 recipz(nr)
      real*8 cfl,psi,wtrans
      real*8 fac,delip,deli,phi,rp1h,rp1h_cfl
      real*8 recip_delip, recip_deli
      real*8 qippp,qipp,qip,qi,qim,qimm,qimmm
      real*8 mskipp,mskip,mski,mskim,mskimm,mskimmm
      real*8 d2,d2p1,d2m1,a,b,c,d
      real*8 dp1h,dm1h, phimd,philc,phimin,phnx
      real*8 delm,delp,delmm,delpp,delmmm,delppp
      real*8 del2mm,del2m,del2,del2p,del2pp
      real*8 del3mm,del3m,del3p,del3pp
      real*8 del4m,del4,del4p
      real*8 del5m,del5p
      real*8 del6
      
      do k=1,nr
       recipz(k)=1/dz(k)
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

         do j=1,ny
          do i=1,nx
           wtrans= wvel(i,j,k,t)*rac(i,j)
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
         phnx = min(max(2.d0/(1.d0-cfl),phimd),
     &        max(0.d0,2.d0*rp1h_cfl,philc))
         phi = max(phimin,min(phi,phnx))

         psi = phi * 0.5d0 * (1.d0 - cfl)
         wflux(i,j,k,t) = wtrans*( qi + psi*delip )
        else
         wflux(i,j,k,t) = 0.d0
        endif

       enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine advz_uw_3rd(
     &           nx,ny,nr,nt,advectionscheme,rac,
     &           maskc,deltat,dz,wvel,tracer,wflux)

     
C calculates the third-order upwind biased area integrated vertical flux due to 
C advection of a tracer
      implicit none
      include "offline_consts.h"
      include "offline_vert_advection.h"
C Local variable definitions
      integer i,j,k,t,km1,km2,kp1
      real*8 deltat
      real*8 rjm,rj,rjp,rjjp,rjjm
      real*8 rtrans

      do t=1,nt    
       do k=1,nr
        do j=1,ny
         do i=1,nx
          wflux(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo             
      
      do t=1,nt    
       do k=1,nr
        km2=max(1,k-2)
        km1=max(1,k-1)
        kp1=min(nr,k+1)

         do j=1,ny
          do i=1,nx
           rjp=(tracer(i,j,kp1,t)-tracer(i,j, k ,t))
     &         *maskc(i,j,kp1)*maskc(i,j,k)
           rj =(tracer(i,j, k ,t)-tracer(i,j,km1,t))
     &         *maskc(i,j,k)*maskc(i,j,km1)
           rjm=(tracer(i,j,km1,t)-tracer(i,j,km2,t))
     &         *maskc(i,j,km1)*maskc(i,j,km2)
           
           rjjp = rjp-rj
           rjjm = rj-rjm

           rtrans= wvel(i,j,k,t)*rac(i,j)
         
           wflux(i,j,k,t) = maskc(i,j,k)
C     &                  *maskc(i,j,km1)
C            wflux(i,j,k,t) =
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

      subroutine advz_uwcntr_2nd(
     &           nx,ny,nr,nt,advectionscheme,rac,
     &           maskc,deltat,dz,wvel,tracer,wflux)

C calculates the center second upwind biased area integrated vertical flux due to 
C advection of a tracer
      implicit none
      include "offline_consts.h"
      include "offline_vert_advection.h"
C Local variable definitions
      integer i,j,k,t,km1,kp1
      real*8 deltat
              
      do t=1,nt    
       do k=1,nr
        km1=max(1,k-1)
        kp1=min(k+1,nr)
         do j=1,ny
          do i=1,nx
            if ( k .eq. 1 ) then
              wflux(i,j,1,t) = 0.d0
            endif
            
            if ( wvel(i,j,kp1,t) .ge. 0.d0 ) then
C              wflux(i,j,kp1,t) = maskc(i,j,km1)*maskc(i,j,kp1)*
              wflux(i,j,kp1,t) = 
     &                        wvel(i,j,kp1,t)*tracer(i,j,kp1,t)*rac(i,j)
            else
C              wflux(i,j,kp1,t) = maskc(i,j,km1)*maskc(i,j,kp1)*
              wflux(i,j,kp1,t) = 
     &                        wvel(i,j,kp1,t)*tracer(i,j, k ,t)*rac(i,j)
            endif
          enddo
         enddo
       enddo
      enddo
      
      return
      end
      
      subroutine calc_flux_divergence(nx,ny,nr,nt,klev,dz,volc,
     &                       maskc,fzon,fmer,fvert,hfint,rfint)
      
      implicit none
      include "offline_consts.h"
      include "offline_divergence.h"
C Local variable definitions
      integer i,j,k,t
      integer ip1,jp1,km1,kp1
      real*8 recip_volc(nx,ny,nr)
      real*8 hdiv (nx,ny,nr)
      real*8 rdiv (nx,ny,nr)
C
C#ifdef DEBUG_MESSAGES 
C      integer, external :: mexPrintf
C      integer msg
C      character(len=250) :: line
C#endif
C
CC From config_summary
CC rkSign =   /* index orientation relative to vertical coordinate */
CC -1.000000000000000E+00
C      real*8 rkSign
C      PARAMETER(rkSign=-1.D0)
C      
CC     Initialize fluxes to zero and calculate reciprocals
C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt('zero the outputs')    
C      write(line,*)'dimensions:\ni:',nx,'\nj:',ny,
C     &                        '\nk',nr,'\nt:',nt 
C      msg=mexPrintf(line//achar(13)) 
C#endif
      do t=1,nt
       do j=1,ny
        do i=1,nx  
         hfint(i,j,t)=0.D0
         rfint(i,j,t)=0.D0
        enddo
       enddo  
      enddo

C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt('calculate 1/volume')    
C#endif
      do j=1,ny
       do i=1,nx
        do k=1,nr 
         recip_volc(i,j,k)=1.D0/volc(i,j,k)
        enddo
       enddo
      enddo

CC--   Divergence of fluxes
C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt('calculate divergence')    
C#endif
      do t=1,nt
       do k=1,nr
        do j=1,ny
         do i=1,nx
          if (i .eq. nx) THEN
           ip1=1
          else
           ip1=i+1
          endif
          
          jp1=MIN(j+1,ny)
          
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

CC     Integrate vertically 
C#ifdef DEBUG_MESSAGES 
C      if (t .eq. 1) then
C         call mexWarnMsgTxt('integrate vertically')  
C      endif
C#endif  
       do j=1,ny
        do i=1,nx  
         do k=1,INT(klev(i,j,t))    
          hfint(i,j,t)=hfint(i,j,t)
     &      +(hdiv(i,j,k)*dz(k)*maskc(i,j,k))
     
          rfint(i,j,t)=rfint(i,j,t)
     &      +(rdiv(i,j,k)*dz(k)*maskc(i,j,k))
         enddo
        enddo
       enddo
      enddo
            
      return
      end
    
      subroutine calc_flux_divergence_smooth(nx,ny,nr,nt,klev,dz,
     &                       volc,maskc,fzon,fmer,fvert,hfint,rfint)
      
      implicit none
      include "offline_consts.h"
      include "offline_divergence.h"
C Local variable definitions
      integer i,j,k,t  
      integer im2,im1,ip1,ip2,ip3
      integer jm2,jm1,jp1,jp2,jp3
      integer km1,kp1
      real*8 recip_volc(nx,ny,nr)
      real*8 hdiv (nx,ny,nr)
      real*8 rdiv (nx,ny,nr)
      real*8 fzon1,fmer2,fzon3,fmer4
      real*8 fverp0, fverp1
C From config_summary
C rkSign =   /* index orientation relative to vertical coordinate */
C -1.000000000000000E+00
C      real*8 rkSign
C      PARAMETER(rkSign=-1.D0)
C            
C#ifdef DEBUG_MESSAGES 
C      integer, external :: mexPrintf
C      integer msg
C      character(len=250) :: line
Cendif
     
C     Initialize fluxes to zero and calculate reciprocals
C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt('zero the outputs')    
C      write(line,*)'dimensions:\ni:',nx,'\nj:',ny,
C     &                        '\nk',nr,'\nt:',nt 
C      msg=mexPrintf(line//achar(13)) 
C#endif
      do t=1,nt
       do j=1,ny
        do i=1,nx  
         hfint(i,j,t)=0.D0
         rfint(i,j,t)=0.D0
        enddo
       enddo  
      enddo

C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt('calculate 1/volume')    
C#endif
      do j=1,ny
       do i=1,nx
        do k=1,nr 
         recip_volc(i,j,k)=1.D0/volc(i,j,k)
        enddo
       enddo
      enddo

C--   Divergence of fluxes
C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt('calculate divergence')    
C#endif
      do t=1,nt    
       do j=1,ny
        do i=1,nx
         do k=1,nr 
          hdiv(i,j,k)=0.D0
          rdiv(i,j,k)=0.D0
         enddo
        enddo
       enddo
       
       do k=1,nr
        do j=1,ny
         do i=1,nx
           if (maskc(i,j,k).gt.0.D0) then

          if (i .eq. nx) then
           ip1=1
           ip2=2
           ip3=3
          elseif (i .eq. nx-1) then
           ip1=nx
           ip2=1
           ip3=2
          elseif (i .eq. nx-2) then
           ip1=nx-1
           ip2=nx
           ip3=1
          else
           ip1=i+1
           ip2=i+2
           ip3=i+3
          endif
          
          if (i .eq. 1) then
           im1=nx
           im2=nx-1
          elseif (i .eq. 2) then
           im1=1
           im2=nx
          else
           im1=i-1
           im2=i-2
          endif
          
          jm2=MAX(j-2,1   )
          jm1=MAX(j-1,1   )
          jp1=MIN(j+1,ny)
          jp2=MIN(j+2,ny)
          jp3=MIN(j+3,ny)
          
          km1  = MAX(1 ,k-1)
          kp1  = MIN(nr,k+1)
        
C Algorithm:  The addresses for position refer to the following 9-point
C stencil (after HW80, Fig. 5.10):
C
C    6 --- 2 --- 5
C    |     |     |
C    3 --- 0 --- 1
C    |     |     |
C    7 --- 4 --- 8
C
C where the point "0" is the point divergence is being calculated
          fzon1=( fzon(im1,j  ,k,t)*maskc(im1,j  ,k)
     &          + fzon(i  ,jm1,k,t)*maskc(i  ,jm1,k) 
     &          + fzon(i  ,j  ,k,t)*maskc(i  ,j  ,k)
     &          + fzon(i  ,jp1,k,t)*maskc(i  ,jp1,k)
     &          + fzon(ip1,jm2,k,t)*maskc(ip1,jm2,k)
     &          + fzon(ip1,jm1,k,t)*maskc(ip1,jm1,k)
     &          + fzon(ip1,j  ,k,t)*maskc(ip1,j  ,k)
     &          + fzon(ip1,jp1,k,t)*maskc(ip1,jp1,k)
     &          + fzon(ip1,jp2,k,t)*maskc(ip1,jp2,k)
     &          + fzon(ip2,jm1,k,t)*maskc(ip2,jm1,k)
     &          + fzon(ip2,j  ,k,t)*maskc(ip2,j  ,k)
     &          + fzon(ip2,jp1,k,t)*maskc(ip2,jp1,k)
     &          + fzon(ip3,j  ,k,t)*maskc(ip3,j  ,k)
     &          )/(
     &          + maskc(im1,j  ,k)
     &          + maskc(i  ,jm1,k)
     &          + maskc(i  ,j  ,k)
     &          + maskc(i  ,jp1,k)
     &          + maskc(ip1,jm2,k)
     &          + maskc(ip1,jm1,k)
     &          + maskc(ip1,j  ,k)
     &          + maskc(ip1,jp1,k)
     &          + maskc(ip1,jp2,k)
     &          + maskc(ip2,jm1,k)
     &          + maskc(ip2,j  ,k)
     &          + maskc(ip2,jp1,k)
     &          + maskc(ip3,j  ,k)
     &          )

          fzon3=( fzon(im2,j  ,k,t)*maskc(im2,j  ,k)
     &          + fzon(im1,jm1,k,t)*maskc(im1,jm1,k)
     &          + fzon(im1,j  ,k,t)*maskc(im1,j  ,k)
     &          + fzon(im1,jp1,k,t)*maskc(im1,jp1,k)
     &          + fzon(i  ,jm2,k,t)*maskc(i  ,jm2,k)
     &          + fzon(i  ,jm1,k,t)*maskc(i  ,jm1,k)
     &          + fzon(i  ,j  ,k,t)*maskc(i  ,j  ,k)
     &          + fzon(i  ,jp1,k,t)*maskc(i  ,jp1,k)
     &          + fzon(i  ,jp2,k,t)*maskc(i  ,jp2,k)
     &          + fzon(ip1,jm1,k,t)*maskc(ip1,jm1,k)
     &          + fzon(ip1,j  ,k,t)*maskc(ip1,j  ,k)
     &          + fzon(ip1,jp1,k,t)*maskc(ip1,jp1,k)
     &          + fzon(ip2,j  ,k,t)*maskc(ip2,j  ,k)
     &          )/(
     &          + maskc(im2,j  ,k)
     &          + maskc(im1,jm1,k)
     &          + maskc(im1,j  ,k)
     &          + maskc(im1,jp1,k)
     &          + maskc(i  ,jm2,k)
     &          + maskc(i  ,jm1,k)
     &          + maskc(i  ,j  ,k)
     &          + maskc(i  ,jp1,k)
     &          + maskc(i  ,jp2,k)
     &          + maskc(ip1,jm1,k)
     &          + maskc(ip1,j  ,k)
     &          + maskc(ip1,jp1,k)
     &          + maskc(ip2,j  ,k)
     &          )

          fmer2=( fmer(i  ,jm1,k,t)*maskc(i  ,jm1,k)
     &          + fmer(im1,j  ,k,t)*maskc(im1,j  ,k)
     &          + fmer(i  ,j  ,k,t)*maskc(i  ,j  ,k)
     &          + fmer(ip1,j  ,k,t)*maskc(ip1,j  ,k)
     &          + fmer(im2,jp1,k,t)*maskc(im2,jp1,k)
     &          + fmer(im1,jp1,k,t)*maskc(im1,jp1,k)
     &          + fmer(i  ,jp1,k,t)*maskc(i  ,jp1,k)
     &          + fmer(ip1,jp1,k,t)*maskc(ip1,jp1,k)
     &          + fmer(ip2,jp1,k,t)*maskc(ip2,jp1,k)
     &          + fmer(im1,jp2,k,t)*maskc(im1,jp2,k)
     &          + fmer(i  ,jp2,k,t)*maskc(i  ,jp2,k)
     &          + fmer(ip1,jp2,k,t)*maskc(ip1,jp2,k)
     &          + fmer(i  ,jp3,k,t)*maskc(i  ,jp3,k)
     &          )/(
     &          + maskc(i  ,jm1,k)
     &          + maskc(im1,j  ,k)
     &          + maskc(i  ,j  ,k)
     &          + maskc(ip1,j  ,k)
     &          + maskc(im2,jp1,k)
     &          + maskc(im1,jp1,k)
     &          + maskc(i  ,jp1,k)
     &          + maskc(ip1,jp1,k)
     &          + maskc(ip2,jp1,k)
     &          + maskc(im1,jp2,k)
     &          + maskc(i  ,jp2,k)
     &          + maskc(ip1,jp2,k)
     &          + maskc(i  ,jp3,k)
     &          )

          fmer4=( fmer(i  ,jm2,k,t)*maskc(i  ,jm2,k)
     &          + fmer(im1,jm1,k,t)*maskc(im1,jm1,k)
     &          + fmer(i  ,jm1,k,t)*maskc(i  ,jm1,k)
     &          + fmer(ip1,jm1,k,t)*maskc(ip1,jm1,k)
     &          + fmer(im2,j  ,k,t)*maskc(im2,j  ,k)
     &          + fmer(im1,j  ,k,t)*maskc(im1,j  ,k)
     &          + fmer(i  ,j  ,k,t)*maskc(i  ,j  ,k)
     &          + fmer(ip1,j  ,k,t)*maskc(ip1,j  ,k)
     &          + fmer(ip2,j  ,k,t)*maskc(ip2,j  ,k)
     &          + fmer(im1,jp1,k,t)*maskc(im1,jp1,k)
     &          + fmer(i  ,jp1,k,t)*maskc(i  ,jp1,k)
     &          + fmer(ip1,jp1,k,t)*maskc(ip1,jp1,k)
     &          + fmer(i  ,jp2,k,t)*maskc(i  ,jp2,k)
     &          )/(
     &          + maskc(i  ,jm2,k)
     &          + maskc(im1,jm1,k)
     &          + maskc(i  ,jm1,k)
     &          + maskc(ip1,jm1,k)
     &          + maskc(im2,j  ,k)
     &          + maskc(im1,j  ,k)
     &          + maskc(i  ,j  ,k)
     &          + maskc(ip1,j  ,k)
     &          + maskc(ip2,j  ,k)
     &          + maskc(im1,jp1,k)
     &          + maskc(i  ,jp1,k)
     &          + maskc(ip1,jp1,k)
     &          + maskc(i  ,jp2,k)
     &          )

         fverp0=( fvert(im2,j  ,k,t)*maskc(im2,j  ,k)
     &          + fvert(im1,jm1,k,t)*maskc(im1,jm1,k)
     &          + fvert(im1,j  ,k,t)*maskc(im1,j  ,k)
     &          + fvert(im1,jp1,k,t)*maskc(im1,jp1,k)
     &          + fvert(i  ,jm2,k,t)*maskc(i  ,jm2,k)
     &          + fvert(i  ,jm1,k,t)*maskc(i  ,jm1,k)
     &          + fvert(i  ,j  ,k,t)*maskc(i  ,j  ,k)
     &          + fvert(i  ,jp1,k,t)*maskc(i  ,jp1,k)
     &          + fvert(i  ,jp2,k,t)*maskc(i  ,jp2,k)
     &          + fvert(ip1,jm1,k,t)*maskc(ip1,jm1,k)
     &          + fvert(ip1,j  ,k,t)*maskc(ip1,j  ,k)
     &          + fvert(ip1,jp1,k,t)*maskc(ip1,jp1,k)
     &          + fvert(ip2,j  ,k,t)*maskc(ip2,j  ,k)
     &          )/( 
     &          + maskc(im2,j  ,k)
     &          + maskc(im1,jm1,k)
     &          + maskc(im1,j  ,k)
     &          + maskc(im1,jp1,k)
     &          + maskc(i  ,jm2,k)
     &          + maskc(i  ,jm1,k)
     &          + maskc(i  ,j  ,k)
     &          + maskc(i  ,jp1,k)
     &          + maskc(i  ,jp2,k)
     &          + maskc(ip1,jm1,k)
     &          + maskc(ip1,j  ,k)
     &          + maskc(ip1,jp1,k)
     &          + maskc(ip2,j  ,k)
     &          )

         fverp1=( fvert(im2,j  ,kp1,t)*maskc(im2,j  ,kp1)
     &          + fvert(im1,jm1,kp1,t)*maskc(im1,jm1,kp1)
     &          + fvert(im1,j  ,kp1,t)*maskc(im1,j  ,kp1)
     &          + fvert(im1,jp1,kp1,t)*maskc(im1,jp1,kp1)
     &          + fvert(i  ,jm2,kp1,t)*maskc(i  ,jm2,kp1)
     &          + fvert(i  ,jm1,kp1,t)*maskc(i  ,jm1,kp1)
     &          + fvert(i  ,j  ,kp1,t)*maskc(i  ,j  ,kp1)
     &          + fvert(i  ,jp1,kp1,t)*maskc(i  ,jp1,kp1)
     &          + fvert(i  ,jp2,kp1,t)*maskc(i  ,jp2,kp1)
     &          + fvert(ip1,jm1,kp1,t)*maskc(ip1,jm1,kp1)
     &          + fvert(ip1,j  ,kp1,t)*maskc(ip1,j  ,kp1)
     &          + fvert(ip1,jp1,kp1,t)*maskc(ip1,jp1,kp1)
     &          + fvert(ip2,j  ,kp1,t)*maskc(ip2,j  ,kp1)
     &          )/(
     &          + maskc(im2,j  ,kp1)
     &          + maskc(im1,jm1,kp1)
     &          + maskc(im1,j  ,kp1)
     &          + maskc(im1,jp1,kp1)
     &          + maskc(i  ,jm2,kp1)
     &          + maskc(i  ,jm1,kp1)
     &          + maskc(i  ,j  ,kp1)
     &          + maskc(i  ,jp1,kp1)
     &          + maskc(i  ,jp2,kp1)
     &          + maskc(ip1,jm1,kp1)
     &          + maskc(ip1,j  ,kp1)
     &          + maskc(ip1,jp1,kp1)
     &          + maskc(ip2,j  ,kp1)
     &          )

           hdiv(i,j,k) = -recip_volc(i,j,k)
     &          *((fzon1-fzon3)+(fmer2-fmer4))
     &          *maskc(i,j,k)*maskc(ip1,j,k)*maskc(i,jp1,k)

          rdiv(i,j,k) = -recip_volc(i,j,k)
     &          *(fverp1-fverp0)*rkSign
     &          *maskc(i,j,k)*maskc(i,j,kp1)
     
         endif
         
         enddo
        enddo
       enddo

C     Integrate vertically          
C#ifdef DEBUG_MESSAGES 
C      if (t .eq. 1) then
C         call mexWarnMsgTxt('integrate vertically')  
C      endif
C#endif  
       do j=1,ny
        do i=1,nx  
         do k=1,INT(klev(i,j,t))    
          hfint(i,j,t)=hfint(i,j,t)
     &      +(hdiv(i,j,k)*dz(k)*maskc(i,j,k))
     
          rfint(i,j,t)=rfint(i,j,t)
     &      +(rdiv(i,j,k)*dz(k)*maskc(i,j,k))
         enddo
        enddo
       enddo
      enddo
            
      return
      end


CCC----------------DIFFUSION IN X DIRECTION----------------
C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt(
C     & 'Calculating explicit diffusion in the X direction') 
C#endif  
C         call diffx(dims(1),dims(2),dims(3),dims(4),
C     &           mskc,dz,dx,rax,kux,kuz,tracer,uflux)
C
CCC----------------DIFFUSION IN Y DIRECTION----------------
C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt(
C     & 'Calculating explicit diffusion in the Y direction') 
C#endif  
C         call diffy(dims(1),dims(2),dims(3),dims(4),
C     &           mskc,dz,dy,ray,kvy,kvz,tracer,vflux)
C
CCC----------------DIFFUSION IN Z DIRECTION----------------
C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt(
C     & 'Calculating explicit diffusion in the Z direction') 
C#endif  
C         call diffz(dims(1),dims(2),dims(3),dims(4),
C     &           msku,mskv,mskc,dx,dy,rac,kwx,
C     &           kwy,tracer,wflux)
C     
CC------------------------COMPUTATIONAL SUBROUTINES------------------------

      subroutine diffx(nx,ny,nr,nt,
     &     maskc,dz,dx,rax,kux,kuz,tracer,udiff )

      implicit none
      include "offline_consts.h"
C Local variable definitions
      integer i, j, k, t, kp1, km1, im1
      real*8, intent(in), dimension(nr)              :: dz    
      real*8, intent(in), dimension(nx, ny)          :: dx
      real*8, intent(in), dimension(nx, ny, nr)      :: rax
      real*8, intent(in), dimension(nx, ny, nr)      :: maskc
      real*8, intent(in), dimension(nx, ny, nr, nt)  :: kux
      real*8, intent(in), dimension(nx, ny, nr, nt)  :: kuz
      real*8, intent(in), dimension(nx, ny, nr, nt)  :: tracer
      real*8, intent(out),dimension(nx, ny, nr, nt)  :: udiff
      real*8 recip_dz(nr)
      real*8 recip_dx(nx,ny)
      real*8 dtdz

C     initialise flux to zero
      do t=1,nt
       do k=1,nr
         recip_dz(k)=1.d0/dz(k)  
        do j=1,ny
         do i=1,nx
          udiff(i,j,k,t) = 0.d0
          
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
        do j=1,ny
         do i=1,nx
          if (i .eq. 1) then
           im1=nx
          else
           im1=i-1
          endif
          udiff(i,j,k,t) = udiff(i,j,k,t)
     &      -rax(i,j,k)
     &      *kux(i,j,k,t)
     &      *recip_dx(i,j)
     &      *(tracer(i,j,k,t)-tracer(im1,j,k,t))
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
        do j=1,ny
         do i=1,nx
           if (i .eq. 1) then
            im1=nx
           else
            im1=i-1
           endif
             
           dtdz =  half*(
     &      half*recip_dz(k)*
     &        (maskc(im1,j,k)*(tracer(im1,j,km1,t)-tracer(im1,j,k,t))
     &        +maskc( i ,j,k)*(tracer( i ,j,km1,t)-tracer( i ,j,k,t))
     &        )
     &      +half*recip_dz(kp1)*
     &        (maskc(im1,j,kp1)*(tracer(im1,j,k,t)-tracer(im1,j,kp1,t))
     &        +maskc( i ,j,kp1)*(tracer( i ,j,k,t)-tracer( i ,j,kp1,t))
     &        )
     &      )
  
           udiff(i,j,k,t) = udiff(i,j,k,t) 
     &                    - rax(i,j,k)*kuz(i,j,k,t)*dtdz
         enddo
        enddo
       enddo
      enddo
       
      return
      end

      subroutine diffy(nx, ny, nr, nt,
     &     maskc, dz, dy, ray, kvy, kvz,tracer,vdiff )

      implicit none
      include "offline_consts.h"
C Local variable definitions
      integer i, j, k, t, kp1, km1, jm1
      real*8, intent(in), dimension(nr)              :: dz    
      real*8, intent(in), dimension(nx, ny)          :: dy
      real*8, intent(in), dimension(nx, ny, nr)      :: ray
      real*8, intent(in), dimension(nx, ny, nr)      :: maskc
      real*8, intent(in), dimension(nx, ny, nr, nt)  :: kvy
      real*8, intent(in), dimension(nx, ny, nr, nt)  :: kvz
      real*8, intent(in), dimension(nx, ny, nr, nt)  :: tracer
      real*8, intent(out),dimension(nx, ny, nr, nt)  :: vdiff
      real*8 recip_dz(nr)
      real*8 recip_dy(nx,ny)
      real*8 dtdz

C     initialise flux to zero
      do t=1,nt
       do k=1,nr
         recip_dz(k)=1.d0/dz(k)  
        do j=1,ny
         do i=1,nx
          vdiff(i,j,k,t) = 0.d0
       
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
        do j=1,ny
         do i=1,nx
          jm1=max(j-1, 1 )
          vdiff(i,j,k,t) = vdiff(i,j,k,t)
     &      -ray(i,j,k)
     &      *kvy(i,j,k,t)
     &      *recip_dy(i,j)
     &      *(tracer(i,j,k,t)-tracer(i,jm1,k,t))
         enddo
        enddo
       enddo
      enddo
      
      do t=1,nt
       do k=1,nr
       km1 = max(k-1,1)
       kp1 = min(k+1,nr)
C    vertical gradients interpolated to u points
        do j=1,ny
         do i=1,nx

          jm1=max(j-1, 1 )

C    vertical gradients interpolated to v points
           dtdz =  half*(
     &     half*recip_dz(k)*
     &        (maskc(i,jm1,k)*(tracer(i,jm1,km1,t)-tracer(i,jm1,k,t))
     &        +maskc(i, j ,k)*(tracer(i, j ,km1,t)-tracer(i, j ,k,t))
     &        )
     &     +half*recip_dz(kp1)*
     &        (maskc(i,jm1,kp1)*(tracer(i,jm1,k,t)-tracer(i,jm1,kp1,t))
     &        +maskc(i, j ,kp1)*(tracer(i, j ,k,t)-tracer(i, j ,kp1,t))
     &        )
     &     )

           vdiff(i,j,k,t) = vdiff(i,j,k,t) 
     &                    - ray(i,j,k)*kvz(i,j,k,t)*dtdz
         enddo
        enddo
       enddo
      enddo
      
      return
      end

      subroutine diffz(nx, ny, nr, nt,
     &     masku, maskv, maskc,dx, dy, rac, kwx, kwy,
     &     tracer, wediff )

      implicit none
      include "offline_consts.h"
C Local variable definitions
      integer i, j, k, t
      integer kp1, km1, jm1, im1, ip1, jp1
      real*8 recip_dy(nx,ny)
      real*8 recip_dx(nx,ny)
      real*8 tmp1, tmp2, dtdx, dtdy
      real*8, intent(in), dimension(nx, ny)          :: dx    
      real*8, intent(in), dimension(nx, ny)          :: dy
      real*8, intent(in), dimension(nx, ny)          :: rac
      real*8, intent(in), dimension(nx, ny, nr)      :: masku
      real*8, intent(in), dimension(nx, ny, nr)      :: maskv
      real*8, intent(in), dimension(nx, ny, nr)      :: maskc
      real*8, intent(in), dimension(nx, ny, nr, nt)  :: kwx
      real*8, intent(in), dimension(nx, ny, nr, nt)  :: kwy
      real*8, intent(in), dimension(nx, ny, nr, nt)  :: tracer
      real*8, intent(out),dimension(nx, ny, nr, nt)  :: wediff

C#ifdef DEBUG_MESSAGES 
C      integer, external :: mexPrintf
C      integer msg
C      character(len=250) :: line
C#endif

C     initialise flux to zero
C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt('DIFFZ: set flux to zero') 
C#endif  
      do t=1,nt
       do k=1,nr
        do j=1,ny
         do i=1,nx
          wediff(i,j,k,t) = 0.d0
         enddo
        enddo
       enddo
      enddo

C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt('DIFFZ: calculate 1/dx and 1/dy') 
C#endif  
      
      do t=1,nt
       do k=1,nr
        do j=1,ny
         do i=1,nx
          if (maskc(i,j,1) .ne. 0.d0) then
           recip_dx(i,j)=(1.d0/dx(i,j))
           recip_dy(i,j)=(1.d0/dy(i,j))
          else
           recip_dx(i,j)=0.d0
           recip_dy(i,j)=0.d0
          endif
C          if (mxIsNaN(recip_dx(i,j)) .eq. 1) then
C          catch a NAN?
C           call mexErrMsgTxt('ERROR: NAN in recip_dx during DIFFZ.')
C          endif
C          if (mxIsNaN(recip_dy(i,j)) .eq. 1) then
C          catch a NAN?
C           call mexErrMsgTxt('ERROR: NAN in recip_dy during DIFFZ.')
C          endif
         enddo
        enddo
       enddo
      enddo

C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt('DIFFZ: calculate gradients and flux') 
C#endif  
C    horizontal gradients interpolated to w points
      do t=1,nt
       do k=2,nr
        do j=1,ny
         do i=1,nx
       
         if (i .eq. 1) then
          im1=nx
         else
          im1=i-1
         endif
          
         if (i .eq. nx) then
          ip1=1
         else
          ip1=i+1
         endif

         km1 = max(k-1,1)
         kp1 = min(k+1,nr)
              
         if ( tracer( i ,j,k,t) .ne. 0 .AND. 
     &        tracer(im1,j,k,t) .ne. 0) then
            tmp1=half*(masku(ip1,j,k)
     &            *recip_dx(ip1,j)*(tracer(ip1,j,k,t)-tracer( i ,j,k,t))
     &           +masku(i,j,k)
     &            *recip_dx( i ,j)*(tracer( i ,j,k,t)-tracer(im1,j,k,t))
     &            )
             
            tmp2=half*(masku(ip1,j,km1)
     &        *recip_dx(ip1,j)*(tracer(ip1,j,km1,t)-tracer( i ,j,km1,t))
     &           +masku(i,j,km1)
     &        *recip_dx( i ,j)*(tracer( i ,j,km1,t)-tracer(im1,j,km1,t)) 
     &            )     
   
            dtdx = half*(tmp1+tmp2)
         else
            dtdx = 0.d0
         endif
                
         jm1=max(j-1, 1 )
         jp1=min(j+1,ny)

         if ( tracer(i, j ,k,t) .ne. 0 .AND. 
     &        tracer(i,jp1,k,t) .ne. 0) then

            tmp1=half*(maskv(i,j,k)
     &            *recip_dy(i,j)*(tracer(i,j,k,t)-tracer(i,jm1,k,t))
     &           +maskv(i,j+1,k)
     &            *recip_dy(i,jp1)*(tracer(i,jp1,k,t)-tracer(i,j,k,t))
     &            )

            tmp2=half*(maskv(i,j,km1)
     &         *recip_dy(i,j)*(tracer(i,j,km1,t)-tracer(i,jm1,km1,t))
     &           +maskv(i,jp1,km1)
     &         *recip_dy(i,jp1)*(tracer(i,jp1,km1,t)-tracer(i,j,km1,t))
     &            )
          
            dtdy = half*(tmp1+tmp2)
         else
            dtdy = 0.d0
         endif
         
C    off-diagonal components of vertical flux
         wediff(i,j,k,t) = wediff(i,j,k,t)
     &       -rac(i,j)*maskc(i,j,k)
     &       *( kwx(i,j,k,t)*dtdx
     &        + kwy(i,j,k,t)*dtdy )
         enddo
        enddo
       enddo
      enddo 

      return
      end

       subroutine implicit_mixing(nx, ny, nr, nt,
     &     deltat,niters,maskc,rac,dz,dzc,diffkz,
     &     tracer,widiff)

c     solve implicitly vertical advection and diffusion
      implicit none
      include "offline_consts.h"
C Local variable definitions
      integer i,j,k,t,iter
      real*8 deltat
      real*8 niters
      real*8, intent(in), dimension(nr)              :: dz    
      real*8, intent(in), dimension(nr)              :: dzc    
      real*8, intent(in), dimension(nx, ny)          :: rac  
      real*8, intent(in), dimension(nx, ny, nr)      :: maskc
      real*8, intent(in), dimension(nx, ny, nr, nt)  :: diffkz        
      real*8, intent(in), dimension(nx, ny, nr, nt)  :: tracer        
      real*8, intent(out), dimension(nx, ny, nr, nt) :: widiff
      real*8 recip_drf(nr)
      real*8 recip_drc(nr)
      real*8 tmptracer   (nx,ny,nr)
      real*8 itertracer  (nx,ny,nr)
      real*8 b5d      (nx,ny,nr)
      real*8 c5d      (nx,ny,nr)
      real*8 d5d      (nx,ny,nr)    

      do t=1,nt
C     initialise
       do k=1,nr        
        recip_drf(k)=1.d0/dz(k)
        recip_drc(k)=1.d0/dzc(k)
        do j=1,ny
         do i=1,nx
          b5d(i,j,k) = 0.d0
          c5d(i,j,k) = 1.d0
          d5d(i,j,k) = 0.d0
         enddo
        enddo
       enddo

C     set the tri-diagonal matrix to solve the implicit diffusion problem
C      1rst lower diagonal :
        do k=2,nr
         do j=1,ny
          do i=1,nx
            b5d(i,j,k) =-deltat
     &                  *maskc(i,j,k-1)*maskc (i,j,k)  
     &                  *recip_drf(k)
     &                  *diffkz(i,j,k,t)*recip_drc(k)
          enddo
         enddo
        enddo
C       1rst upper diagonal :
        do k=1,nr-1
         do j=1,ny
          do i=1,nx
            d5d(i,j,k) = -deltat
     &                 *maskc(i,j,k+1)*maskc (i,j,k)    
     &                 *recip_drf(k)
     &                 *diffkz(i,j,k+1,t)*recip_drc(k+1)
          enddo
         enddo
        enddo
C       main diagonal :
        do k=1,nr
         do j=1,ny
          do i=1,nx
            c5d(i,j,k) = 1.0000 - b5d(i,j,k) - d5d(i,j,k)
          enddo
         enddo
        enddo

C       get local copy of tracer for the tridiag solver
        do k=1,nr
         do j=1,ny
          do i=1,nx
           tmptracer(i,j,k)=tracer(i,j,k,t)
          enddo
         enddo
        enddo
      
C iterate a few times
      do iter=1,int(niters)
C       solve tri-diagonal system :
         call solve_matrix( nx, ny, nr, b5d, c5d, d5d, tmptracer)

        do k=1,nr
         do j=1,ny
          do i=1,nx
           itertracer(i,j,k)=tmptracer(i,j,k)
          enddo
         enddo
        enddo
C end of the iteration loop
       enddo
       
C      calculate diffusive flux
       do k=nr,1,-1
        if ( k.ge.2 ) then
         do j=1,ny
          do i=1,nx
           widiff(i,j,k,t) =
     &         -rac(i,j)*diffkz(i,j,k,t)*recip_drc(k)*rksign
     &         * (tmptracer(i,j,k) - tmptracer(i,j,k-1))
     &         * maskc(i,j,k)
     &         * maskc(i,j,k-1)
          enddo
         enddo
        else
         do j=1,ny
          do i=1,nx
           widiff(i,j,k,t) = 0.d0
          enddo
         enddo
        endif
       enddo

C     end of time loop
      enddo
      
      return
      end
      
      subroutine solve_matrix(nx, ny, nr, a3d, b3d, c3d, y3d)

C      used to solve implicitly vertical advection & diffusion
C     a3d :: matrix lower diagnonal
C     b3d :: matrix main  diagnonal
C     c3d :: matrix upper diagnonal
C     y3d :: input = y vector ; output = x = solution of a*x=y
      implicit none
C Local variable definitions
      integer nx, ny, nr, i,j,k
      real*8, intent(in), dimension(nx, ny, nr)      :: a3d
      real*8, intent(in), dimension(nx, ny, nr)      :: b3d
      real*8, intent(in), dimension(nx, ny, nr)      :: c3d
      real*8, intent(inout), dimension(nx, ny, nr)   :: y3d
      real*8 tmpvar, recvar 
      real*8 c3d_prime (nx,ny,nr)
      real*8 y3d_prime (nx,ny,nr)      
      real*8 y3d_update(nx,ny,nr)

      do k=1,nr
       do j=1,ny
        do i=1,nx
         c3d_prime(i,j,k) = 0.d0
         y3d_prime(i,j,k) = 0.d0
         y3d_update(i,j,k) = 0.d0
        enddo
       enddo
      enddo

C     forward sweep
      do k=1,nr
       do j=1,ny
        do i=1,nx
         if ( k .eq. 1 ) then
           if ( b3d(i,j,1) .ne. 0.d0 ) then
             recvar = 1.d0 / b3d(i,j,1)
             c3d_prime(i,j,1) = c3d(i,j,1)*recvar
             y3d_prime(i,j,1) = y3d(i,j,1)*recvar
           else
             c3d_prime(i,j,1) = 0.d0
             y3d_prime(i,j,1) = 0.d0
C             call mexErrMsgTxt(
C     &         'ERROR: forward sweep of implicit mixing matrix solver.')
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
             c3d_prime(i,j,k) = 0.d0
             y3d_prime(i,j,k) = 0.d0
C             call mexErrMsgTxt(
C     &         'ERROR: forward sweep of implicit mixing matrix solver.')
           endif
         endif
        enddo
       enddo
      enddo

C     backward sweep
      do k=nr,1,-1
       do j=1,ny
        do i=1,nx
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
       do j=1,ny
        do i=1,nx
         y3d(i,j,k) = y3d_update(i,j,k)
        enddo
       enddo
      enddo
      return
      end

      subroutine calc_vertical_int(nx,ny,nr,nt,klev,dz,
     &                       maskc,tracer,rfint)

      implicit none
      include "offline_consts.h"
C Local variable definitions
      integer i,j,k,t
      real*8, intent(in), dimension(nr)              :: dz    
      real*8, intent(in), dimension(nx, ny, nr)      :: maskc
      real*8, intent(in), dimension(nx, ny, nt)      :: klev
      real*8, intent(in), dimension(nx, ny, nr, nt)  :: tracer
      real*8, intent(out), dimension(nx, ny, nt)     :: rfint

C#ifdef DEBUG_MESSAGES 
C      integer, external :: mexPrintf
C      integer msg
C      character(len=250) :: line
C#endif

C     Initialize fluxes to zero
C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt('zero the outputs')    
C      write(line,*)'dimensions:\ni:',nx,'\nj:',ny,
C     &                        '\nk',nr,'\nt:',nt 
C      msg=mexPrintf(line//achar(13)) 
C#endif
      do t=1,nt
       do j=1,ny
        do i=1,nx  
         rfint(i,j,t)=0.D0
        enddo
       enddo  
      enddo

C     Integrate vertically
C#ifdef DEBUG_MESSAGES 
C      call mexWarnMsgTxt('integrate vertically')  
C#endif  
      do t=1,nt
       do j=1,ny
        do i=1,nx  
         do k=1,INT(klev(i,j,t))    
          rfint(i,j,t)=rfint(i,j,t)
     &      +(tracer(i,j,k,t)*dz(k)*maskc(i,j,k))
         enddo
        enddo
       enddo
      enddo
                  
      return
      end
