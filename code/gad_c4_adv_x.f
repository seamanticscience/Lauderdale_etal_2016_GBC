C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_c4_adv_x.F,v 1.7 2011/03/29 15:50:30 jmc Exp $

      SUBROUTINE GAD_C4_ADV_X(
     I           iMax,jMax,Nr,Nt,
     I           rA, maskC, maskU,
     I           uVel,tracer,
     O           uT )

C !DESCRIPTION:
C Calculates the area integrated zonal flux due to advection of a tracer
C using centered fourth-order interpolation:
C \begin{equation*}
C F^x_{adv} = U \overline{ \theta - \frac{1}{6} \delta_{ii} \theta }^i
C \end{equation*}
C Near boundaries, the scheme reduces to a second if the flow is away
C from the boundary and to third order if the flow is towards
C the boundary.

      IMPLICIT NONE

      integer*4 i,j,k,t,im1,im2,ip1
      integer*4 iMax,jMax,Nr,Nt

      real*8 maskC(iMax,jMax,Nr)
      real*8 maskU(iMax,jMax,Nr)
      real*8 rA(iMax,jMax,Nr)
      real*8 uVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 uT(iMax,jMax,Nr,Nt)

      real*8 Rjm,Rj,Rjp,Rjjm,Rjjp
      real*8 uTrans,half,oneSixth
      PARAMETER(half=1.D0/2.D0)
      PARAMETER(oneSixth=1.D0/6.D0)
      

      DO t=1,Nt    
       DO k=1,Nr
        DO j=1,jMax
         DO i=1,iMax           
          IF (i .EQ. 2) THEN
           im1=1
           im2=iMax
          ELSEIF (i .EQ. 1) THEN
           im1=iMax
           im2=iMax-1
          ELSE
           im1=i-1
           im2=i-2
          ENDIF

          IF (i .EQ. iMax) THEN
           ip1=1
          ELSE
           ip1=i+1
          ENDIF     
                    
          Rjp = (tracer(ip1,j,k,t)-tracer( i ,j,k,t))
C     &           *maskU(ip1,j,k)*maskC(ip1,j,k)
          Rj  = (tracer( i ,j,k,t)-tracer(im1,j,k,t))
C     &           *maskU( i ,j,k)*maskC( i ,j,k)
          Rjm = (tracer(im1,j,k,t)-tracer(im2,j,k,t))
C     &           *maskU(im1,j,k)*maskC(im1,j,k)
          Rjjp=(Rjp-Rj)
          Rjjm=(Rj-Rjm)
        
          uTrans=uVel(i,j,k,t)*rA(i,j,k)

C          uT(i,j,k,t) = maskC( i ,j,k)*maskC(im1,j,k)*
C     &                 *maskU( i ,j,k)*maskU(im1,j,k)*     
          uT(i,j,k,t) = 
     &     uTrans*(
     &       Tracer(i,j,k,t)+Tracer(im1,j,k,t)-oneSixth*( Rjjp+Rjjm )
     &             )*half
     &       +ABS( uTrans )*half*oneSixth*( Rjjp-Rjjm )
     &       *( 1.D0 - maskU(im1,j,k)*maskU(ip1,j,k) )
         ENDDO
        ENDDO
       ENDDO
      ENDDO  

      RETURN
      END
