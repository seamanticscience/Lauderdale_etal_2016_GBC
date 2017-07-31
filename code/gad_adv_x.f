C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_c2_adv_x.F,v 1.4 2008/02/29 030:59 mlosch Exp $

      SUBROUTINE GAD_ADV_X( iMax,jMax,Nr,Nt,rA,maskC,maskU,
     &           uVel,tracer,uT )

C !DESCRIPTION:
C Calculates the area integrated zonal flux due to advection of a tracer using
C centered second-order interpolation:
C \begin{equation*}
C F^x_{adv} = U \overline{\theta}^i
C \end{equation*}

      IMPLICIT NONE

      integer*4 i,j,k,t,im1,ip1
      integer*4 iMax,jMax,Nr,Nt
      real*8 maskC(iMax,jMax,Nr)
      real*8 maskU(iMax,jMax,Nr)
      real*8 rA(iMax,jMax,Nr)
      real*8 uVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 uT(iMax,jMax,Nr,Nt)
C      real*8 af(iMax,jMax,Nr,Nt)
      real*8 half
      PARAMETER(half=1.D0/2.D0)
      
         DO t=1,Nt
          DO k=1,Nr
           DO j=1,jMax
            DO i=1,iMax
             IF (i .EQ. 1) THEN
                im1=iMax
             ELSE
                im1=i-1
             ENDIF
             
             IF (i .EQ. iMax) THEN
                ip1=1
             ELSE
                ip1=i+1
             ENDIF
             
             uT(i,j,k,t) = uVel(i,j,k,t)*rA(i,j,k)
     &        *(Tracer(i,j,k,t)+Tracer(im1,j,k,t))*half             

C             uT(i,j,k,t) = maskC(i,j,k)*maskC(im1,j,k )*
C     &                     maskU(i,j,k)*maskU(im1,j,k )*         
C     &                     uVel(i,j,k,t)*rA(i,j,k)
c     &        *(Tracer(i,j,k,t)+Tracer(im1,j,k,t))*half
            ENDDO
           ENDDO
          ENDDO
         ENDDO

      RETURN
      END
