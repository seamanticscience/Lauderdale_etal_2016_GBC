C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_c2_adv_x.F,v 1.4 2008/02/29 030:59 mlosch Exp $

      SUBROUTINE GAD_UPWIND_X( iMax,jMax,Nr,Nt,rA,maskC,
     &           maskU,uVel,tracer,uT )

C !DESCRIPTION:
C Calculates the area integrated zonal flux due to advection of a tracer using
C centered second-order interpolation:
C \begin{equation*}
C F^x_{adv} = U \overline{\theta}^i
C \end{equation*}

      IMPLICIT NONE

      integer*4 i,j,k,t,im1
      integer*4 iMax,jMax,Nr,Nt
      real*8 maskC(iMax,jMax,Nr)
      real*8 maskU(iMax,jMax,Nr)
      real*8 rA(iMax,jMax,Nr)
      real*8 uVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 uT(iMax,jMax,Nr,Nt)
C      real*8 af(iMax,jMax,Nr,Nt)
      
         DO t=1,Nt
          DO k=1,Nr
           DO j=1,jMax
            DO i=1,iMax
             IF (i .EQ. 1) THEN
                im1=iMax
             ELSE
                im1=i-1
             ENDIF
             
             IF (uVel(i,j,k,t) .GE. 0.D0) THEN
C                uT(i,j,k,t) = maskU(i,j,k)*maskU(im1,j,k )*
C     &                        maskC(i,j,k)*maskC(im1,j,k )* 
                uT(i,j,k,t) = 
     &                        uVel(i,j,k,t)*tracer(im1,j,k,t)*rA(i,j,k)
             ELSE
C                uT(i,j,k,t) = maskU(i,j,k )*maskC(i,j,k)*
                uT(i,j,k,t) = 
     &                        uVel(i,j,k,t)*tracer( i ,j,k,t)*rA(i,j,k)
             ENDIF

            ENDDO
           ENDDO
          ENDDO
         ENDDO

      RETURN
      END
