C $Header: /u/gcmpack/MITgcm/pkg/gmredi/gmredi_ytransport.F,v 1.21 2014/09/09 22:34:06 jmc Exp $
      SUBROUTINE MIT_DIFF_Y(
     I     iMax, jMax, Nr, Nt,
     I     maskC, dz, dy, yA, 
     I     Kvy, Kvz, Tracer,
     U     df )

C     *==========================================================*
C     | o SUBROUTINE GMREDI_YTRANSPORT
C     |   Add horizontal y transport terms from GM/Redi
C     |   parameterization.
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     !INPUT/OUTPUT PARAMETERS:
C     1,iMax    :: Range of 1rst index where results will be set
C     1,jMax    :: Range of 2nd  index where results will be set
C     yA           :: Area of Y face
C     Tracer       :: 3D Tracer field
C     df           :: Diffusive flux component work array.
      INTEGER iMax, jMax, Nr, Nt
      INTEGER i, j, k, t, kp1, km1, jm1

      real*8 dz      (Nr)
      real*8 recip_dz(Nr)
      real*8 dy      (iMax,jMax)
      real*8 recip_dy(iMax,jMax)
      real*8 yA      (iMax,jMax,Nr)
      real*8 maskC   (iMax,jMax,Nr)
      real*8 Kvy     (iMax,jMax,Nr,Nt)
      real*8 Kvz     (iMax,jMax,Nr,Nt)
      real*8 Tracer  (iMax,jMax,Nr,Nt)
      real*8 df      (iMax,jMax,Nr,Nt)

      real*8 dTdz
      real*8 half
      PARAMETER(half=1.D0/2.D0)

C     Initialise flux to zero
      DO t=1,Nt
       DO k=1,Nr
         recip_dz(k)=1.D0/dz(k)  
        DO j=1,jMax
         DO i=1,iMax
          df(i,j,k,t) = 0.D0
       
          IF (maskC(i,j,1) .NE. 0.D0) THEN
           recip_dy(i,j)=1.D0/dy(i,j)
          ELSE
           recip_dy(i,j)=0.D0
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      
C--   Area integrated zonal flux
      DO t=1,Nt
       DO k=1,Nr
        DO j=1,jMax
         DO i=1,iMax
          jm1=MAX(j-1, 1 )
          df(i,j,k,t) = df(i,j,k,t)
     &      -yA(i,j,k)
     &      *Kvy(i,j,k,t)
     &      *recip_dy(i,j)
     &      *(Tracer(i,j,k,t)-Tracer(i,jm1,k,t))
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      
      DO t=1,Nt
       DO k=1,Nr
       km1 = MAX(k-1,1)
       kp1 = MIN(k+1,Nr)
C-    Vertical gradients interpolated to U points
        DO j=1,jMax
         DO i=1,iMax

          jm1=MAX(j-1, 1 )

C-    Vertical gradients interpolated to V points
           dTdz =  half*(
     &     half*recip_dz(k)*
     &        (maskC(i,jm1,k)*(Tracer(i,jm1,km1,t)-Tracer(i,jm1,k,t))
     &        +maskC(i, j ,k)*(Tracer(i, j ,km1,t)-Tracer(i, j ,k,t))
     &        )
     &     +half*recip_dz(kp1)*
     &        (maskC(i,jm1,kp1)*(Tracer(i,jm1,k,t)-Tracer(i,jm1,kp1,t))
     &        +maskC(i, j ,kp1)*(Tracer(i, j ,k,t)-Tracer(i, j ,kp1,t))
     &        )
     &     )

           df(i,j,k,t) = df(i,j,k,t) - yA(i,j,k)*Kvz(i,j,k,t)*dTdz
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      
      RETURN
      END
