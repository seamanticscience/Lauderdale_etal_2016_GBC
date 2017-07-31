C $Header: /u/gcmpack/MITgcm/pkg/gmredi/gmredi_xtransport.F,v 1.21 2014/09/09 22:34:06 jmc Exp $
      SUBROUTINE MIT_DIFF_X(
     I     iMax, jMax, Nr, Nt,
     I     maskC, dz, dx, xA, 
     I     Kux, Kuz, Tracer,
     U     df )

C     *==========================================================*
C     | o SUBROUTINE GMREDI_XTRANSPORT
C     |   Add horizontal x transport terms from GM/Redi
C     |   parameterization.
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     !INPUT/OUTPUT PARAMETERS:
C     1,iMax    :: Range of 1rst index where results will be set
C     1,jMax    :: Range of 2nd  index where results will be set
C     xA           :: Area of X face
C     Tracer       :: 3D Tracer field
C     df           :: Diffusive flux component work array.
      INTEGER iMax, jMax, Nr, Nt
      INTEGER i, j, k, t, kp1, km1, im1

      real*8 dz      (Nr)
      real*8 recip_dz(Nr)
      real*8 dx      (iMax,jMax)
      real*8 recip_dx(iMax,jMax)
      real*8 maskC   (iMax,jMax,Nr)
      real*8 xA      (iMax,jMax,Nr)
      real*8 Kux     (iMax,jMax,Nr,Nt)
      real*8 Kuz     (iMax,jMax,Nr,Nt)
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
           recip_dx(i,j)=1.D0/dx(i,j)
          ELSE
           recip_dx(i,j)=0.D0
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
          IF (i .EQ. 1) THEN
           im1=iMax
          ELSE
           im1=i-1
          ENDIF
          df(i,j,k,t) = df(i,j,k,t)
     &      -xA(i,j,k)
     &      *Kux(i,j,k,t)
     &      *recip_dx(i,j)
     &      *(Tracer(i,j,k,t)-Tracer(im1,j,k,t))
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      
C-    Off-diagonal components of horizontal flux
C     If not GM_EXTRA_DIAGONAL (GM_ExtraDiag) then will be zero
      DO t=1,Nt
       DO k=1,Nr
       km1 = MAX(k-1,1)
       kp1 = MIN(k+1,Nr)
C-    Vertical gradients interpolated to U points
        DO j=1,jMax
         DO i=1,iMax
           IF (i .EQ. 1) THEN
            im1=iMax
           ELSE
            im1=i-1
           ENDIF
             
           dTdz =  half*(
     &      half*recip_dz(k)*
     &        (maskC(im1,j,k)*(Tracer(im1,j,km1,t)-Tracer(im1,j,k,t))
     &        +maskC( i ,j,k)*(Tracer( i ,j,km1,t)-Tracer( i ,j,k,t))
     &        )
     &      +half*recip_dz(kp1)*
     &        (maskC(im1,j,kp1)*(Tracer(im1,j,k,t)-Tracer(im1,j,kp1,t))
     &        +maskC( i ,j,kp1)*(Tracer( i ,j,k,t)-Tracer( i ,j,kp1,t))
     &        )
     &      )
  
           df(i,j,k,t) = df(i,j,k,t) - xA(i,j,k)*Kuz(i,j,k,t)*dTdz
         ENDDO
        ENDDO
       ENDDO
      ENDDO
       
      RETURN
      END
