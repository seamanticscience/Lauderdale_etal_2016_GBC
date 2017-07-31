C $Header: /u/gcmpack/MITgcm/pkg/gmredi/gmredi_rtransport.F,v 1.19 2014/09/09 22:34:06 jmc Exp $
      SUBROUTINE MIT_DIFF_Z(
     I     iMax, jMax, Nr, Nt,
     I     maskU, maskV, maskC, 
     I     dx, dy, rA, Kwx, Kwy,
     I     Tracer,
     U     df )

C     *==========================================================*
C     | o SUBROUTINE GMREDI_RTRANSPORT                           |
C     |   Add vertical transport terms from GM/Redi              |
C     |   parameterization.                                      |
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     !INPUT/OUTPUT PARAMETERS:
C     1,iMax    :: Range of 1rst index where results will be set
C     1,jMax    :: Range of 2nd  index where results will be set
C     Tracer       :: 3D Tracer field
C     df           :: Diffusive flux component work array.

      INTEGER iMax, jMax, Nr, Nt
      INTEGER i, j, k, t
      INTEGER kp1, km1, jm1, im1, ip1, jp1
      real*8 dx      (iMax,jMax)
      real*8 dy      (iMax,jMax)
      real*8 recip_dy(iMax,jMax)
      real*8 recip_dx(iMax,jMax)
      real*8 rA      (iMax,jMax)
      real*8 maskU   (iMax,jMax,Nr)
      real*8 maskV   (iMax,jMax,Nr)
      real*8 maskC   (iMax,jMax,Nr)
      real*8 Kwx     (iMax,jMax,Nr,Nt)
      real*8 Kwy     (iMax,jMax,Nr,Nt)
      real*8 Tracer  (iMax,jMax,Nr,Nt)
      real*8 df      (iMax,jMax,Nr,Nt)

      real*8 dTdx, dTdy
      real*8 half
      PARAMETER(half=1.D0/2.D0)

C     Initialise flux to zero
      DO t=1,Nt
       DO k=1,Nr
        DO j=1,jMax
         DO i=1,iMax
          df(i,j,k,t) = 0.D0
          
          IF (maskC(i,j,1) .NE. 0.D0) THEN
           recip_dx(i,j)=(1.D0/dx(i,j))
           recip_dy(i,j)=(1.D0/dy(i,j))
          ELSE
           recip_dx(i,j)=0.D0
           recip_dy(i,j)=0.D0
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C-    Horizontal gradients interpolated to W points
      DO t=1,Nt
       DO k=2,Nr
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

         km1 = MAX(k-1,1)
         kp1 = MIN(k+1,Nr)
                
         dTdx = half*(
     &    half*(maskU(ip1,j,k)
     &         *recip_dx(ip1,j)*(Tracer(ip1,j,k,t)-Tracer(i,j,k,t))
     &        +maskU(i,j,k)
     &         *recip_dx(i,j)*(Tracer(i,j,k,t)-Tracer(im1,j,k,t))
     &         )
     &   +half*(maskU(ip1,j,km1)
     &         *recip_dx(ip1,j)*(Tracer(ip1,j,km1,t)-Tracer(i,j,km1,t))
     &        +maskU(i,j,km1)
     &         *recip_dx(i,j)*(Tracer(i,j,km1,t)-Tracer(im1,j,km1,t)) 
     &         )
     &         )

         jm1=MAX(j-1, 1 )
         jp1=MIN(j+1,jMax)

         dTdy = half*(
     &    half*(maskV(i,j,k)
     &         *recip_dy(i,j)*(Tracer(i,j,k,t)-Tracer(i,jm1,k,t))
     &        +maskV(i,j+1,k)
     &         *recip_dy(i,jp1)*(Tracer(i,jp1,k,t)-Tracer(i,j,k,t))
     &         )
     &   +half*(maskV(i,j,km1)
     &         *recip_dy(i,j)*(Tracer(i,j,km1,t)-Tracer(i,jm1,km1,t))
     &        +maskV(i,jp1,km1)
     &         *recip_dy(i,jp1)*(Tracer(i,jp1,km1,t)-Tracer(i,j,km1,t))
     &         )
     &         )

C-    Off-diagonal components of vertical flux
         df(i,j,k,t) = df(i,j,k,t)
     &       -rA(i,j)*maskC(i,j,k)
     &       *( Kwx(i,j,k,t)*dTdx
     &        + Kwy(i,j,k,t)*dTdy )
         ENDDO
        ENDDO
       ENDDO
      ENDDO 

c     IF (.NOT.implicitDiffusion) THEN
C This vertical diffusion term is currently implemented
C by adding the VisbeckK*Kwz diffusivity to KappaRT/S
C See calc_diffusivity.F and calc_gt.F (calc_gs.F)
c      DO j=1,jMax
c       DO i=1,iMax
c        df(i,j) = df(i,j) - _rA(i,j,bi,bj)
c    &    *maskUp(i,j)*VisbeckK(i,j,bi,bj)*Kwz(i,j,k,bi,bj)
c    &    *recip_drC(k)*rkfac
c    &    *(Tracer(i,j,k-1)-Tracer(i,j,k))
c    &    *maskInC(i,j,bi,bj)
c       ENDDO
c      ENDDO
c     ENDIF

      RETURN
      END
