C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_implicit_r.F,v 1.23 2015/06/03 13:39:22 rpa Exp $
      SUBROUTINE MIT_IMPLICIT_RMIX(
     I     iMax, jMax, Nr, Nt,
     I     deltaT, iterret, maskC, rA, dzf, dzc, 
     I     kappaZ, Tracer,
     O     df, 
C     O     iterTracer, 
     O     errCode )

C     Solve implicitly vertical advection and diffusion for one tracer.
      IMPLICIT NONE

C !INPUT/OUTPUT PARAMETERS:
C == Routine Arguments ==
C kappaZ   :: 3-D array for vertical diffusion coefficient
C maskC    :: inverse of cell open-depth factor
C Tracer  &
C tmpTracer:: tracer field at current time step
C newTracer:: future tracer field
C a5d      :: 2nd  lower diagonal of the pentadiagonal matrix
C b5d      :: 1rst lower diagonal of the pentadiagonal matrix
C c5d      :: main diagonal       of the pentadiagonal matrix
C d5d      :: 1rst upper diagonal of the pentadiagonal matrix
C e5d      :: 2nd  upper diagonal of the pentadiagonal matrix
C errCode  :: > 0 if singular matrix

      integer*4 iMax,jMax,Nr,Nt
      integer*4 i,j,k,t,iter
      integer*4 errCode 
      real*8 iterret
      real*8 deltaT
      real*8 recip_drF(Nr)
      real*8 recip_drC(Nr)
      real*8 dzf      (Nr)
      real*8 dzc      (Nr)
      real*8 rA       (iMax,jMax)
      real*8 maskC    (iMax,jMax,Nr)
      real*8 tmpTracer(iMax,jMax,Nr)
      real*8 iterTracer(iMax,jMax,Nr)
      real*8 Tracer   (iMax,jMax,Nr,Nt)
      real*8 kappaZ   (iMax,jMax,Nr,Nt)
C      real*8 iterTracer(iMax,jMax,Nr,Nt,15)
      real*8 df       (iMax,jMax,Nr,Nt)
      real*8 b5d      (iMax,jMax,Nr)
      real*8 c5d      (iMax,jMax,Nr)
      real*8 d5d      (iMax,jMax,Nr)    

C From config_summary
C rkSign =   /* index orientation relative to vertical coordinate */
C -1.000000000000000E+00
      real*8 rkSign
      PARAMETER(rkSign=-1.D0)

      DO t=1,Nt

C--   Initialise
       DO k=1,Nr        
        recip_drF(k)=1.D0/dzf(k)
        recip_drC(k)=1.D0/dzc(k)
        DO j=1,jMax
         DO i=1,iMax
          b5d(i,j,k) = 0.D0
          c5d(i,j,k) = 1.D0
          d5d(i,j,k) = 0.D0
         ENDDO
        ENDDO
       ENDDO

C--   set the tri-diagonal matrix to solve the implicit diffusion problem
C-     1rst lower diagonal :
        DO k=2,Nr
         DO j=1,jMax
          DO i=1,iMax
            b5d(i,j,k) =-deltaT
     &                  *maskC(i,j,k-1)*maskC (i,j,k)  
     &                  *recip_drF(k)
     &                  *kappaZ(i,j,k,t)*recip_drC(k)
          ENDDO
         ENDDO
        ENDDO
C-     1rst upper diagonal :
        DO k=1,Nr-1
         DO j=1,jMax
          DO i=1,iMax
            d5d(i,j,k) = -deltaT
     &                 *maskC(i,j,k+1)*maskC (i,j,k)    
     &                 *recip_drF(k)
     &                 *kappaZ(i,j,k+1,t)*recip_drC(k+1)
          ENDDO
         ENDDO
        ENDDO
C-     Main diagonal :
        DO k=1,Nr
         DO j=1,jMax
          DO i=1,iMax
            c5d(i,j,k) = 1.0000 - b5d(i,j,k) - d5d(i,j,k)
          ENDDO
         ENDDO
        ENDDO

C       Get local copy of Tracer for the tridiag solver
        DO k=1,Nr
         DO j=1,jMax
          DO i=1,iMax
           tmpTracer(i,j,k)=Tracer(i,j,k,t)
          ENDDO
         ENDDO
        ENDDO
      
C Iterate a few times
      DO iter=1,int(iterret)
C       Solve tri-diagonal system :
         CALL SOLVE_TRIDIAGONAL( iMax, jMax, Nr,
     I                           b5d, c5d, d5d,
     U                           tmpTracer,
     O                           errCode)

        DO k=1,Nr
         DO j=1,jMax
          DO i=1,iMax
           iterTracer(i,j,k)=tmpTracer(i,j,k)
          ENDDO
         ENDDO
        ENDDO
C end of the iteration loop
       ENDDO
       
C      Calculate diffusive flux
       DO k=Nr,1,-1
        IF ( k.GE.2 ) THEN
         DO j=1,jMax
          DO i=1,iMax
           df(i,j,k,t) =
     &         -rA(i,j)*KappaZ(i,j,k,t)*recip_drC(k)*rkSign
     &         * (tmpTracer(i,j,k) - tmpTracer(i,j,k-1))
     &         * maskC(i,j,k)
     &         * maskC(i,j,k-1)
          ENDDO
         ENDDO
        ELSE
         DO j=1,jMax
          DO i=1,iMax
           df(i,j,k,t) = 0.D0
          ENDDO
         ENDDO
        ENDIF
       ENDDO

C     End of time loop
      ENDDO
      
      RETURN
      END