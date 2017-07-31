C----------------------------------------------------------------------------------------------
C $Header: /u/gcmpack/MITgcm/model/src/solve_tridiagonal.F,v 1.11 2014/08/14 16:49:19 jmc Exp $

      SUBROUTINE SOLVE_TRIDIAGONAL(
     I                     iMax,jMax,Nr,
     I                     a3d, b3d, c3d,
     U                     y3d, 
     O                     errCode )
C     *==========================================================*
C     | S/R SOLVE_TRIDIAGONAL
C     | o Solve a tri-diagonal system A*X=Y (dimension Nr)
C     *==========================================================*
C     | o Used to solve implicitly vertical advection & diffusion
C     *==========================================================*
      IMPLICIT NONE

C     !INPUT/OUTPUT PARAMETERS:
C     a3d :: matrix lower diagnonal
C     b3d :: matrix main  diagnonal
C     c3d :: matrix upper diagnonal
C     y3d :: Input = Y vector ; Output = X = solution of A*X=Y
C     errCode :: > 0 if singular matrix
      integer*4 iMax,jMax,Nr,i,j,k
      real*8 a3d(iMax,jMax,Nr)
      real*8 b3d(iMax,jMax,Nr)
      real*8 c3d(iMax,jMax,Nr)
      real*8 y3d(iMax,jMax,Nr)
      real*8 tmpVar, recVar 
      integer*4 errCode
C# ifdef SOLVE_DIAGONAL_KINNER
C      real*8 c3d_prime (Nr)
C      real*8 y3d_prime (Nr)
C      real*8 y3d_update(Nr)
C# else
      real*8 c3d_prime (iMax,jMax,Nr)
      real*8 y3d_prime (iMax,jMax,Nr)      
      real*8 y3d_update(iMax,jMax,Nr)
C# endif
      
      errCode = 0

C#ifdef SOLVE_DIAGONAL_KINNER
C
CC--   Main loop
C      DO j=1,jMax
C       DO i=1,iMax
C
C        DO k=1,Nr
C         c3d_prime(k) = 0.D0
C         y3d_prime(k) = 0.D0
C         y3d_update(k) = 0.D0
C        ENDDO
C
CC--   Forward sweep
C        DO k=1,Nr
C         IF ( k.EQ.1 ) THEN
C           IF ( b3d(i,j,1).NE.0.D0 ) THEN
C             c3d_prime(1) = c3d(i,j,1) / b3d(i,j,1)
C             y3d_prime(1) = y3d(i,j,1) / b3d(i,j,1)
C           ELSE
C             c3d_prime(1) = 0.D0
C             y3d_prime(1) = 0.D0
C             errCode = 1
C           ENDIF
C         ELSE
C           tmpVar = b3d(i,j,k) - a3d(i,j,k)*c3d_prime(k-1)
C           IF ( tmpVar .NE. 0.D0 ) THEN
C             recVar = 1.D0 / tmpVar
C             c3d_prime(k) = c3d(i,j,k)*recVar
C             y3d_prime(k) = (y3d(i,j,k) - y3d_prime(k-1)*a3d(i,j,k))
C     &                      *recVar
C           ELSE
C             c3d_prime(k) = 0.D0
C             y3d_prime(k) = 0.D0
C             errCode = 1
C           ENDIF
C         ENDIF
C        ENDDO
C
CC--   Backward sweep
C        DO k=Nr,1,-1
C         IF ( k.EQ.Nr ) THEN
C          y3d_update(k)=y3d_prime(k)
C         ELSE
C          y3d_update(k)=y3d_prime(k)-c3d_prime(k)*y3d_update(k+1)
C         ENDIF
C        ENDDO
C
CC--   Update array
C        DO k=1,Nr
C         y3d(i,j,k) = y3d_update(k)
C        ENDDO
C
C       ENDDO
C      ENDDO
C
C#else  /* ndef SOLVE_DIAGONAL_KINNER */

C--   Init.
      DO k=1,Nr
       DO j=1,jMax
        DO i=1,iMax
         c3d_prime(i,j,k) = 0.D0
         y3d_prime(i,j,k) = 0.D0
         y3d_update(i,j,k) = 0.D0
        ENDDO
       ENDDO
      ENDDO

C--   Forward sweep
      DO k=1,Nr
       DO j=1,jMax
        DO i=1,iMax
         IF ( k .EQ. 1 ) THEN
           IF ( b3d(i,j,1) .NE. 0.D0 ) THEN
             recVar = 1.D0 / b3d(i,j,1)
             c3d_prime(i,j,1) = c3d(i,j,1)*recVar
             y3d_prime(i,j,1) = y3d(i,j,1)*recVar
           ELSE
             c3d_prime(i,j,1) = 0.D0
             y3d_prime(i,j,1) = 0.D0
             errCode = 1
           ENDIF
         ELSE
           tmpVar = b3d(i,j,k) - a3d(i,j,k)*c3d_prime(i,j,k-1)
           IF ( tmpVar .NE. 0.D0 ) THEN
             recVar = 1.D0 / tmpVar
             c3d_prime(i,j,k) = c3d(i,j,k)*recVar
             y3d_prime(i,j,k) = ( y3d(i,j,k)
     &                          - a3d(i,j,k)*y3d_prime(i,j,k-1)
     &                          )*recVar
           ELSE
             c3d_prime(i,j,k) = 0.D0
             y3d_prime(i,j,k) = 0.D0
             errCode = 2
           ENDIF
         ENDIF
        ENDDO
       ENDDO
      ENDDO

C--   Backward sweep
      DO k=Nr,1,-1
       DO j=1,jMax
        DO i=1,iMax
         IF ( k .EQ. Nr ) THEN
          y3d_update(i,j,k) = y3d_prime(i,j,k)
         ELSE
          y3d_update(i,j,k) = y3d_prime(i,j,k)
     &                      - c3d_prime(i,j,k)*y3d_update(i,j,k+1)
         ENDIF
        ENDDO
       ENDDO
      ENDDO

CC--   Update array
      DO k=1,Nr
       DO j=1,jMax
        DO i=1,iMax
         y3d(i,j,k) = y3d_update(i,j,k)
        ENDDO
       ENDDO
      ENDDO
C#endif /* SOLVE_DIAGONAL_KINNER */
      RETURN
      END
