C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_dst3fl_adv_y.F,v 1.15 2014/04/04 20:29:08 jmc Exp $

      SUBROUTINE GAD_DST3FL_ADV_Y( iMax,jMax,Nr,Nt,rA, 
     &              maskC,maskV,deltaT,deltaY,vVel,tracer,vT )

C     /==========================================================\
C     | SUBROUTINE GAD_DST3FL_ADV_Y                              |
C     | o Compute Meridional advective Flux of Tracer using      |
C     |   3rd Order DST Sceheme with flux limiting               |
C     |==========================================================|
     
      IMPLICIT NONE

      integer*4 i,j,k,t,jm1,jm2,jp1
      integer*4 iMax,jMax,Nr,Nt
      
      real*8 deltaT
      real*8 deltaY(iMax,jMax)
      real*8 recipY(iMax,jMax)
      real*8 rA(iMax,jMax,Nr)
      real*8 maskC(iMax,jMax,Nr)
      real*8 maskV(iMax,jMax,Nr)
      real*8 vVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 vT(iMax,jMax,Nr,Nt)
C      real*8 af(iMax,jMax,Nr,Nt)

      real*8 Rjm,Rj,Rjp,d0,d1,psiP,psiM
      real*8 thetaMin,thetaMax,thetaP,thetaM
      real*8 vTrans,vCFL,oneSixth,half
      PARAMETER(half=1.D0/2.D0)
      PARAMETER(thetaMax = 1.D+20 )
      PARAMETER(thetaMin = 1.D-20 )
      PARAMETER(oneSixth=1.D0/6.D0)

      DO t=1,Nt    
       DO k=1,Nr
        DO j=1,jMax
         DO i=1,iMax
          vT(i,j,k,t) = 0.d0
         ENDDO
        ENDDO
       ENDDO
      ENDDO        

      DO j=1,jMax
       DO i=1,iMax
          recipY(i,j)=1/deltaY(i,j)
       ENDDO
      ENDDO
      
      DO t=1,Nt    
       DO k=1,Nr
        DO j=1,jMax
         DO i=1,iMax           
          jp1=MIN(j+1,jMax)
          jm1=MAX(j-1, 1  )
          jm2=MAX(j-2, 1  )
          
          Rjp=(tracer(i,jp1,k,t)-tracer(i, j ,k,t))
     &     *maskC(i, j ,k)*maskC(i,jp1,k)
          Rj =(tracer(i, j ,k,t)-tracer(i,jm1,k,t))
     &     *maskC(i, j ,k)*maskC(i,jm1,k)
          Rjm=(tracer(i,jm1,k,t)-tracer(i,jm2,k,t))
     &     *maskC(i,jm1,k)*maskC(i,jm2,k)

          vCFL = ABS( vVel(i,j,k,t)*deltaT*recipY(i,j) )
          vTrans=vVel(i,j,k,t)*rA(i,j,k)
          
          d0=(2.D0-vCFL)*(1.D0-vCFL)*oneSixth
          d1=(1.D0-vCFL*vCFL)*oneSixth

C-      the old version: can produce overflow, division by zero,
c       and is wrong for tracer with low concentration:
c       thetaP=Rjm/(1.D-20+Rj)
c       thetaM=Rjp/(1.D-20+Rj)
C-      the right expression, but not bounded:
c       thetaP=0.D0
c       thetaM=0.D0
c       IF (Rj.NE.0.D0) thetaP=Rjm/Rj
c       IF (Rj.NE.0.D0) thetaM=Rjp/Rj
C-      prevent |thetaP,M| to reach too big value:
           IF ( ABS(Rj)*thetaMax .LE. ABS(Rjm) ) THEN
             thetaP=SIGN(thetaMax,Rjm*Rj)
           ELSE
             thetaP=Rjm/Rj
           ENDIF
           IF ( ABS(Rj)*thetaMax .LE. ABS(Rjp) ) THEN
             thetaM=SIGN(thetaMax,Rjp*Rj)
           ELSE
             thetaM=Rjp/Rj
           ENDIF

           psiP=d0+d1*thetaP
           psiP=MAX(0.D0,MIN(MIN(1.D0,psiP),
     &          thetaP*(1.D0 -vCFL)/(vCFL+thetaMin) ))
           psiM=d0+d1*thetaM
           psiM=MAX(0.D0,MIN(MIN(1.D0,psiM),
     &          thetaM*(1.D0 -vCFL)/(vCFL+thetaMin) ))

           vT(i,j,k,t)=maskC(i, j ,k)*maskC(i,jm1,k)
C     &                *maskV(i, j ,k)*maskV(i,jm1,k)
C           vT(i,j,k,t)=(      
     &   *( half*(vTrans+ABS(vTrans))
     &      *( tracer(i, j ,k,t) + psiM*Rj )
     &  +half*(vTrans-ABS(vTrans))
     &      *( tracer(i,jm1,k,t) - psiP*Rj ))
         ENDDO
        ENDDO
       ENDDO
      ENDDO  
      
      RETURN
      END
