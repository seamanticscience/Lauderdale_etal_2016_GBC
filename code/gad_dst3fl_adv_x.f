C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_dst3fl_adv_x.F,v 1.16 2014/04/04 20:29:08 jmc Exp $
      SUBROUTINE GAD_DST3FL_ADV_X( iMax,jMax,Nr,Nt,rA, 
     &              maskC,maskU,deltaT,deltaX,uVel,tracer,uT )
     
C     /==========================================================\
C     | SUBROUTINE GAD_DST3FL_ADV_X                              |
C     | o Compute Zonal advective Flux of Tracer using           |
C     |   3rd Order DST Sceheme with flux limiting               |
C     |==========================================================|

      IMPLICIT NONE

      integer*4 i,j,k,t,im1,im2,ip1
      integer*4 iMax,jMax,Nr,Nt
      
      real*8 deltaT
      real*8 deltaX(iMax,jMax)
      real*8 recipX(iMax,jMax)
      real*8 rA(iMax,jMax,Nr)
      real*8 maskC(iMax,jMax,Nr)
      real*8 maskU(iMax,jMax,Nr)
      real*8 uVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 uT(iMax,jMax,Nr,Nt)
C      real*8 af(iMax,jMax,Nr,Nt)

      real*8 Rjm,Rj,Rjp,d0,d1,psiP,psiM
      real*8 thetaMin,thetaMax,thetaP,thetaM
      real*8 uTrans,uCFL,oneSixth,half
      PARAMETER(half=1.D0/2.D0)
      PARAMETER(thetaMax = 1.D+20 )
      PARAMETER(thetaMin = 1.D-20 )
      PARAMETER(oneSixth=1.D0/6.D0)

      DO t=1,Nt    
       DO k=1,Nr
        DO j=1,jMax
         DO i=1,iMax
          uT(i,j,k,t) = 0.d0
         ENDDO
        ENDDO
       ENDDO
      ENDDO        

      DO j=1,jMax
       DO i=1,iMax
          recipX(i,j)=1/deltaX(i,j)
       ENDDO
      ENDDO
      
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
          
          Rjp=(tracer(im1,j,k,t)-tracer( i ,j,k,t))
     &     *maskC( i ,j,k)*maskC(im1,j,k)
          Rj =(tracer( i ,j,k,t)-tracer(im1,j,k,t))
     &     *maskC( i ,j,k)*maskC(im1,j,k)
          Rjm=(tracer(im1,j,k,t)-tracer(im2,j,k,t))
     &     *maskC(im1,j,k)*maskC(im2,j,k)

          uCFL = ABS( uVel(i,j,k,t)*deltaT*recipX(i,j) )
          uTrans=uVel(i,j,k,t)*rA(i,j,k)
          
          d0=(2.D0-uCFL)*(1.D0-uCFL)*oneSixth
          d1=(1.D0-uCFL*uCFL)*oneSixth
          
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
     &          thetaP*(1.D0 -uCFL)/(uCFL+thetaMin) ))
           psiM=d0+d1*thetaM
           psiM=MAX(0.D0,MIN(MIN(1.D0,psiM),
     &          thetaM*(1.D0 -uCFL)/(uCFL+thetaMin) ))

           uT(i,j,k,t)=maskC( i ,j,k)*maskC(im1,j,k)
C     &                *maskU( i ,j,k)*maskU(im1,j,k)
C           uT(i,j,k,t)=(
     &   *(half*(uTrans+ABS(uTrans))
     &      *( tracer(i,j, k ,t) + psiM*Rj )
     &  +half*(uTrans-ABS(uTrans))
     &      *( tracer(im1,j,k,t) - psiP*Rj ))
         ENDDO
        ENDDO
       ENDDO
      ENDDO  

      RETURN
      END
