C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_dst3_adv_x.F,v 1.15 2012/07/04 20:22:26 jmc Exp $

      SUBROUTINE GAD_DST3_ADV_X( iMax,jMax,Nr,Nt,rA, 
     &              maskC,maskU,deltaT,deltaX,uVel,tracer,uT )
   
C !DESCRIPTION:
C  Calculates the area integrated zonal flux due to advection of a
C  tracer using 3rd-order Direct Space and Time (DST-3) Advection Scheme

      IMPLICIT NONE

      integer*4 i,j,k,t,im1,im2
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

      real*8 Rjm,Rj,Rjp,d0,d1,oneSixth
      real*8 uTrans,uCFL,half
      PARAMETER(half=1.D0/2.D0)
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
     &     *maskU( i ,j,k)*maskC( i ,j,k)
          Rjm=(tracer(im1,j,k,t)-tracer(im2,j,k,t))
     &     *maskC(im1,j,k)*maskC(im2,j,k)

          uCFL = ABS( uVel(i,j,k,t)*deltaT*recipX(i,j) )
          uTrans=uVel(i,j,k,t)*rA(i,j,k)
          
          d0=(2.D0-uCFL)*(1.D0-uCFL)*oneSixth
          d1=(1.D0-uCFL*uCFL)*oneSixth

          uT(i,j,k,t)=maskC( i ,j,k)*maskC(im1,j,k)
C     &              *maskU( i ,j,k)*maskU(im1,j,k)
C          uT(i,j,k,t)= (
     &   *(half*(uTrans+ABS(uTrans))
     &      *( Tracer(im1,j,k,t) + (d0*Rj+d1*Rjm) )
     &    +half*(uTrans-ABS(uTrans))
     &      *( Tracer( i ,j,k,t) - (d0*Rj+d1*Rjp) ))

         ENDDO
        ENDDO
       ENDDO
      ENDDO  

      RETURN
      END
