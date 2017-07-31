C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_dst3_adv_y.F,v 1.14 2012/07/04 20:22:26 jmc Exp $

      SUBROUTINE GAD_DST3_ADV_Y( iMax,jMax,Nr,Nt,rA, 
     &              maskC,maskV,deltaT,deltaY,vVel,tracer,vT )

C !DESCRIPTION:
C  Calculates the area integrated Meridional flux due to advection of a
C  tracer using 3rd-order Direct Space and Time (DST-3) Advection Scheme

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

      real*8 Rjm,Rj,Rjp,d0,d1,oneSixth
      real*8 vTrans,vCFL,half
      PARAMETER(half=1.D0/2.D0)
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

          vT(i,j,k,t)=maskC(i, j ,k)*maskC(i,jm1,k)
C     &               *maskV(i, j ,k)*maskV(i,jm1,k)   
C          vT(i,j,k,t)= (  
     &   *( half*(vTrans+ABS(vTrans))
     &      *( Tracer(i,jm1,k,t) + (d0*Rj+d1*Rjm) )
     &  +half*(vTrans-ABS(vTrans))
     &      *( Tracer(i, j ,k,t) - (d0*Rj+d1*Rjp) ))
         ENDDO
        ENDDO
       ENDDO
      ENDDO  
      
      RETURN
      END
