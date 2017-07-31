C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_os7mp_adv_r.F,v 1.7 2008/06/16 13:40:25 jmc Exp $

      SUBROUTINE GAD_OS7MP_ADV_Z(iMax,jMax,Nr,Nt,rA, 
     &              maskC,deltaT,deltaZ,wVel,tracer,wT )
C     /==========================================================\
C     | SUBROUTINE GAD_OS7MP_ADV_Z                               |
C     | o Compute Vertical advective Flux of tracer Q using      |
C     |   7th Order DST Sceheme with monotone preserving limiter |
C     |==========================================================|
      IMPLICIT NONE


      integer*4 i,j,k,t,kp3,kp2,kp1,km1,km2,km3,km4
      integer*4 iMax,jMax,Nr,Nt
      
      real*8 deltaT
      real*8 deltaZ(Nr)
      real*8 recipZ(Nr)
      real*8 rA(iMax,jMax)
      real*8 maskC(iMax,jMax,Nr)
      real*8 maskUp(iMax,jMax,Nr)
      real*8 wVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 wT(iMax,jMax,Nr,Nt)
      
      real*8 cfl,Psi,wTrans
      real*8 wLoc,Fac,DelIp,DelI,Phi,Eps,rp1h,rp1h_cfl
      real*8 recip_DelIp, recip_DelI
      real*8 Qippp,Qipp,Qip,Qi,Qim,Qimm,Qimmm
      real*8 MskIpp,MskIp,MskI,MskIm,MskImm,MskImmm
      real*8 d2,d2p1,d2m1,A,B,C,D
      real*8 dp1h,dm1h, PhiMD,PhiLC,PhiMin,PhiMax
      real*8 DelM,DelP,DelMM,DelPP,DelMMM,DelPPP
      real*8 Del2MM,Del2M,Del2,Del2P,Del2PP
      real*8 Del3MM,Del3M,Del3P,Del3PP
      real*8 Del4M,Del4,Del4P
      real*8 Del5M,Del5P
      real*8 Del6

      Eps = 1.D-20
      
      DO k=1,Nr
       recipZ(k)=1/deltaZ(k)
      ENDDO
      
      DO t=1,Nt    
       DO k=1,Nr
       
        km4=MAX(1,k-4)
        km3=MAX(1,k-3)
        km2=MAX(1,k-2)
        km1=MAX(1,k-1)
        kp1=MIN(Nr,k+1)
        kp2=MIN(Nr,k+2)
        kp3=MIN(Nr,k+3)

         DO j=1,jMax
          DO i=1,iMax
           wTrans= wVel(i,j,k,t)*rA(i,j)
           cfl=ABS( wVel(i,j,k,t)*deltaT*recipZ(k) )

			IF (wTrans.LT.0.D0) THEN
			 Qippp = tracer(i,j,kp2,t)
			 Qipp  = tracer(i,j,kp1,t)
			 Qip   = tracer(i,j,k  ,t)
			 Qi    = tracer(i,j,km1,t)
			 Qim   = tracer(i,j,km2,t)
			 Qimm  = tracer(i,j,km3,t)
			 Qimmm = tracer(i,j,km4,t)

			 MskIpp  = maskC(i,j,kp2) * float(kp2-kp1)
			 MskIp   = maskC(i,j,kp1) * float(kp1-k)
			 MskI    = maskC(i,j,k  ) * float(k-km1)
			 MskIm   = maskC(i,j,km1) * float(km1-km2)
			 MskImm  = maskC(i,j,km2) * float(km2-km3)
			 MskImmm = maskC(i,j,km3) * float(km3-km4)
			ELSEIF (wTrans.GT.0.D0) THEN
			 Qippp = tracer(i,j,km3,t)
			 Qipp  = tracer(i,j,km2,t)
			 Qip   = tracer(i,j,km1,t)
			 Qi    = tracer(i,j,k  ,t)
			 Qim   = tracer(i,j,kp1,t)
			 Qimm  = tracer(i,j,kp2,t)
			 Qimmm = tracer(i,j,kp3,t)

			 MskIpp  = maskC(i,j,km2) * float(km2-km3)
			 MskIp   = maskC(i,j,km1) * float(km1-km2)
			 MskI    = maskC(i,j,k  ) * float(k-km1)
			 MskIm   = maskC(i,j,kp1) * float(kp1-k)
			 MskImm  = maskC(i,j,kp2) * float(kp2-kp1)
			 MskImmm = maskC(i,j,kp3) * float(kp3-kp2)
			ELSE
			 Qippp = 0.D0
			 Qipp  = 0.D0
			 Qip   = 0.D0
			 Qi    = 0.D0
			 Qim   = 0.D0
			 Qimm  = 0.D0
			 Qimmm = 0.D0

			 MskIpp  = 0.D0
			 MskIp   = 0.D0
			 MskI    = 0.D0
			 MskIm   = 0.D0
			 MskImm  = 0.D0
			 MskImmm = 0.D0
			ENDIF

        IF (wTrans.NE.0.D0) THEN
C        2nd order correction [i i-1]
         Fac = 1.D0
         DelP = (Qip-Qi)*MskI
         Phi = Fac * DelP
C        3rd order correction [i i-1 i-2]
         Fac = Fac * ( cfl + 1.D0 )/3.D0
         DelM = (Qi-Qim)*MskIm
         Del2 = DelP - DelM
         Phi = Phi - Fac * Del2
C        4th order correction [i+1 i i-1 i-2]
         Fac = Fac * ( cfl - 2.D0 )/4.D0
         DelPP = (Qipp-Qip)*MskIp*MskI
         Del2P = DelPP - DelP
         Del3P = Del2P - Del2
         Phi = Phi + Fac * Del3p
C        5th order correction [i+1 i i-1 i-2 i-3]
         Fac = Fac * ( cfl - 3.D0 )/5.D0
         DelMM = (Qim-Qimm)*MskImm*MskIm
         Del2M = DelM - DelMM
         Del3M = Del2 - Del2M
         Del4 = Del3P - Del3M
         Phi = Phi + Fac * Del4
C        6th order correction [i+2 i+1 i i-1 i-2 i-3]
         Fac = Fac * ( cfl + 2.D0 )/6.D0
         DelPPP = (Qippp-Qipp)*MskIpp*MskIp*MskI
         Del2PP = DelPP - DelP
         Del3PP = Del2PP - Del2P
         Del4P = Del3PP - Del3P
         Del5P = Del4P - Del4
         Phi = Phi + Fac * Del5P
C        7th order correction [i+2 i+1 i i-1 i-2 i-3 i-4]
         Fac = Fac * ( cfl + 2.D0 )/7.D0
         DelMMM = (Qimm-Qimmm)*MskImmm*MskImm*MskIm
         Del2MM = DelMM - DelMMM
         Del3MM = Del2M - Del2MM
         Del4M = Del3M - Del3MM
         Del5M = Del4 - Del4M
         Del6 = Del5P - Del5M
         Phi = Phi - Fac * Del6

         DelIp = ( Qip - Qi ) * MskI
c        Phi = sign(1.D0,Phi)*sign(1.D0,DelIp)
c    &        *abs(Phi+Eps)/abs(DelIp+Eps)
C--   simplify and avoid division by zero
         recip_DelIp = sign(1.D0,DelIp)/max(abs(DelIp),Eps)
         Phi = Phi*recip_DelIp

         DelI = ( Qi - Qim ) * MskIm
c        rp1h =sign(1.D0,DelI)*sign(1.D0,DelIp)
c    &        *abs(DelI+Eps)/abs(DelIp+Eps)
C--   simplify and avoid division by zero
         recip_DelI = sign(1.D0,DelI)/max(abs(DelI),Eps)
         rp1h = DelI*recip_DelIp
         rp1h_cfl = rp1h/(cfl+Eps)

C        TVD limiter
c        Phi = max(0.D0, min( 2./(1-cfl), Phi, 2.*rp1h_cfl ) )

C        MP limiter
         d2   = Del2 !( ( Qip + Qim ) - 2.*Qi  ) * MskI * MskIm
         d2p1 = Del2P !( ( Qipp + Qi ) - 2.*Qip ) * MskIp * MskI
         d2m1 = Del2M !( ( Qi + Qimm ) - 2.*Qim ) * MskIm * MskImm
         A = 4.D0*d2 - d2p1
         B = 4.D0*d2p1 - d2
         C = d2
         D = d2p1
         dp1h = max(min(A,B,C,D),0.D0)+min(max(A,B,C,D),0.D0)
         A = 4.D0*d2m1 - d2
         B = 4.D0*d2 - d2m1
         C = d2m1
         D = d2
         dm1h = max(min(A,B,C,D),0.D0)+min(max(A,B,C,D),0.D0)
c        qMD = 0.5*( ( Qi + Qip ) - dp1h )
c        qMD = 0.5D0*( ( 2.D0*Qi + DelIp ) - dp1h )
c        qUL = Qi + (1.D0-cfl)/(cfl+Eps)*DelI
c        qLC = Qi + 0.5D0*( 1.D0+dm1h/(DelI+Eps) )*(qUL-Qi)
c        PhiMD = 2.D0/(1.D0-cfl)*(qMD-Qi+Eps)/(DelIp+Eps)
c        PhiLC = 2.D0*rp1h_cfl*(qLC-Qi+Eps)/(qUL-Qi+Eps)
C--   simplify and avoid division by zero
         PhiMD = 1.D0/(1.D0-cfl)*(DelIp-dp1h)*recip_DelIp
         PhiLC = rp1h_cfl*( 1.D0+dm1h*recip_DelI )
C--
         PhiMin = max(min(0.D0,PhiMD),
     &        min(0.D0,2.D0*rp1h_cfl,PhiLC))
         PhiMax = min(max(2.D0/(1.D0-cfl),PhiMD),
     &        max(0.D0,2.D0*rp1h_cfl,PhiLC))
         Phi = max(PhiMin,min(Phi,PhiMax))

         Psi = Phi * 0.5D0 * (1.D0 - cfl)
         wT(i,j,k,t) = wTrans*( Qi + Psi*DelIp )
        ELSE
         wT(i,j,k,t) = 0.D0
        ENDIF

       ENDDO
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END
