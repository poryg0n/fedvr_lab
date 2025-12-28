      MODULE PROPAGATORS
      IMPLICIT NONE
      CONTAINS

       SUBROUTINE SPLIT_OP(ITIME,DT,NX,VEXTV,HAM,EIG,C0,C1)
       IMPLICIT NONE
       INTEGER :: NX ,II, JJ, ITIME
       REAL*8 :: DT, DT2
       REAL*8, DIMENSION(NX) :: XX, EIG
       REAL*8, DIMENSION(NX) :: VEXTV
       COMPLEX*16 :: CE
       COMPLEX*16, DIMENSION(NX) :: C0, C2, C1
       REAL*8, DIMENSION(NX,NX) :: HAM

     
C **** * Split operator ********************************
       DT2=DT/2.D0
       IF (ITIME.eq.0) THEN
          C1=C0
       ELSE
C      exp(-i E dt/2)*Psi
          DO II=1,NX
             CE=DCMPLX(0.d0,eig(II))
             C1(II)=CDEXP(-CE*DT2)*C1(II)
          ENDDO

C      to coordinate representation
          C2=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                C2(II)=C2(II)+HAM(II,JJ)*C1(JJ)
             ENDDO
          ENDDO
C
C      exp(-i V dt)*Psi
          DO II=1,NX
!            Vext=XX(II)+0.5d0   ! external potential
!            CE=DCMPLX(0.d0,Vext)
             CE=DCMPLX(0.d0,VEXTV(II))
             C2(II)=CDEXP(-CE*DT)*C2(II)
          ENDDO

C      back to energy representation
          C1=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                C1(II)=C1(II)+HAM(JJ,II)*C2(JJ)
             ENDDO
          ENDDO
C          
C      exp(-i E dt/2)*Psi
          DO II=1,NX
             CE=DCMPLX(0.d0,eig(II))
             C1(II)=CDEXP(-CE*DT2)*C1(II)
          ENDDO
       END IF
       END SUBROUTINE


       SUBROUTINE IFAB2(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,CP,CN,C3)
       IMPLICIT NONE
       INTEGER :: NX ,II, JJ, ITIME
       REAL*8 :: DT, DT2 
       REAL*8, DIMENSION(NX) :: XX, EIG
       REAL*8, DIMENSION(NX) :: VEXTV, NPSIV
       COMPLEX*16 :: CE
       COMPLEX*16, DIMENSION(NX) :: CP, CPP, CN, CNN, C5, C4, C3 
       REAL*8, DIMENSION(NX,NX) :: HAM

C **** ****** IFAB2 ****************
       DT2=DT/2.0D0
       IF (ITIME.eq.0) THEN
          C3=CP
       ELSE IF (ITIME.eq.1) THEN
          C3=CN
       ELSE 
C **** *  (n>=2) ******
          C5=C3
C      to coordinate representation
          C4=DCMPLX(0.d0,0.d0)
          CPP=DCMPLX(0.d0,0.d0)
          CNN=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CPP(II)=CPP(II)+HAM(II,JJ)*CP(JJ)
                CNN(II)=CNN(II)+HAM(II,JJ)*CN(JJ)
             ENDDO
          ENDDO

C      N(Psi,t) = n*Psi
          DO II=1,NX
             CPP(II)=NPSIV(II)*CPP(II)
             CNN(II)=NPSIV(II)*CNN(II)
          ENDDO

C      back to energy representation
          CP=DCMPLX(0.d0,0.d0)
          CN=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+HAM(JJ,II)*CPP(JJ)
                CN(II)=CN(II)+HAM(JJ,II)*CNN(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt/2)*Psi
          DO II=1,NX
             CE=DCMPLX(0.d0,eig(II))
             CP(II)=CDEXP(-CE*DT)*CP(II)
             CN(II)=CDEXP(-CE*DT2)*CN(II)
             C3(II)=CDEXP(-CE*DT2)*C3(II)
          ENDDO

C      to coordinate representation
          CPP=DCMPLX(0.d0,0.d0)
          CNN=DCMPLX(0.d0,0.d0)
          C4=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CPP(II)=CPP(II)+HAM(II,JJ)*CP(JJ)
                CNN(II)=CNN(II)+HAM(II,JJ)*CN(JJ)
                C4(II)=C4(II)+HAM(II,JJ)*C3(JJ)
             ENDDO
          ENDDO

C      exp(-i dt Vext) * Psi
             DO II=1,NX
                CE=DCMPLX(0.d0,VEXTV(II))
                CPP(II)=CDEXP(-2*CE*DT)*CPP(II)
                CNN(II)=CDEXP(-CE*DT)*CNN(II)
                C4(II)=CDEXP(-CE*DT)*C4(II)
          ENDDO

C      back to energy representation
          CP=DCMPLX(0.d0,0.d0)
          CN=DCMPLX(0.d0,0.d0)
          C3=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+HAM(JJ,II)*CPP(JJ)
                CN(II)=CN(II)+HAM(JJ,II)*CNN(JJ)
                C3(II)=C3(II)+HAM(JJ,II)*C4(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt)*Psi
          DO II=1,NX
             CE=DCMPLX(0.d0,eig(II))
             CP(II)=DCMPLX(0.0d0,DT2)*CDEXP(-CE*DT)*CP(II)
             CN(II)=DCMPLX(0.0d0,-3*DT2)*CDEXP(-CE*DT2)*CN(II)
             C3(II)=CDEXP(-CE*DT2)*C3(II)
             C3(II)=C3(II)+CN(II)+CP(II)
          ENDDO

          CP=C5
          CN=C3
       END IF
       END SUBROUTINE



       SUBROUTINE IFRK1(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,C0,C3)
       IMPLICIT NONE
       INTEGER :: NX ,II, JJ, ITIME
       REAL*8 :: DT, DT2, VEXT, NPSI
       REAL*8, DIMENSION(NX), INTENT(IN) :: EIG, VEXTV
       COMPLEX*16 :: CE
       COMPLEX*16, DIMENSION(NX), INTENT(IN) :: C0
       COMPLEX*16, DIMENSION(NX) :: CN, CNN, C4, C3 
       COMPLEX*16, DIMENSION(NX) :: CK1
       REAL*8, DIMENSION(NX,NX), INTENT(IN) :: HAM, NPSIV

       DT2=DT/2.0D0
       IF (ITIME.eq.0) THEN
          C3=C0
       ELSE 
          CALL IFRK_S1(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,0,0.D0,CK1,C3)
          CALL START_VEC(ITIME,DT,NX,VEXTV,HAM,EIG,0.0D0,0,CN,C3)

          C3=CN+DT*CK1

       END IF
       END SUBROUTINE


       SUBROUTINE IFRK2(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,C0,C3)
       IMPLICIT NONE
       INTEGER :: NX ,II, JJ, ITIME
       REAL*8 :: DT, DT2, BETA, ALPHA
       REAL*8, DIMENSION(NX), INTENT(IN) :: EIG, VEXTV
       COMPLEX*16 :: CE
       COMPLEX*16, DIMENSION(NX), INTENT(IN) :: C0
       COMPLEX*16, DIMENSION(NX) :: CP, CN, CPP, CNN, C3
       COMPLEX*16, DIMENSION(NX) :: CK1, CK2
       REAL*8, DIMENSION(NX,NX), INTENT(IN) :: HAM, NPSIV


       ALPHA =  1.D0/1
       BETA =  1.D0/(2.D0*ALPHA)
       DT2=DT/2.0D0
       IF (ITIME.eq.0) THEN
          C3=C0
       ELSE 
C **** *  (n>=2) ******
! **** ******* Compute K1(t_n) ****************
          CALL IFRK_S1(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,0,0.D0,CK1,C3)
          CALL IFRK_S2(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,0,ALPHA,1.D0,
     &  CK2,C3)

          CALL START_VEC(ITIME,DT,NX,VEXTV,HAM,EIG,0.0D0,0,CN,C3)

          C3 = CN + DT*((1-BETA)*CK1 + BETA*CK2)


       END IF
       END SUBROUTINE


       SUBROUTINE IFRK3(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,C0,C3)
       IMPLICIT NONE
       INTEGER :: NX ,II, JJ, ITIME
       REAL*8 :: DT, DT2, ALPHA, BETA, DZETA
       REAL*8, DIMENSION(NX), INTENT(IN) :: EIG, VEXTV
       COMPLEX*16 :: CE
       COMPLEX*16, DIMENSION(NX), INTENT(IN) :: C0
       COMPLEX*16, DIMENSION(NX) :: CP, CN, CPP, CNN, C3
       COMPLEX*16, DIMENSION(NX) :: CK1, CK2, CK3
       REAL*8, DIMENSION(NX,NX), INTENT(IN) :: HAM, NPSIV

C      Generic
       ALPHA =  1.D0/2             
       BETA = 1.D0
       DZETA = 1.D0

C      Kutta3
C      ALPHA =  1.D0/2
C      BETA =  2*ALPHA
C      DZETA=1.D0

C      Heun3
C      ALPHA =  1.D0/3
C      BETA =  2*ALPHA
C      DZETA=1.D0

C      Ralston3
C      ALPHA=1.D0/2
C      BETA=3.D0/4
C      DZETA=1.D0



       DT2=DT/2.0D0
       IF (ITIME.eq.0) THEN
          C3=C0
       ELSE 
C **** *  (n>=2) ******
! **** ******* Compute K1(t_n) ****************
          CALL IFRK_S1(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,0,0.D0,CK1,C3)
          CALL IFRK_S2(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,0,ALPHA,BETA,
     &  CK2,C3)
C         CALL IFRK_S3(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,0,
C    &  ALPHA,BETA,DZETA,CK3,C3)
          CALL IFRK_S3G(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,0,
     &  ALPHA,BETA,DZETA,CK3,C3)

          CALL START_VEC(ITIME,DT,NX,VEXTV,HAM,EIG,0.0D0,0,CN,C3)

!         C3 = CN + DT*((1.D0/4)*CK1 + (3.D0/4)*CK3)   ! Heun3
!         C3 = CN + DT*((2.D0/9)*CK1 + (1.D0/3)*CK2 + (4.D0/9)*CK3) ! Ralsont3
!         C3 = CN + DT*((1.D0/6)*CK1 + (2.D0/3)*CK2 + (1.D0/6)*CK3) ! Kutta3

!         If Generic (IFRK_S3G), use this one
          C3 = CN + DT*((.5D0-1.D0/(6*ALPHA))*CK1 
     &  + (1.D0/(6*ALPHA*(1.D0-ALPHA)))*CK2
     &  + ((2.D0-3*ALPHA)/(6*(1.D0-ALPHA)))*CK3)


       END IF
       END SUBROUTINE


       SUBROUTINE IFRK4(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,C0,C3)
       IMPLICIT NONE
       INTEGER :: NX ,II, JJ, ITIME
       REAL*8 :: DT, DT2, ALPHA, BETA, DZETA
       REAL*8, DIMENSION(NX), INTENT(IN) :: EIG, VEXTV
       COMPLEX*16 :: CE
       COMPLEX*16, DIMENSION(NX), INTENT(IN) :: C0
       COMPLEX*16, DIMENSION(NX) :: CP, CN, CPP, CNN, C3
       COMPLEX*16, DIMENSION(NX) :: CK1, CK2, CK3, CK4
       REAL*8, DIMENSION(NX,NX), INTENT(IN) :: HAM, NPSIV

!      Classic
       ALPHA =  1.D0/2             
       BETA = ALPHA
       DZETA = 1.D0


!      38_rule
!      ALPHA =  1.D0/3            
!      BETA = 2*ALPHA
!      DZETA = 1.D0


       DT2=DT/2.0D0
       IF (ITIME.eq.0) THEN
          C3=C0
       ELSE 
C **** *  (n>=2) ******
! **** ******* Compute K1(t_n) ****************
          CALL IFRK_S1(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,0,0.D0,CK1,C3)
          CALL IFRK_S2(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,0,ALPHA,BETA,
     &  CK2,C3)
          CALL IFRK_S3(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,0,
     &  ALPHA,BETA,DZETA,CK3,C3)
          CALL IFRK_S4(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,0,
     &  ALPHA,BETA,DZETA,CK4,C3)

          CALL START_VEC(ITIME,DT,NX,VEXTV,HAM,EIG,0.0D0,0,CN,C3)

          C3 = CN + (DT/6)*(CK1 + 2.D0*CK2 + 2.D0*CK3 + CK4)    ! Classic_RK4
!         C3 = CN + (DT/8)*(CK1 + 3.D0*CK2 + 3.D0*CK3 + CK4)    ! 38_RK4


       END IF
       END SUBROUTINE







       SUBROUTINE IFRK_S1(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,MIDS,ALPHA,
     &  CK1,C0)
       IMPLICIT NONE
       INTEGER :: NX, II, JJ, ITIME, MIDS
       REAL*8 :: DT, DT2, ALPHA, GAMMA
       COMPLEX*16 :: CE
       REAL*8, DIMENSION(NX), INTENT(IN) :: EIG, VEXTV
       COMPLEX*16, DIMENSION(NX) :: CP, CPP, C0,
     &                                CK1
       REAL*8, DIMENSION(NX,NX), INTENT(IN) :: HAM, NPSIV

       IF (MIDS.EQ.0) THEN
          GAMMA = (1.D0-ALPHA)
       ELSE
          GAMMA = ALPHA
       END IF

       DT2=DT/2.0D0

          CP=C0

C      to coordinate representation
          CPP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CPP(II)=CPP(II)+HAM(II,JJ)*CP(JJ)
             ENDDO
          ENDDO

C      N(Psi,t) = n*Psi
!         CPP=NPSIV*CPP
!         CPP=MATMUL(NPSIV,CPP)
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+NPSIV(II,JJ)*CPP(JJ)
             ENDDO
          ENDDO
          CPP=CP


C      back to energy representation
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+HAM(JJ,II)*CPP(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt/2)*Psi
          DO II=1,NX
!            CE=(1.D0-ALPHA)*DCMPLX(0.d0,eig(II))
             CE=GAMMA*DCMPLX(0.d0,eig(II))
             CP(II)=CDEXP(-CE*DT2)*CP(II)
          ENDDO

C      to coordinate representation
          CPP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CPP(II)=CPP(II)+HAM(II,JJ)*CP(JJ)
             ENDDO
          ENDDO

C      exp(-i dt Vext) * Psi
          DO II=1,NX
!            CE=(1.D0-ALPHA)*DCMPLX(0.d0,VEXTV(II))
             CE=GAMMA*DCMPLX(0.d0,VEXTV(II))
             CPP(II)=CDEXP(-CE*DT)*CPP(II)
          ENDDO

C      back to energy representation
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+HAM(JJ,II)*CPP(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt)*Psi
          DO II=1,NX
!            CE=(1.D0-ALPHA)*DCMPLX(0.d0,eig(II))
             CE=GAMMA*DCMPLX(0.d0,eig(II))
             CK1(II)=DCMPLX(0.0d0,-1.D0)*CDEXP(-CE*DT2)*CP(II)
          ENDDO

       END SUBROUTINE



       SUBROUTINE IFRK_S2(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,MIDS,
     &  ALPHA,BETA,CK2,C0)
       IMPLICIT NONE
       INTEGER :: NX, II, JJ, ITIME, MIDS
       REAL*8 :: DT, DT2, ALPHA, BETA, GAMMA
       COMPLEX*16 :: CE
       REAL*8, DIMENSION(NX), INTENT(IN) :: EIG, VEXTV
       COMPLEX*16, DIMENSION(NX) :: CP, CPP, CN, CNN, C0,
     &                               CK1, CK2
       REAL*8, DIMENSION(NX,NX), INTENT(IN) :: HAM, NPSIV

       IF (MIDS.EQ.0) THEN
          GAMMA = (1.D0-ALPHA)
       ELSE
          GAMMA = (BETA-ALPHA)
       END IF

       DT2 = DT/2.D0
! **** ******* Compute K1(t_n + 1/2) ****************
C      to coordinate representation

          CALL IFRK_S1(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,1,ALPHA,CK1,C0)
          CALL START_VEC(ITIME,DT,NX,VEXTV,HAM,EIG,ALPHA,1,CN,C0)

          CP = CN + ALPHA*DT*CK1

! **** ******* Compute K2(t_n) ****************
C      to coordinate representation
          CPP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CPP(II)=CPP(II)+HAM(II,JJ)*CP(JJ)
             ENDDO
          ENDDO

C      N(Psi,t) = n*Psi
!         CPP=NPSIV*CPP
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+NPSIV(II,JJ)*CPP(JJ)
             ENDDO
          ENDDO
          CPP=CP

C      back to energy representation
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+HAM(JJ,II)*CPP(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt/2)*Psi
          DO II=1,NX
!            CE=DCMPLX(0.d0,eig(II))
!            CP(II)=CDEXP(-CE*(1.D0-ALPHA)*DT2)*CP(II)
             CE=GAMMA*DCMPLX(0.d0,eig(II))
             CP(II)=CDEXP(-CE*DT2)*CP(II)
          ENDDO

C      to coordinate representation
          CPP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CPP(II)=CPP(II)+HAM(II,JJ)*CP(JJ)
             ENDDO
          ENDDO

C      exp(-i dt Vext) * Psi
          DO II=1,NX
C            CE=DCMPLX(0.d0,VEXTV(II))
             CE=GAMMA*DCMPLX(0.d0,VEXTV(II))
             CPP(II)=CDEXP(-CE*DT)*CPP(II)
          ENDDO

C      back to energy representation
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+HAM(JJ,II)*CPP(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt)*Psi
          DO II=1,NX
!            CE=DCMPLX(0.d0,eig(II))
!            CK2(II)=DCMPLX(0.0d0,-1.D0)*CDEXP(-CE*(1-ALPHA)*DT2)*CP(II)
             CE=GAMMA*DCMPLX(0.d0,eig(II))
             CK2(II)=DCMPLX(0.0d0,-1.D0)*CDEXP(-CE*DT2)*CP(II)
          ENDDO


       END SUBROUTINE


       SUBROUTINE IFRK_S3(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,MIDS,
     &  ALPHA,BETA,DZETA,CK3,C0)
       IMPLICIT NONE
       INTEGER :: NX ,II, JJ, ITIME, MIDS
       REAL*8 :: DT, DT2, BETA, ALPHA, DZETA, GAMMA
       REAL*8, DIMENSION(NX), INTENT(IN) :: EIG, VEXTV
       COMPLEX*16 :: CE
       COMPLEX*16, DIMENSION(NX), INTENT(IN) :: C0
       COMPLEX*16, DIMENSION(NX) :: CP, CN, CPP, CNN, C3
       COMPLEX*16, DIMENSION(NX) :: CK1, CK2, CK3
       COMPLEX*16, DIMENSION(NX) :: CKP, CKPP
       REAL*8, DIMENSION(NX,NX), INTENT(IN) :: HAM, NPSIV


       IF (MIDS.EQ.0) THEN
          GAMMA = (1.D0-BETA)
       ELSE
          GAMMA = (DZETA-BETA)
       END IF

       DT2=DT/2.0D0
C **** *  (n>=2) ******
! **** ******* Compute K1(t_n) ****************
          CALL IFRK_S1(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,1,BETA,CK1,C0)
          CALL IFRK_S2(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,1,ALPHA,BETA,
     &  CK2,C0)
          CALL START_VEC(ITIME,DT,NX,VEXTV,HAM,EIG,BETA,1,CN,C0)

          CP = CN + BETA*DT*CK2                 ! Heun3, Ralston3 or RK4
!         CP = CN + DT*(-1.D0*CK1 + 2.D0*CK2)   ! Kutta3
!         CP = CN + DT*((-1.D0/3)*CK1 + 1.D0*CK2)   ! 38_RK4


! **** ******* Compute K3(t_n) ****************
C      to coordinate representation
          CPP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CPP(II)=CPP(II)+HAM(II,JJ)*CP(JJ)
             ENDDO
          ENDDO

C      N(Psi,t) = n*Psi
!         CPP=NPSIV*CPP
!         CPP=MATMUL(NPSIV,CPP)
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+NPSIV(II,JJ)*CPP(JJ)
             ENDDO
          ENDDO
          CPP=CP


C      back to energy representation
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+HAM(JJ,II)*CPP(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt/2)*Psi
          DO II=1,NX
             CE=GAMMA*DCMPLX(0.d0,eig(II))
             CP(II)=CDEXP(-CE*DT2)*CP(II)
          ENDDO

C      to coordinate representation
          CPP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CPP(II)=CPP(II)+HAM(II,JJ)*CP(JJ)
             ENDDO
          ENDDO

C      exp(-i dt Vext) * Psi
          DO II=1,NX
             CE=GAMMA*DCMPLX(0.d0,VEXTV(II))
             CPP(II)=CDEXP(-CE*DT)*CPP(II)
          ENDDO

C      back to energy representation
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+HAM(JJ,II)*CPP(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt)*Psi
          DO II=1,NX
             CE=GAMMA*DCMPLX(0.d0,eig(II))
             CK3(II)=DCMPLX(0.0d0,-1.D0)*CDEXP(-CE*DT2)*CP(II)
          ENDDO


       END SUBROUTINE



       SUBROUTINE IFRK_S3G(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,MIDS,
     &  ALPHA,BETA,DZETA,CK3,C0)
       IMPLICIT NONE
       INTEGER :: NX ,II, JJ, ITIME, MIDS
       REAL*8 :: DT, DT2, BETA, ALPHA, DZETA, GAMMA
       REAL*8, DIMENSION(NX), INTENT(IN) :: EIG, VEXTV
       COMPLEX*16 :: CE
       COMPLEX*16, DIMENSION(NX), INTENT(IN) :: C0
       COMPLEX*16, DIMENSION(NX) :: CP, CN, CPP, CNN, C3
       COMPLEX*16, DIMENSION(NX) :: CK1, CK2, CK3
       COMPLEX*16, DIMENSION(NX) :: CKP, CKPP
       REAL*8, DIMENSION(NX,NX), INTENT(IN) :: HAM, NPSIV


       IF (MIDS.EQ.0) THEN
          GAMMA = (1.D0-BETA)
       ELSE
          GAMMA = (DZETA-BETA)
       END IF

       DT2=DT/2.0D0
C **** *  (n>=2) ******
! **** ******* Compute K1(t_n) ****************
          CALL IFRK_S1(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,1,BETA,CK1,C0)
          CALL IFRK_S2(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,1,ALPHA,BETA,
     &  CK2,C0)
          CALL START_VEC(ITIME,DT,NX,VEXTV,HAM,EIG,BETA,1,CN,C0)

          CP = CN + DT*((1.D0+(1.D0-ALPHA)/(ALPHA*(3*ALPHA-2.D0)))*CK1 
     &  - ((1.D0-ALPHA)/(ALPHA*(3*ALPHA-2.D0)))*CK2)   ! Generic


! **** ******* Compute K3(t_n) ****************
C      to coordinate representation
          CPP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CPP(II)=CPP(II)+HAM(II,JJ)*CP(JJ)
             ENDDO
          ENDDO

C      N(Psi,t) = n*Psi
!         CPP=NPSIV*CPP
!         CPP=MATMUL(NPSIV,CPP)
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+NPSIV(II,JJ)*CPP(JJ)
             ENDDO
          ENDDO
          CPP=CP


C      back to energy representation
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+HAM(JJ,II)*CPP(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt/2)*Psi
          DO II=1,NX
             CE=GAMMA*DCMPLX(0.d0,eig(II))
             CP(II)=CDEXP(-CE*DT2)*CP(II)
          ENDDO

C      to coordinate representation
          CPP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CPP(II)=CPP(II)+HAM(II,JJ)*CP(JJ)
             ENDDO
          ENDDO

C      exp(-i dt Vext) * Psi
          DO II=1,NX
             CE=GAMMA*DCMPLX(0.d0,VEXTV(II))
             CPP(II)=CDEXP(-CE*DT)*CPP(II)
          ENDDO

C      back to energy representation
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+HAM(JJ,II)*CPP(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt)*Psi
          DO II=1,NX
             CE=GAMMA*DCMPLX(0.d0,eig(II))
             CK3(II)=DCMPLX(0.0d0,-1.D0)*CDEXP(-CE*DT2)*CP(II)
          ENDDO


       END SUBROUTINE



       SUBROUTINE IFRK_S4(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,MIDS,
     &  ALPHA,BETA,DZETA,CK4,C0)
       IMPLICIT NONE
       INTEGER :: NX ,II, JJ, ITIME, MIDS
       REAL*8 :: DT, DT2, BETA, ALPHA, DZETA, GAMMA
       REAL*8, DIMENSION(NX), INTENT(IN) :: EIG, VEXTV
       COMPLEX*16 :: CE
       COMPLEX*16, DIMENSION(NX), INTENT(IN) :: C0
       COMPLEX*16, DIMENSION(NX) :: CP, CN, CPP, CNN, C3
       COMPLEX*16, DIMENSION(NX) :: CK1, CK2, CK3, CK4
       COMPLEX*16, DIMENSION(NX) :: CKP, CKPP
       REAL*8, DIMENSION(NX,NX), INTENT(IN) :: HAM, NPSIV


       IF (MIDS.EQ.0) THEN
          GAMMA = (1.D0-DZETA)
       ELSE
          GAMMA = (1.D0-DZETA)
       END IF

       DT2=DT/2.0D0
C **** *  (n>=2) ******
! **** ******* Compute K1(t_n) ****************
C         CALL IFRK_S1(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,1,DZETA,CK1,C0)
C         CALL IFRK_S2(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,1,ALPHA,BETA,
C    &  CK2,C0)

          CALL IFRK_S3(ITIME,DT,NX,NPSIV,VEXTV,HAM,EIG,1,
     &  ALPHA,BETA,DZETA,CK3,C0)
!         CALL START_VEC(ITIME,DT,NX,VEXTV,HAM,EIG,0.D0,0,CN,C0)
          CALL START_VEC(ITIME,DT,NX,VEXTV,HAM,EIG,DZETA,1,CN,C0)

          CP = CN + DZETA*DT*CK3            
!         CP = CN + DZETA*DT*(CK1 - CK2 + CK3)            ! 38_RK4


! **** ******* Compute K4(t_n) ****************
C      to coordinate representation
          CPP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CPP(II)=CPP(II)+HAM(II,JJ)*CP(JJ)
             ENDDO
          ENDDO

C      N(Psi,t) = n*Psi
!         CPP=NPSIV*CPP
!         CPP=MATMUL(NPSIV,CPP)
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+NPSIV(II,JJ)*CPP(JJ)
             ENDDO
          ENDDO
          CPP=CP


C      back to energy representation
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+HAM(JJ,II)*CPP(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt/2)*Psi
          DO II=1,NX
             CE=GAMMA*DCMPLX(0.d0,eig(II))
             CP(II)=CDEXP(-CE*DT2)*CP(II)
          ENDDO

C      to coordinate representation
          CPP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CPP(II)=CPP(II)+HAM(II,JJ)*CP(JJ)
             ENDDO
          ENDDO

C      exp(-i dt Vext) * Psi
          DO II=1,NX
             CE=GAMMA*DCMPLX(0.d0,VEXTV(II))
             CPP(II)=CDEXP(-CE*DT)*CPP(II)
          ENDDO

C      back to energy representation
          CP=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CP(II)=CP(II)+HAM(JJ,II)*CPP(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt)*Psi
          DO II=1,NX
             CE=GAMMA*DCMPLX(0.d0,eig(II))
             CK4(II)=DCMPLX(0.0d0,-1.D0)*CDEXP(-CE*DT2)*CP(II)
          ENDDO


       END SUBROUTINE



       SUBROUTINE START_VEC(ITIME,DT,NX,VEXTV,HAM,EIG,ALPHA,MIDS,
     &  CN,C0)
       IMPLICIT NONE
       INTEGER :: NX, II, JJ, ITIME, MIDS
       REAL*8 :: DT, DT2, ALPHA, GAMMA
       REAL*8, DIMENSION(NX), INTENT(IN) :: EIG, VEXTV
       COMPLEX*16 :: CE
       COMPLEX*16, DIMENSION(NX) :: CNN, CN, C0
       REAL*8, DIMENSION(NX,NX), INTENT(IN) :: HAM

       IF (MIDS.EQ.0) THEN
          GAMMA = (1.D0-ALPHA)
       ELSE
          GAMMA = ALPHA
       END IF

       DT2=DT/2.0D0
C      exp(-i E dt/2)*Psi
          DO II=1,NX
!            CE=(1.D0-ALPHA)*DCMPLX(0.d0,eig(II))
             CE=GAMMA*DCMPLX(0.d0,eig(II))
             CN(II)=CDEXP(-CE*DT2)*C0(II)
          ENDDO

CC      to coordinate representation
          CNN=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CNN(II)=CNN(II)+HAM(II,JJ)*CN(JJ)
             ENDDO
          ENDDO

C      exp(-i dt Vext) * Psi
          DO II=1,NX
!            CE=(1.D0-ALPHA)*DCMPLX(0.d0,VEXTV(II))
             CE=GAMMA*DCMPLX(0.d0,VEXTV(II))
             CNN(II)=CDEXP(-CE*DT)*CNN(II)
          ENDDO

C      back to energy representation
          CN=DCMPLX(0.d0,0.d0)
          DO II=1,NX
             DO JJ=1,NX
                CN(II)=CN(II)+HAM(JJ,II)*CNN(JJ)
             ENDDO
          ENDDO

C      exp(-i E dt)*Psi
          DO II=1,NX
!            CE=(1.D0-ALPHA)*DCMPLX(0.d0,eig(II))
             CE=GAMMA*DCMPLX(0.d0,eig(II))
             CN(II)=CDEXP(-CE*DT2)*CN(II)
          ENDDO
       END SUBROUTINE
      END MODULE

