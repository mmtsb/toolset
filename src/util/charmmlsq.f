      SUBROUTINE FROTU(R,U)
C-----------------------------------------------------------------------
C     THIS ROUTINE SOLVES THE CONSTRAINED MINIMIZATION EQUATION
C     USING LAGRANGE MULTIPLIERS.
C     BERNARD R. BROOKS
C
      implicit none

      REAL*8 R(3,3),W(6),A(3,3),B(3,3),U(3,3),SCR(24)
      INTEGER I,I1,I2,IPT,J,K,JP,JQ,KP,KQ
      REAL*8 TRACE,EVS,DET
C

      DET=0.0
      DO I=1,3
        I1=I+1
        IF(I1.GT.3) I1=I1-3
        I2=I+2
        IF(I2.GT.3) I2=I2-3
        DET=DET+R(I,1)*(R(I1,2)*R(I2,3)-R(I2,2)*R(I1,3))
      ENDDO
C
      IPT=0
      DO I=1,3
         DO J=I,3
            IPT=IPT+1
            W(IPT)=0.0
            DO K=1,3
               W(IPT)=W(IPT)+R(J,K)*R(I,K)
            ENDDO
         ENDDO
      ENDDO
C
      TRACE=W(1)+W(4)+W(6)

      IF(TRACE.LT.3.D-6) THEN
        DO I=1,3
          DO J=1,3
            U(I,J)=0.D0
          ENDDO
          U(I,I)=1.D0
        ENDDO
        RETURN
      ENDIF
C
      CALL DIAGQ(3,3,W,A,SCR(4),SCR(7),SCR(10),SCR(13),SCR(1),
     $           SCR(16),SCR(19),SCR(22),0)
C
      DO I=1,3
        SCR(I)=SQRT(ABS(SCR(I)))
        IF(SCR(I).LT.1.D-6) SCR(I)=1.D-6
      ENDDO
C
      IF(DET.LT.0.D0) SCR(1)=-SCR(1)
C
      DO J=1,3
        EVS=SCR(J)
        DO I=1,3
           B(I,J)=0.0
           DO K=1,3
              B(I,J)=B(I,J)+R(K,I)*A(K,J)/EVS
           ENDDO
         ENDDO
      ENDDO
C
      DET=0.0
      DO I=1,3
        I1=I+1
        IF(I1.GT.3) I1=I1-3
        I2=I+2
        IF(I2.GT.3) I2=I2-3
        DET=DET+A(I,1)*(A(I1,2)*A(I2,3)-A(I2,2)*A(I1,3))
      ENDDO

      DO J=1,3
        IF(ABS(SCR(J)).LE.1.D-6) THEN
          JP=J+1
          JQ=J+2
          IF(JP.GT.3) JP=JP-3
          IF(JQ.GT.3) JQ=JQ-3
          DO K=1,3
            KP=K+1
            KQ=K+2
            IF(KP.GT.3) KP=KP-3
            IF(KQ.GT.3) KQ=KQ-3
            B(K,J)=B(KP,JP)*B(KQ,JQ)-B(KP,JQ)*B(KQ,JP)
            IF(DET.LT.0.D0) B(K,J)=-B(K,J)
          ENDDO
        ENDIF
        CALL NORMAL(B(1,J),3)
      ENDDO
C
      DO 30 J=1,3
        DO 30 I=1,3
          U(I,J)=0.D0
          DO 30 K=1,3
            U(I,J)=U(I,J)+A(I,K)*B(J,K)
 30   CONTINUE
C
      DO J=1,3
        CALL NORMAL(U(1,J),3)
      ENDDO
C
C     CHECK TO INSURE UNITY (AS OPPOSED TO ANTI-UNITARY)
      DET=0.0
      DO I=1,3
        I1=I+1
        IF(I1.GT.3) I1=I1-3
        I2=I+2
        IF(I2.GT.3) I2=I2-3
        DET=DET+U(I,1)*(U(I1,2)*U(I2,3)-U(I2,2)*U(I1,3))
      ENDDO
      
c      IF(ABS(DET-1.0).GT.1.D-4) WRITE(*,55) DET

 55   FORMAT(/' ***** WARNING ***** FROM FROTU. ROTATION MATRIX IS ',
     $       'NOT UNITARY.'/,'  THE DETERMINANT IS',F14.8/)
C
      RETURN
      END


      SUBROUTINE DIAGQ(NX,NFRQX,DD,VEC,A,B,P,W,EV,TA,TB,Y,NADD)
C
C   THIS ROUTINE IS A CONGLOMERATION OF GIVEN, HOUSEC, AND EIGEN
C   WHERE THE BEST FEATURES OF EACH WERE KEPT AND SEVERAL OTHER
C   MODIFICATIONS HAVE BEEN MADE TO INCREASE EFFICIENCY AND ACCURACY.
C
C   By Bernard R. Brooks   1981
C
C   NX      - ORDER OF MATRIX
C   NFRQX   - NUMBER OF ROOTS DESIRED
C   DD      - SECOND DERIVATIVE MATRIX IN UPPER TRIANGULAR FORM
C   VEC     - EIGENVECTORS RETURNED (NX,NFRQX)
C   EV      - EIGENVALUES RETURNED (NX)
C   A,B,P,W,TA,TB,Y - ALL SCRATCH VECTORS (NX+1)
C   NADD    - NUMBER OF LOWEST ROOTS TO SKIP
C
C
C
      implicit none
      INTEGER NX,NFRQX,NADD
      REAL*8 DD(*),VEC(*),A(*),B(*),P(*),W(*),EV(*),TA(*),TB(*),Y(*)
C
      REAL*8 ETA,THETA,DEL1,DELTA,SMALL,DELBIG,THETA1,TOLER
      REAL*8 RPOWER,RPOW1,RAND1,DUNITY,FACTOR,ANORM,U,ANORMR
      REAL*8 SUM1,BX,S,SGN,TEMP,XKAP,EXPR,ALIMIT,ROOTL,ROOTX,TRIAL,F0
      REAL*8 AROOT,ELIM1,ELIM2,T,EPR,XNORM,XNORM1
      INTEGER N,NEV,NEVADD,NTOT,I,IPT,J,IJ,NN,MI,MI1,JI,JI2,II
      INTEGER ML,ML1,L,M,K,MJ,MJ1,NOMTCH,NOM,IA,ITER
      INTEGER J1,MK,MK1,KK
      REAL*8 vt,rl,ei
C
C
      ETA=2.22045D-16
      THETA=4.4923D+307
C
      N=NX
      NEV=NFRQX
      NEVADD=NEV+NADD
C
      DEL1=ETA/100.0
      DELTA=ETA**2*100.0
      SMALL=ETA**2/100.0
      DELBIG=THETA*DELTA/1000.0
      THETA1=1000.0/THETA
      TOLER=100.0*ETA
      RPOWER=8388608.0
      RPOW1=RPOWER*0.50
      RAND1=RPOWER-3.0
      DUNITY=1.0

C
      FACTOR=0.0
      NTOT=(N*(N+1))/2
      DO I=1,NTOT
         FACTOR=MAX(FACTOR,ABS(DD(I)))
      ENDDO
      IF(FACTOR.GT.THETA1) GOTO 818
      WRITE(*,811)
 811  FORMAT(' WARNING FROM <DIAGQ>. ZERO MATRIX PASSED.',
     1     ' IDENTITY MATRIX RETURNED.')
      DO I=1,NEV
         EV(I)=0.0
         IPT=(I-1)*N
         DO J=1,N
            IPT=IPT+1
            VEC(IPT)=0.0
            IF(I+NADD.EQ.J) VEC(IPT)=1.0
         ENDDO
      ENDDO
      RETURN
C
 818  CONTINUE

      IJ=0
      ANORM=0.0
      DO I=1,N
         DO J=I,N
            IJ=IJ+1
            U=(DD(IJ)/FACTOR)**2
            IF(I.EQ.J) U=U*0.5
            ANORM=ANORM+U
         ENDDO
      ENDDO
      ANORM=SQRT(ANORM+ANORM)*FACTOR
      ANORMR=DUNITY/ANORM
      DO I=1,NTOT
         DD(I)=DD(I)*ANORMR
      ENDDO

C
      NN=N-1
      MI=0
      MI1=N-1
C *
C * LOOP THRU 70.
C *
      DO 70 I=1,NN
         SUM1=0.0
         B(I)=0.0
         JI=I+1
         IPT=MI+I
         A(I)=DD(IPT)
         IPT=IPT+1
         BX=DD(IPT)
         JI2=JI+1
         IF(JI.EQ.N) GOTO 11
         DO J=JI2,N
            IPT=IPT+1
            SUM1=SUM1+DD(IPT)*DD(IPT)
         ENDDO
         IF(SUM1.GT.SMALL) GOTO 15
C
 11      B(I)=BX
         DD(MI+JI)=0.0
         MI=MI+MI1
         MI1=MI1-1
         GOTO 70
C
 15      CONTINUE
         
         S=SQRT(SUM1+BX**2)
         SGN=SIGN(DUNITY,BX)
         TEMP=ABS(BX)
         W(JI)=SQRT(.5*(DUNITY+(TEMP/S)))
         IPT=MI+JI
         DD(IPT)=W(JI)
         II=I+2
         IF(II.GT.N) GOTO 22
         TEMP=SGN/(2.0*W(JI)*S)
         DO J=II,N
            IPT=IPT+1
            W(J)=TEMP * DD(IPT)
            DD(IPT)=W(J)
         ENDDO
 22      B(I)=-SGN*S
C
C
         DO J=JI,N
            P(J)=0.0
         ENDDO
         ML=MI+MI1
         ML1=MI1-1
         DO L=JI,N
            IPT=ML+L
            DO M=L,N
               BX=DD(IPT)
               P(L)=P(L)+BX*W(M)
               IF(L.NE.M) P(M)=P(M)+BX*W(L)
               IPT=IPT+1
            ENDDO
            ML=ML+ML1
            ML1=ML1-1
         ENDDO
C
C
         XKAP=0.0
         DO K=JI,N
            XKAP=XKAP+W(K)*P(K)
         ENDDO
         DO L=JI,N
            P(L)=P(L)-XKAP*W(L)
         ENDDO
         MJ=MI+MI1
         MJ1=MI1-1
         DO J=JI,N
            DO K=J,N
               EXPR=(P(J)*W(K))+(P(K)*W(J))
               DD(MJ+K)=DD(MJ+K)-EXPR-EXPR
            ENDDO
            MJ=MJ+MJ1
            MJ1=MJ1-1
         ENDDO
         MI=MI+MI1
         MI1=MI1-1
 70   CONTINUE

C *
C * END OF FIRST MAJOR PART. NEXT BEGIN STURM METHOD.
C *

      A(N)=DD(MI+N)
      B(N)=0.0
C
      ALIMIT=1.0
      DO I=1,N
         W(I)=B(I)
         B(I)=B(I)*B(I)
      ENDDO
      DO I=1,NEVADD
         EV(I)=ALIMIT
      ENDDO
      ROOTL=-ALIMIT

C
      DO 200 I=1,NEVADD
         ROOTX=ALIMIT
         DO J=I,NEVADD
            ROOTX=MIN(ROOTX,EV(J))
         ENDDO
         EV(I)=ROOTX
 130     TRIAL=(ROOTL+EV(I))*0.5

c         IF(TRIAL.EQ.ROOTL.OR.TRIAL.EQ.EV(I)) GOTO 200
         if (abs(trial-rootl).lt.1E-15 .or.
     $        abs(trial-ev(i)).lt.1E-15) goto 200

         NOMTCH=N
         J=1
 150     F0=A(J)-TRIAL
 160     CONTINUE
         IF(ABS(F0).LT.THETA1) GOTO 170
         IF(F0.GE.0.0) NOMTCH=NOMTCH-1
         J=J+1
         IF(J.GT.N) GOTO 180
         F0=A(J)-TRIAL-B(J-1)/F0
         GOTO160
 170     J=J+2
         NOMTCH=NOMTCH-1
         IF(J.LE.N) GOTO 150
 180     CONTINUE
         IF(NOMTCH.GE.I) GOTO 190
         ROOTL=TRIAL
         GOTO 130
 190     EV(I)=TRIAL
         NOM=MIN0(NEVADD,NOMTCH)
         EV(NOM)=TRIAL
         GOTO 130
 200  CONTINUE

      DO I=1,NEV
         EV(I)=EV(I+NADD)
      ENDDO
C
      DO 600 I=1,NEV
         AROOT=EV(I)
         DO J=1,N
            Y(J)=1.0
         ENDDO
         IF(I.EQ.1) GOTO 250
         IF(ABS(EV(I-1)-AROOT).LT.TOLER) GOTO 260
 250     IA=-1
 260     IA=IA+1
         ELIM1=A(1)-AROOT
         ELIM2=W(1)

         DO J=1,NN
            IF(ABS(ELIM1).LE.ABS(W(J))) GOTO 270
            TA(J)=ELIM1
            TB(J)=ELIM2
            P(J)=0.0
            TEMP=W(J)/ELIM1
            ELIM1=A(J+1)-AROOT-TEMP*ELIM2
            ELIM2=W(J+1)
            GOTO 280
 270        TA(J)=W(J)
            TB(J)=A(J+1)-AROOT
            P(J)=W(J+1)
            TEMP=1.0
            IF(ABS(W(J)).GT.THETA1) TEMP=ELIM1/W(J)
            ELIM1=ELIM2-TEMP*TB(J)
            ELIM2=-TEMP*W(J+1)
 280        B(J)=TEMP
         ENDDO
C
         TA(N)=ELIM1
         TB(N)=0.0
         P(N)=0.0
         P(NN)=0.0
         ITER=1

         IF(IA.NE.0) GOTO 460
C     
 320     L=N+1

         DO 400 J=1,N
            L=L-1
 330        CONTINUE
            IF(N-L-1) 340,350,360
 340        ELIM1=Y(L)
            GOTO 370
 350        ELIM1=Y(L)-Y(L+1)*TB(L)
            GOTO 370
 360        ELIM1=Y(L)-Y(L+1)*TB(L)-Y(L+2)*P(L)
 370        CONTINUE

C
C  OVERFLOW CHECK
            IF(ELIM1.GT.DELBIG .or. ELIM1.lt.-DELBIG) GOTO 380
            TEMP=TA(L)
            IF(ABS(TEMP).LT.DELTA) TEMP=DELTA
            Y(L)=ELIM1/TEMP
            GOTO 400
 380        DO K=1,N
               Y(K)=Y(K)/DELBIG
            ENDDO
            GOTO 330
 400     CONTINUE
C
         GOTO (410,500),ITER
 410     ITER=ITER+1
 420     ELIM1=Y(1)
         DO 450 J=1,NN
            IF(ABS(TA(J)-W(J)).LT.1E-15) GOTO 440
            Y(J)=ELIM1
            ELIM1=Y(J+1)-ELIM1*B(J)
            GOTO 450
 440        Y(J)=Y(J+1)
            ELIM1=ELIM1-Y(J+1)*B(J)
 450     CONTINUE
         Y(N)=ELIM1
         GOTO 320
C     
 460     CONTINUE
         DO J=1,N
            RAND1=MOD(4099.0*RAND1,RPOWER)
            Y(J)=RAND1/RPOW1-1.0
         ENDDO
         GOTO 320


C     
C     ORTHOG TO PREVIOUS

 500     IF(IA.EQ.0) GOTO 550
         DO J1=1,IA
            K=I-J1
            TEMP=0.0
            IPT=(K-1)*N
            DO J=1,N
               IPT=IPT+1
               TEMP=TEMP+Y(J)*VEC(IPT)
            ENDDO
            IPT=(K-1)*N
            DO J=1,N
               IPT=IPT+1
               Y(J)=Y(J)-TEMP*VEC(IPT)
            ENDDO
         ENDDO
 550     GOTO(420,560),ITER
C
C  NORMALIZE

 560     ELIM1=0.0

         DO J=1,N
            ELIM1=MAX(ELIM1,ABS(Y(J)))
         ENDDO
         TEMP=0.0
         DO J=1,N
            ELIM2=Y(J)/ELIM1
            TEMP=TEMP+ELIM2*ELIM2
         ENDDO
         TEMP=DUNITY/(SQRT(TEMP)*ELIM1)
         DO J=1,N
            Y(J)=Y(J)*TEMP
            IF(ABS(Y(J)).LT.DEL1) Y(J)=0.0
         ENDDO
         IPT=(I-1)*N
         DO J=1,N
            IPT=IPT+1
            VEC(IPT)=Y(J)
         ENDDO
 600  CONTINUE
C

      DO I=1,NEV
         IPT=(I-1)*N
         DO J=1,N
            IPT=IPT+1
            Y(J)=VEC(IPT)
         ENDDO
C
         L=N-2
         MK=(N*(N-1))/2-3
         MK1=3
C     
         DO J=1,L
            T=0.0
            K=N-J-1
            M=K+1
            DO KK=M,N
               T=T+DD(MK+KK)*Y(KK)
            ENDDO
            DO KK=M,N
               EPR=T*DD(MK+KK)
               Y(KK)=Y(KK)-EPR-EPR
            ENDDO
            MK=MK-MK1
            MK1=MK1+1
         ENDDO
C
         T=0.0
         DO J=1,N
            T=T+Y(J)*Y(J)
         ENDDO
         XNORM=SQRT(T)
         XNORM1=DUNITY/XNORM
         DO J=1,N
            Y(J)=Y(J)*XNORM1
         ENDDO
C
         IPT=(I-1)*N
         DO J=1,N
            IPT=IPT+1
            VEC(IPT)=Y(J)
         ENDDO
      ENDDO
C
      DO I=1,N
         EV(I)=EV(I)*ANORM
      ENDDO

      RETURN
      END


      SUBROUTINE NORMAL(V,N)
C-----------------------------------------------------------------------
C     NORMALIZES VECTOR V OF LENGTH N
C
      implicit none
      INTEGER N,I
      REAL*8 V(N),C
C
      C=0.0
      DO I=1,N
        C=C+V(I)*V(I)
      ENDDO
      IF(C.LT.1.0D-12) THEN
c         WRITE(*,25) C
         DO I=1,N
            V(I)=0.0
         ENDDO
         RETURN
      ENDIF
 25   FORMAT(' **** WARNING **** TRYING TO NORMALIZE A ZERO VECTOR',
     $       ' NORM=',E12.4/' IT WILL BE ZEROED.')
C
      C=1.0/SQRT(C)
      DO I=1,N
        V(I)=V(I)*C
      ENDDO
C
      RETURN
      END
    
      SUBROUTINE FROTUMEM(R,U)
C-----------------------------------------------------------------------
C     THIS ROUTINE SOLVES THE CONSTRAINED MINIMIZATION EQUATION
C     USING LAGRANGE MULTIPLIERS.
C     BERNARD R. BROOKS
C
      implicit none

C      REAL*8 R(3,3),W(6),A(3,3),B(3,3),U(3,3),SCR(24)
C      INTEGER I,I1,I2,IPT,J,K,JP,JQ,KP,KQ
C      REAL*8 TRACE,EVS,DET
C
C Afra
C I chaged the above to the following to solve a 2x2 matrix instead of 3x3 
      REAL*8 R(2,2),W(4),A(2,2),B(2,2),U(3,3),SCR(24)
      INTEGER I,I1,I2,IPT,J,K,JP,JQ,KP,KQ
      REAL*8 TRACE,EVS,DET
      DET=0.0


C      DET=0.0
C      DO I=1,3
C       I1=I+1
C        IF(I1.GT.3) I1=I1-3
C        I2=I+2
C        IF(I2.GT.3) I2=I2-3
C        DET=DET+R(I,1)*(R(I1,2)*R(I2,3)-R(I2,2)*R(I1,3))
C      ENDDO
C
C Afra
C determinat of a 2x2 matrix
      DET=R(1,1)*R(2,2)-R(1,2)*R(2,1)
      IPT=0
C      DO I=1,3
C         DO J=I,3
C            IPT=IPT+1
C            W(IPT)=0.0
C            DO K=1,3
C               W(IPT)=W(IPT)+R(J,K)*R(I,K)
C            ENDDO
C         ENDDO
C      ENDDO
C
C Afra
C solve for 2x2 matrix
      DO I=1,2
         DO J=I,2
            IPT=IPT+1
            W(IPT)=0.0
            DO K=1,2
               W(IPT)=W(IPT)+R(J,K)*R(I,K)
               WRITE(*,*) W(IPT)
            ENDDO
         ENDDO
      ENDDO

C      TRACE=W(1)+W(4)+W(6)
C Afra
C trace of a 2x2 matrix
      TRACE=W(1)+W(3)

      IF(TRACE.LT.3.D-6) THEN
        DO I=1,3
          DO J=1,3
            U(I,J)=0.D0
          ENDDO
          U(I,I)=1.D0
        ENDDO
        RETURN
      ENDIF
C
C      CALL DIAGQ(3,3,W,A,SCR(4),SCR(7),SCR(10),SCR(13),SCR(1),
C     $           SCR(16),SCR(19),SCR(22),0)
C
C Afra
C diagonalizes a 2x2 matrix instead of 3x3 matrix
      CALL DIAGQ(2,2,W,A,SCR(4),SCR(7),SCR(10),SCR(13),SCR(1),
     $           SCR(16),SCR(19),SCR(22),0)

C      DO I=1,3
C        SCR(I)=SQRT(ABS(SCR(I)))
C        IF(SCR(I).LT.1.D-6) SCR(I)=1.D-6
C      ENDDO
C Afra
C changed to 2x2
      DO I=1,2
        SCR(I)=SQRT(ABS(SCR(I)))
        IF(SCR(I).LT.1.D-6) SCR(I)=1.D-6
      ENDDO
C
      IF(DET.LT.0.D0) SCR(1)=-SCR(1)
C
C      DO J=1,3
C        EVS=SCR(J)
C        DO I=1,3
C           B(I,J)=0.0
C           DO K=1,3
C              B(I,J)=B(I,J)+R(K,I)*A(K,J)/EVS
C           ENDDO
C         ENDDO
C      ENDDO
C
C Afra
C changed to 2x2
      DO J=1,2
        EVS=SCR(J)
        DO I=1,2
           B(I,J)=0.0
           DO K=1,2
              B(I,J)=B(I,J)+R(K,I)*A(K,J)/EVS
           ENDDO
         ENDDO
      ENDDO

      DET=0.0
C      DO I=1,3
C        I1=I+1
C        IF(I1.GT.3) I1=I1-3
C        I2=I+2
C        IF(I2.GT.3) I2=I2-3
C        DET=DET+A(I,1)*(A(I1,2)*A(I2,3)-A(I2,2)*A(I1,3))
C      ENDDO
C Afra
C changed to 2x2
      DET=A(1,1)*A(2,2)-A(2,1)*A(1,2)

C      DO J=1,3
C        IF(ABS(SCR(J)).LE.1.D-6) THEN
C          JP=J+1
C          JQ=J+2
C          IF(JP.GT.3) JP=JP-3
C          IF(JQ.GT.3) JQ=JQ-3
C          DO K=1,3
C            KP=K+1
C            KQ=K+2
C            IF(KP.GT.3) KP=KP-3
C            IF(KQ.GT.3) KQ=KQ-3
C            B(K,J)=B(KP,JP)*B(KQ,JQ)-B(KP,JQ)*B(KQ,JP)
C            IF(DET.LT.0.D0) B(K,J)=-B(K,J)
C          ENDDO
C        ENDIF
C        CALL NORMAL(B(1,J),3)
C      ENDDO

C Afra
C rewrite to 2x2
      DO J=1,2
        IF(ABS(SCR(J)).LE.1.D-6) THEN
          JP=J+1
          JQ=J+2
          IF(JP.GT.2) JP=JP-2
          IF(JQ.GT.2) JQ=JQ-2
          DO K=1,2
            KP=K+1
            KQ=K+2
            IF(KP.GT.2) KP=KP-2
            IF(KQ.GT.2) KQ=KQ-2
            B(K,J)=B(KP,JP)*B(KQ,JQ)-B(KP,JQ)*B(KQ,JP)
            IF(DET.LT.0.D0) B(K,J)=-B(K,J)
          ENDDO
        ENDIF
        CALL NORMAL(B(1,J),2)
      ENDDO

C
      DO J=1,3
        DO I=1,3
          U(I,J)=0.D0
        ENDDO
      ENDDO
      DO J=1,3
        DO I=1,3
C          DO 30 K=1,3
C Afra
C changed to 2x2
          DO K=1,2
            U(I,J)=U(I,J)+A(I,K)*B(J,K)
          ENDDO
        ENDDO
      ENDDO
C 30   CONTINUE
C
C      DO J=1,3
C Afra
      DO J=1,2
C        CALL NORMAL(U(1,J),3)
C Afra
        CALL NORMAL(U(1,J),2)
      ENDDO
C Afra
            U(1,3)=0.D0
            U(2,3)=0.D0
            U(3,3)=1.D0
            U(3,1)=0.D0
            U(3,2)=0.D0
C
C     CHECK TO INSURE UNITY (AS OPPOSED TO ANTI-UNITARY)
      DET=0.0
      DO I=1,3
        I1=I+1
        IF(I1.GT.3) I1=I1-3
        I2=I+2
        IF(I2.GT.3) I2=I2-3
        DET=DET+U(I,1)*(U(I1,2)*U(I2,3)-U(I2,2)*U(I1,3))
      ENDDO
      
c      IF(ABS(DET-1.0).GT.1.D-4) WRITE(*,55) DET

 55   FORMAT(/' ***** WARNING ***** FROM FROTUMEM. ROTATION MATRIX IS ',
     $       'NOT UNITARY.'/,'  THE DETERMINANT IS',F14.8/)
C
      RETURN
      END
