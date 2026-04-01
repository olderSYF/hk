*USER SUBROUTINES
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

       ! Arguments:
       !          I/O  Type
       !  PROPS    I   R()  : List with model parameters
       !  DSTRAN   I   R()  : Strain increment
       !  DDSDDE   O   R(,) : Material stiffness matrix
       !  STRESS  I/O  R()  : stresses
       !  STATEV  I/O  R()  : state variables
       !
       ! PROPS layout (1:8):
       !  1 : G       Shear modulus (MPa)
       !  2 : ENU     Poisson's ratio
       !  3 : SIGCI   Intact rock UCS (MPa, positive)
       !  4 : XMI     Intact rock material constant mi
       !  5 : GSI     Geological Strength Index (10-100)
       !  6 : DISTD   Disturbance factor D (0-1)
       !  7 : SPSI    sin(dilation angle) (0-1)
       !  8 : TENS    Tensile strength cutoff (MPa, positive; 0=auto)
       !
       ! STATEV layout (1:4):
       !  1 : Accumulated equivalent plastic strain
       !  2 : Yield flag (0=elastic, 1=shear, 2=tension)
       !  3 : Current equivalent friction angle (degrees)
       !  4 : Current equivalent cohesion (MPa)

        !---  Local variables
      DIMENSION dSig(NTENS), Sig(NTENS)
      DIMENSION SigC(NTENS)
      DIMENSION SigE(NTENS),SigEQ(NTENS)
      Parameter (Ip=-1)
      DIMENSION xN1(3),xN2(3),xN3(3),Tmp1(3),Tmp2(3),Tmp3(3)

      IAPEX=0
      ITENS=0
      TauMax=1
      IPL=-1

        ! Read material parameters
        IntGlo= NPT
        G     = PROPS(1)
        ENU   = PROPS(2)
        VNU   = PROPS(2)
        SIGCI = PROPS(3)
        XMI   = PROPS(4)
        GSI   = PROPS(5)
        DISTD = PROPS(6)
        SPSI  = PROPS(7)
        TENS  = PROPS(8)

        one = 1.0d0
        two = 2.0d0

!       Compute Hoek-Brown rock-mass constants
!         mb = mi*exp((GSI-100)/(28-14D))
!         s  = exp((GSI-100)/(9-3D))
!         a  = 0.5+(1/6)*(exp(-GSI/15)-exp(-20/3))
        XMB = XMI * exp( (GSI-100.0d0) / (28.0d0-14.0d0*DISTD) )
        XS  = exp( (GSI-100.0d0) / (9.0d0-3.0d0*DISTD) )
        XA  = 0.5d0 + (one/6.0d0)*( exp(-GSI/15.0d0)
     &        - exp(-20.0d0/3.0d0) )

        ! Auto tensile strength: sigma_t = s*sigci/mb
        if (TENS <= 0.0d0) then
          TENS = XS * SIGCI / XMB
        end if

        ! calculate elastic stress increment
        ! (DSig = D_elastic * strain increment DEps)
        FAC = two * G / ( one - two * ENU )
        D1 = FAC * ( one - ENU )
        D2 = FAC * ENU
        DSTRANVOL = DSTRAN(1) + DSTRAN(2) + DSTRAN(3)
        dSig(1) = (D1 - D2) * DSTRAN(1) + D2 * DSTRANVOL
        dSig(2) = (D1 - D2) * DSTRAN(2) + D2 * DSTRANVOL
        dSig(3) = (D1 - D2) * DSTRAN(3) + D2 * DSTRANVOL
        dSig(4) = G * DSTRAN(4)
        if (NTENS == 6) then
        dSig(5) = G * DSTRAN(5)
        dSig(6) = G * DSTRAN(6)
        end if

        ! elastic stress
        SigE = STRESS + dSig
        Stress=SigE
        SigEQ = SigE
        SigC = SigE

        DDSDDE = 0.0
        DDSDDE(1:3,1:3) = D2
        DDSDDE(1,1) = D1
        DDSDDE(2,2) = D1
        DDSDDE(3,3) = D1
        DDSDDE(4,4) = G
        if (NTENS == 6) then
          DDSDDE(5,5) = G
          DDSDDE(6,6) = G
        end if

!------------ Calculate principle stresses and their direction
!             Note: PrnSig called WITH NTENS (standalone DLL version)

      call PrnSig(1,NTENS,SigE,xN1,xN2,xN3,Sig1,Sig2,Sig3,P,Q)

!       Convert to compression-positive for HB yield evaluation
!         S1C = -Sig3  (most compressive principal, positive)
!         S3C = -Sig1  (least compressive principal, positive)
        S1C = -Sig3
        S3C = -Sig1

!       HB yield function (compression-positive):
!         F = S1C - S3C - sigci*(mb*S3C/sigci + s)^a
!       Tension cutoff (tension-positive): FT3 = Sig3 - TENS
        TENS_ARG = XMB * S3C / SIGCI + XS
        if (TENS_ARG < 0.0d0) TENS_ARG = 0.0d0
        F_HB = S1C - S3C - SIGCI * TENS_ARG**XA

        FT1 = SIG1 - TENS
        FT2 = SIG2 - TENS
        FT3 = SIG3 - TENS
      IF (F_HB > 1d-12 .OR. FT3 > 1d-12) GOTO 240
      IF (IPL  == 0) GOTO 360
!
!         ***POINTS CHANGING FROM PLASTIC TO ELASTIC STATE***
!
      IPL =0
      GOTO 360
!
!                   *** Plastic stress points ***
!
  240 IPL =1

!       Instantaneous MC equivalent parameters from HB tangent:
!         dS1/dS3 = 1 + a*mb*(mb*S3/sigci+s)^(a-1)
!         sin(phi_eq) = (dS1/dS3-1)/(dS1/dS3+1)
        TENS_ARG = XMB * S3C / SIGCI + XS
        if (TENS_ARG < 0.0d0) TENS_ARG = 0.0d0

        if (TENS_ARG > 0.0d0) then
          DERIV = one + XA * XMB * TENS_ARG**(XA - one)
        else
          DERIV = one
        end if

        if (DERIV < one + 1.0d-10) DERIV = one + 1.0d-10

        SPHI = (DERIV - one) / (DERIV + one)
        if (SPHI > 0.999d0) SPHI=0.999d0
        CPHI = sqrt(one - SPHI*SPHI)

        S1C_YIELD = S3C + SIGCI * TENS_ARG**XA
        COHS = 0.5d0 * (S1C_YIELD - DERIV * S3C) * (one - SPHI)
        if (COHS < 0.0d0) COHS = 0.0d0

        ! Re-evaluate MC yield functions with equivalent parameters
      F21 = 0.5d0*(SIG2-SIG1) + 0.5d0*(SIG2+SIG1)*SPHI - COHS
      F32 = 0.5d0*(SIG3-SIG2) + 0.5d0*(SIG3+SIG2)*SPHI - COHS
      F31 = 0.5d0*(SIG3-SIG1) + 0.5d0*(SIG3+SIG1)*SPHI - COHS

      call PrnSig( 0,NTENS,SigEQ,Tmp1,Tmp2,Tmp3,
     &         SigV1,SigV2,SigV3,Dum1,Dum2)

      FF1 = 0.5d0*    (SIGV2-SIGV1) + 0.5d0*(SIGV2+SIGV1)*SPHI - COHS
      FF2 = 0.5d0*DABS(SIGV2-SIGV3) + 0.5d0*(SIGV2+SIGV3)*SPHI - COHS
      FF3 = 0.5d0*DABS(SIGV3-SIGV1) + 0.5d0*(SIGV3+SIGV1)*SPHI - COHS

      IF (FF1 > 1d-12 .OR. FF2 > 1d-12 .OR. FF3 > 1d-12) GOTO 250

!       ****** PLASTIC NEGATIVE POINTS ******
      IPL =-1
!
!         *** STRESSES ARE BROUGHT BACK TO THE YIELD SURFACE ***
!
  250 CONTINUE
      VNUI0 = VNU
      VNUI1 = 1 - 2*VNUI0
      VNUI2 = 1 + 2*VNUI0
      VNUIQ = VNUI2 / VNUI1
!
      DUM = SPSI / VNUI1
      PSIMIN = G *(-1+DUM)
      PSIMET = G *( 1+DUM)
      PSINU = 2*G *VNUI0*DUM
!
      HA= ( 1 - SPHI + SPSI - SPHI * SPSI) * SIG1 +
     &      (-2 - 2/VNUI1*SPHI*SPSI) * SIG2 +
     &     ( 1 - SPHI - SPSI + VNUIQ * SPHI * SPSI) * SIG3 +
     &     2*(1 + SPSI) * COHS

      HB= ( 1 + SPHI + SPSI + VNUIQ * SPHI * SPSI) * SIG1 -
     &     ( 2 + 2/VNUI1*SPHI*SPSI) * SIG2 +
     &     ( 1 + SPHI - SPSI - SPHI * SPSI) * SIG3 -
     &     2*(1 - SPSI) * COHS

      HAB= (VNUI1+SPSI)*SIG1+
     &      (VNUI1-SPSI)*(SIG3-TENS)-
     &      (VNUI1+SPSI)*(TENS*(1+SPHI)-2*COHS)/(1-SPHI)

      HBA=(1-VNUI0)*SIG1-VNUI0*SIG3-(1-VNUI0)*(TENS*(1+SPHI)-
     &      2*COHS)/(1-SPHI)+VNUI0*TENS

      HAO=(1-VNUI0)*SIG2-VNUI0*SIG3-VNUI1*TENS

      HOC=(1-VNUI0)*SIG1-VNUI0*SIG3-VNUI1*TENS

      HAA=(VNUI1+VNUI2*SPSI)*(SIG1-(TENS*(1+SPHI)-2*COHS)/
     &     (1-SPHI))+
     &     (VNUI1-SPSI)*(SIG2-TENS)+(VNUI1-SPSI)*(SIG3-TENS)

      HAAB=VNUI0*(SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))-
     &      (SIG2-TENS)+VNUI0*(SIG3-TENS)

      HAAO=(SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))-
     &      VNUI0*(SIG2-TENS)-VNUI0*(SIG3-TENS)

      HBB=(VNUI1+SPSI)*(SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))+
     &     (VNUI1+SPSI)*(SIG2-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))+
     &     (VNUI1-VNUI2*SPSI)*(SIG3-TENS)

!        ***** DETERMINING RETURN AREA ****

      IAREA=0
      IASIGN=0

      IF (F31 > 0 .AND. HA >= 0 .AND. HB < 0 .AND. HAB < 0) THEN
        IAREA=2
        IASIGN=IASIGN+1
      ENDIF

      IF (F31 > 0 .AND. HB >= 0) THEN
        IASIGN=IASIGN+1
        IF (HBB < 0) THEN
          IAREA=1
          GOTO 260
        ELSE
          IAREA=8
        ENDIF
      ENDIF

      IF (F31 > 0 .AND. HA < 0) THEN
        IASIGN=IASIGN+1
        IF (HAA < 0) THEN
          IAREA=3
        ELSE
          IAREA=7
        ENDIF
      ENDIF

      IF (FT3 > 1d-12 .AND.
     &     HBA >= 0 .AND.
     &     HAO < 0 .AND.
     &     HOC < 0) THEN
        IAREA=4
        IASIGN=IASIGN+1
      ENDIF

      IF (FT3 > 1d-12 .AND. HAB >= 0 .AND. HBA < 0) THEN
        IASIGN=IASIGN+1
        IF (HAAB > 0) THEN
          IAREA=5
        ELSE
          IF (IAREA == 0) IAREA=7
        ENDIF
      ENDIF

      IF (FT3 > 1d-12 .AND. HAO >= 0) THEN
        IASIGN=IASIGN+1
        IF (HAAO >= 0) THEN
          IAREA=6
        ELSE
          IF (IAREA == 0) IAREA=7
        ENDIF
      ENDIF

 260  CONTINUE

      DSP1=0
      DSP2=0
      DSP3=0

!         *** EXTENSION POINT ***

      IF (IAREA == 1) THEN
        A11 =       G *(1+SPHI*SPSI/VNUI1)
        A12 = 0.5d0*G *(1+SPHI+SPSI+VNUIQ*SPHI*SPSI)
        DETER  = A11*A11-A12*A12
        RLAM31 = (F31*A11-F32*A12) / DETER
        RLAM32 = (F32*A11-F31*A12) / DETER
        RLAM21 = 0
        DSP1 = RLAM31*PSIMIN + RLAM32*PSINU
        DSP2 = RLAM31*PSINU  + RLAM32*PSIMIN
        DSP3 = RLAM31*PSIMET + RLAM32*PSIMET
      ENDIF

!         *** REGULAR YIELD SURFACE ***

      IF (IAREA == 2) THEN
        A11 = G *(1+SPHI*SPSI/VNUI1)
        RLAM31 = F31 / A11
        RLAM21 = 0
        RLAM32 = 0
        DSP1 = RLAM31*PSIMIN
        DSP2 = RLAM31*PSINU
        DSP3 = RLAM31*PSIMET
      END IF

!         *** COMPRESSION POINT ***

      IF (IAREA == 3) THEN
        A11 =       G *(1+SPHI*SPSI/VNUI1)
        A12 = 0.5d0*G *(1-SPHI-SPSI+VNUIQ*SPHI*SPSI)
        DETER  = A11*A11-A12*A12
        RLAM31 = (F31*A11-F21*A12) / DETER
        RLAM21 = (F21*A11-F31*A12) / DETER
        RLAM32 = 0
        DSP1 = RLAM31*PSIMIN + RLAM21*PSIMIN
        DSP2 = RLAM31*PSINU  + RLAM21*PSIMET
        DSP3 = RLAM31*PSIMET + RLAM21*PSINU
      END IF

      IF (IAREA == 4) THEN
        DUM=   2*G /VNUI1
        RLAMT3=FT3/(DUM*(1-VNUI0))
        DSP1 = DUM*RLAMT3*VNUI0
        DSP2 = DSP1
        DSP3 = DUM*RLAMT3*(1-VNUI0)
        ITENS=1
        IPL =2
      ENDIF

      IF (IAREA == 5) THEN
        A11 = VNUI1+SPHI*SPSI
        A12 = VNUI1+SPHI
        A21 = VNUI1+SPSI
        A22 = 2*(1-VNUI0)
        DUM = G /VNUI1
        DETER  = A11*A22-A12*A21
        RLAM31 = (F31*A22-FT3*A12) / DETER / DUM
        RLAMT3 = (FT3*A11-F31*A21) / DETER / DUM
        DSP1 = DUM*(RLAM31*(-VNUI1+SPSI)+RLAMT3*2*VNUI0)
        DSP2 = DUM*(RLAM31*2*VNUI0*SPSI+RLAMT3*2*VNUI0)
        DSP3 = DUM*(RLAM31*(VNUI1+SPSI)+RLAMT3*2*(1-VNUI0))
        ITENS=1
        IPL =2
      ENDIF

      IF (IAREA == 6) THEN
        RLAMT2 = (FT2*(1-VNUI0)-FT3*VNUI0)/(2*G)
        RLAMT3 = (FT3*(1-VNUI0)-FT2*VNUI0)/(2*G)
        DUM = 2*G/VNUI1
        DSP1 = DUM*VNUI0*(RLAMT2+RLAMT3)
        DSP2 = DUM*(RLAMT2*(1-VNUI0)+RLAMT3*   VNUI0 )
        DSP3 = DUM*(RLAMT2*   VNUI0 +RLAMT3*(1-VNUI0))
        ITENS=1
        IPL =2
      ENDIF

      IF (IAREA == 7) THEN
        DSP1 = SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI)
        DSP2 = SIG2-TENS
        DSP3 = SIG3-TENS
        ITENS=1
        IPL =2
      ENDIF

      IF (IAREA == 8) THEN
        DSP1 = SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI)
        DSP2 = SIG2-(SIG1-DSP1)
        DSP3 = SIG3-TENS
        ITENS=1
        IPL=2
      ENDIF

      IF (IAREA == 9) THEN
        DSP1 = SIG1-TENS
        DSP2 = SIG2-TENS
        DSP3 = SIG3-TENS
        IAPEX=1
        IPL=2
      ENDIF

!--------Check if principle stresses are on the yield surface
      SIG1R=SIG1-DSP1
      SIG2R=SIG2-DSP2
      SIG3R=SIG3-DSP3

      F21R = 0.5d0*(SIG2R-SIG1R) + 0.5d0*(SIG2R+SIG1R)*SPHI - COHS
      F32R = 0.5d0*(SIG3R-SIG2R) + 0.5d0*(SIG3R+SIG2R)*SPHI - COHS
      F31R = 0.5d0*(SIG3R-SIG1R) + 0.5d0*(SIG3R+SIG1R)*SPHI - COHS
      FT1R = SIG1R - TENS
      FT2R = SIG2R - TENS
      FT3R = SIG3R - TENS

      IF ( F21R > 1d-6.OR.
     &      F32R > 1d-6.OR.
     &      F31R > 1d-6.OR.
     &      FT1R > 1d-6.OR.
     &      FT2R > 1d-6.OR.
     &      FT3R > 1d-6    ) THEN
        ICREC=0

        IF (IAREA == 5.AND.F32R > 1E-6) THEN
!------------ZONE 8
          ICREC=ICREC+1
          DSP1 = SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI)
          DSP2 = SIG2-(SIG1-DSP1)
          DSP3 = SIG3-TENS
        ENDIF

        IF (IAREA == 6.AND.FT1R > 1d-6) THEN
!------------ZONE 9
          ICREC=ICREC+1
          DSP1 = SIG1-TENS
          DSP2 = SIG2-TENS
          DSP3 = SIG3-TENS
          ITENS=0
          IAPEX=1
        ENDIF
      ENDIF

!--------Check again if the calculated stresses are on the yield surface
      SIG1R=SIG1-DSP1
      SIG2R=SIG2-DSP2
      SIG3R=SIG3-DSP3

      F21 = 0.5d0*(SIG2R-SIG1R) + 0.5d0*(SIG2R+SIG1R)*SPHI - COHS
      F32 = 0.5d0*(SIG3R-SIG2R) + 0.5d0*(SIG3R+SIG2R)*SPHI - COHS
      F31 = 0.5d0*(SIG3R-SIG1R) + 0.5d0*(SIG3R+SIG1R)*SPHI - COHS
      FT1 = SIG1R - TENS
      FT2 = SIG2R - TENS
      FT3 = SIG3R - TENS

!         *** Computing Cartesian stress components ***

      Call CarSig(Sig1R,Sig2R,Sig3R,xN1,xN2,xN3,NTENS,SigC)
      STRESS=SigC

!------- Calculate TAUMAX for checking accuracy on plastic point
      TAUMAX=0.5d0*(SIG3-DSP3-SIG1+DSP1)
      TAUMAX = MAX(TAUMAX,COHS,0.5d0)

!------- Update state variables
!        STATEV(1): accumulated equivalent plastic strain
!        STATEV(2): yield flag (0=elastic, 1=shear, 2=tension)
!        STATEV(3): equivalent friction angle (degrees)
!        STATEV(4): equivalent cohesion (MPa)

      DEPS_EQ = SQRT( (2.0d0/3.0d0)*(DSP1**2+DSP2**2+DSP3**2) )
      STATEV(1) = STATEV(1) + DEPS_EQ

      if (IPL == 1) then
        STATEV(2) = 1.0d0
      else if (IPL == 2) then
        STATEV(2) = 2.0d0
      else
        STATEV(2) = 0.0d0
      end if

      if (SPHI > 0.0d0 .AND. SPHI < 1.0d0) then
        STATEV(3) = ASIN(SPHI)*180.0d0/3.14159265358979d0
      else
        STATEV(3) = 0.0d0
      end if
      STATEV(4) = COHS

  360 CONTINUE

      RETURN

      End
!----------------------------------------------------------------------
      Subroutine MatTranspose (A,Ia,AT,IAt,N1,N2)
      implicit none
      integer :: Ia, N1, N2, IAt
      double precision :: A(Ia,*),AT(Iat,*)
      integer :: I, J

      Do I=1,N1
        Do J=1,N2
          AT(I,J)=A(J,I)
        End Do
      End Do

      End
!***********************************************************************
      Subroutine MatMat(A,Ia,B,Ib,N1,N2,N3,C,Ic)
!  Matrix multiplication C(N1,N3) = A(N1,N2) * B(N2,N3)
      implicit none
      integer :: Ia, Ib, Ic, N1, N2, N3
      double precision :: A(Ia,*),B(Ib,*),C(Ic,*)
      integer :: I, J, K
      double precision :: Sum

      Do I=1,N1
        Do J=1,N3
          Sum=0.0d0
          Do K=1,N2
            Sum=Sum+A(I,K)*B(K,J)
          End Do
          C(I,J)=Sum
        End Do
      End Do

      End
!***********************************************************************
      subroutine CarSig(S1, S2, S3, xN1, xN2, xN3, ntens, Stress)
!  Returns Cartesian stresses from principal stresses and directions.
!
!  S1, S2, S3      I   R     principal stress
!  xN1, xN2, xN3   I   R()   principal direction
!  Stress          O   R()   cartesian stress

      implicit none

        double precision, intent(in) :: S1, S2, S3
        double precision, intent(in), dimension(3) :: xN1, xN2, xN3
        integer, intent(in):: ntens
        double precision, intent(out), dimension(ntens) :: Stress

        integer :: I
        integer :: IDim
        double precision, dimension(:, :), allocatable :: SM, T, TT, STT

        IDim = 3
        allocate(SM(IDim,IDim), T(IDim,IDim),
     &           TT(IDim,IDim), STT(IDim,IDim))

        do I = 1,3
            T(I,1) = xN1(I)
            T(I,2) = xN2(I)
            T(I,3) = xN3(I)
            TT(1,I) = T(I,1)
            TT(2,I) = T(I,2)
            TT(3,I) = T(I,3)
        end do

        SM = 0.0
        SM(1,1) = S1
        SM(2,2) = S2
        SM(3,3) = S3

        call MatMat(SM, IDim, TT,  IDim, IDim, IDim, IDim , STT, IDim)
        call MatMat(T,  IDim, STT, IDim, IDim, IDim, IDim , SM,  IDim)

        do I = 1, IDim
          Stress(I) = SM(I, I)
        end do

        Stress(4) = SM(2, 1)
        if (ntens == 6) then
          Stress(5) = SM(3, 2)
          Stress(6) = SM(3, 1)
        end if

      end subroutine CarSig
!***********************************************************************
      Subroutine PrnSig(IOpt,NTENS,S,xN1,xN2,xN3,S1,S2,S3,P,Q)
!  Calculates principal stresses and their directions.
!  IOpt=1: also computes directions; IOpt=0: values only
!
!  S       I   R(NTENS)  stress vector (tension-positive)
!  xN1..3  O   R(3)      principal directions (unit vectors)
!  S1<=S2<=S3            principal stresses in ascending order
!  P                     mean effective stress (positive in compression)
!  Q                     deviatoric stress invariant
      implicit none
      integer, intent(in) :: IOpt, NTENS
      double precision, intent(in)  :: S(NTENS)
      double precision, intent(out) :: xN1(3),xN2(3),xN3(3)
      double precision, intent(out) :: S1,S2,S3,P,Q

      double precision :: A(3,3), V(3,3), D(3)
      integer :: I, J
      double precision :: OffDiag, Tol
      double precision :: CC, SS, TT, Tau, Theta, Apq
      double precision :: Tmp1, Tmp2
      integer :: p_idx, q_idx, r_idx, Iter, MaxIter

!     Build symmetric 3x3 stress tensor from Voigt vector
!     Voigt order: s11,s22,s33,s12(,s23,s13) — tension-positive
      A(1,1) = S(1)
      A(2,2) = S(2)
      A(3,3) = S(3)
      A(1,2) = S(4)
      A(2,1) = S(4)
      if (NTENS == 6) then
        A(2,3) = S(5)
        A(3,2) = S(5)
        A(1,3) = S(6)
        A(3,1) = S(6)
      else
        A(2,3) = 0.0d0
        A(3,2) = 0.0d0
        A(1,3) = 0.0d0
        A(3,1) = 0.0d0
      end if

!     Initialise eigenvector matrix V to identity
      V = 0.0d0
      V(1,1) = 1.0d0
      V(2,2) = 1.0d0
      V(3,3) = 1.0d0

!     Classical Jacobi eigenvalue iteration for 3x3 symmetric matrix
      MaxIter = 100
      Tol = 1.0d-12
      Do Iter = 1, MaxIter
        OffDiag = ABS(A(1,2)) + ABS(A(1,3)) + ABS(A(2,3))
        IF (OffDiag < Tol) EXIT
        Do p_idx = 1, 2
          Do q_idx = p_idx+1, 3
            IF (ABS(A(p_idx,q_idx)) < 1.0d-15) CYCLE
            Theta = 0.5d0*(A(q_idx,q_idx)-A(p_idx,p_idx))
     &              / A(p_idx,q_idx)
            IF (Theta >= 0.0d0) THEN
              TT = 1.0d0/(Theta + SQRT(1.0d0+Theta*Theta))
            ELSE
              TT = 1.0d0/(Theta - SQRT(1.0d0+Theta*Theta))
            END IF
            CC  = 1.0d0/SQRT(1.0d0+TT*TT)
            SS  = TT*CC
            Tau = SS/(1.0d0+CC)
!           Save off-diagonal element before zeroing
            Apq = A(p_idx,q_idx)
!           Update diagonal elements
            A(p_idx,p_idx) = A(p_idx,p_idx) - TT*Apq
            A(q_idx,q_idx) = A(q_idx,q_idx) + TT*Apq
!           Zero the (p,q) and (q,p) elements
            A(p_idx,q_idx) = 0.0d0
            A(q_idx,p_idx) = 0.0d0
!           Update remaining off-diagonal elements and eigenvectors
            Do r_idx = 1, 3
              IF (r_idx /= p_idx .AND. r_idx /= q_idx) THEN
!               Save old values before updating
                Tmp1 = A(r_idx,p_idx)
                Tmp2 = A(r_idx,q_idx)
                A(r_idx,p_idx) = Tmp1 - SS*(Tmp2+Tau*Tmp1)
                A(p_idx,r_idx) = A(r_idx,p_idx)
                A(r_idx,q_idx) = Tmp2 + SS*(Tmp1-Tau*Tmp2)
                A(q_idx,r_idx) = A(r_idx,q_idx)
              END IF
!             Update eigenvectors (use saved old values)
              Tmp1 = V(r_idx,p_idx)
              Tmp2 = V(r_idx,q_idx)
              V(r_idx,p_idx) = Tmp1 - SS*(Tmp2+Tau*Tmp1)
              V(r_idx,q_idx) = Tmp2 + SS*(Tmp1-Tau*Tmp2)
            End Do
          End Do
        End Do
      End Do

      D(1) = A(1,1)
      D(2) = A(2,2)
      D(3) = A(3,3)

!     Sort eigenvalues D(1)<=D(2)<=D(3) and reorder eigenvectors
      Do I = 1, 2
        Do J = I+1, 3
          If (D(J) < D(I)) Then
            TT  = D(I)
            D(I) = D(J)
            D(J) = TT
            Do p_idx = 1, 3
              TT = V(p_idx,I)
              V(p_idx,I) = V(p_idx,J)
              V(p_idx,J) = TT
            End Do
          End If
        End Do
      End Do

      S1 = D(1)
      S2 = D(2)
      S3 = D(3)

      xN1(1) = V(1,1); xN1(2) = V(2,1); xN1(3) = V(3,1)
      xN2(1) = V(1,2); xN2(2) = V(2,2); xN2(3) = V(3,2)
      xN3(1) = V(1,3); xN3(2) = V(2,3); xN3(3) = V(3,3)

      P = -(S1+S2+S3)/3.0d0
      Q = SQRT(0.5d0*((S1-S2)**2+(S2-S3)**2+(S3-S1)**2))

      End Subroutine PrnSig
