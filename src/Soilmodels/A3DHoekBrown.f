    module ModHoekBrown
    contains

    Subroutine ESM_HB(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER,&
     DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)
!**********************************************************************
!
!  Function: ESM wrapper for the Hoek-Brown constitutive model.
!            Follows the same pattern as ESM_MC in A3DMohrCoulombStandard.f.
!
!  PROPS layout (1:8):
!    1 : G       Shear modulus (MPa)
!    2 : ENU     Poisson's ratio
!    3 : SIGCI   Intact rock UCS (MPa, positive)
!    4 : XMI     Intact rock material constant mi
!    5 : GSI     Geological Strength Index (10-100)
!    6 : DISTD   Disturbance factor D (0-1)
!    7 : SPSI    sin(dilation angle) (0-1)
!    8 : TENS    Tensile strength cutoff (MPa, positive; 0 = auto)
!
!  STATEV layout (1:4):
!    1 : Accumulated equivalent plastic strain
!    2 : Yield flag (0=elastic, 1=shear, 2=tension)
!    3 : Current equivalent friction angle (degrees)
!    4 : Current equivalent cohesion (MPa)
!
!**********************************************************************

      implicit double precision (a-h, o-z)
      integer :: NTENS, NSTATEV, NADDVAR, NPROPS, NPT, NOEL, IDSET, NUMBEROFPHASES
      double precision :: EUNLOADING, PLASTICMULTIPLIER
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS), DSTRAN(NTENS), STATEV(NSTATEV), &
                ADDITIONALVAR(NADDVAR), PROPS(NPROPS)

!---  Local variables required in standard UMAT interface
      integer :: IStep, TimeStep
      double precision, dimension(:), allocatable :: ddsddt
      double precision, dimension(:), allocatable :: drplde
      double precision, dimension(:), allocatable :: stran
      double precision, dimension(:), allocatable :: time
      double precision, dimension(:), allocatable :: predef
      double precision, dimension(:), allocatable :: dpred
      double precision, dimension(:), allocatable :: coords
      double precision, dimension(:,:), allocatable :: ddsdde
      double precision, dimension(:,:), allocatable :: drot
      double precision, dimension(:,:), allocatable :: dfgrd0
      double precision, dimension(:,:), allocatable :: dfgrd1
      double precision :: sse, spd, scd
      double precision :: rpl
      double precision :: drpldt
      double precision :: pnewdt, dtime, temp, dtemp, celent
      double precision :: Porosity, WaterPressure, WaterPressure0
      double precision :: GasPressure, GasPressure0, DegreeSaturation
      integer :: ndi, nshr, layer, kspt, kstep, kinc

      allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), &
                predef(1), dpred(1), coords(3), ddsdde(ntens,ntens), &
                drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )

!     Initialization
      Eunloading = 0.0d0
      PlasticMultiplier = 0.0d0

!     Rename additional variables
      Porosity         = AdditionalVar(1)
      WaterPressure    = AdditionalVar(2)
      WaterPressure0   = AdditionalVar(3)
      GasPressure      = AdditionalVar(4)
      GasPressure0     = AdditionalVar(5)
      DegreeSaturation = AdditionalVar(6)
      time(1)          = AdditionalVar(7)   ! TotalRealTime
      time(2)          = AdditionalVar(8)   ! OverallTotalTime
      dtime            = AdditionalVar(9)   ! TimeIncrement
      IStep            = AdditionalVar(10)
      TimeStep         = AdditionalVar(11)  ! Very first step: IStep=1, TimeStep=1

!     Call the UMAT_HB subroutine
      call UMAT_HB(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, &
                   stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname, &
                   ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, &
                   celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)

!     Eunloading: used to define the maximum time step
      Eunloading = max(ddsdde(1,1), ddsdde(2,2), ddsdde(3,3))

!     PlasticMultiplier: yield flag stored in STATEV(2) for result output
      PlasticMultiplier = statev(2)

      return

    end subroutine ESM_HB

!***********************************************************************
    Subroutine UMAT_HB(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
      RPL,DDSDDT,DRPLDE,DRPLDT, &
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!**********************************************************************
!
!  Function: Hoek-Brown constitutive model (built-in version).
!            Algorithm:
!              1. Compute HB parameters: mb, s, a
!              2. Auto-compute tensile strength if TENS=0
!              3. Elastic trial stress
!              4. Decompose into principal stresses via PrnSig
!              5. Convert to compression-positive convention for HB check
!              6. Evaluate HB yield function and tension cutoff
!              7. If plastic: compute instantaneous MC equivalent params
!              8. Apply MC return-mapping (Areas 1-9) with equiv params
!              9. Reconstruct Cartesian stress via CarSig
!             10. Update state variables
!
!**********************************************************************

      implicit double precision (a-h, o-z)

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS), STATEV(NSTATEV), &
                DDSDDE(NTENS,NTENS), DDSDDT(NTENS), DRPLDE(NTENS), &
                STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1), &
                PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)

      ! Local arrays
      DIMENSION dSig(6), SigC(6)
      DIMENSION SigE(6), SigEQ(6)
      Parameter (Ip=-1)
      DIMENSION xN1(3), xN2(3), xN3(3), Tmp1(3), Tmp2(3), Tmp3(3)

      ! Initialise flags
      IAPEX = 0
      ITENS = 0
      TauMax = 1.0d0
      IPL = -1

      ! Read material parameters
      !  PROPS(1) : G       Shear modulus (MPa)
      !  PROPS(2) : ENU     Poisson's ratio
      !  PROPS(3) : SIGCI   Intact rock UCS (MPa, positive)
      !  PROPS(4) : XMI     Intact rock constant mi
      !  PROPS(5) : GSI     Geological Strength Index
      !  PROPS(6) : DISTD   Disturbance factor D
      !  PROPS(7) : SPSI    sin(dilation angle)
      !  PROPS(8) : TENS    Tensile strength cutoff (MPa, positive; 0=auto)

      IntGlo = NPT
      G      = PROPS(1)
      ENU    = PROPS(2)
      VNU    = PROPS(2)
      SIGCI  = PROPS(3)
      XMI    = PROPS(4)
      GSI    = PROPS(5)
      DISTD  = PROPS(6)
      SPSI   = PROPS(7)
      TENS   = PROPS(8)

      one = 1.0d0
      two = 2.0d0

      ! ----------------------------------------------------------------
      ! Compute Hoek-Brown rock-mass constants
      !   mb = mi * exp( (GSI - 100) / (28 - 14*D) )
      !   s  = exp( (GSI - 100) / (9  -  3*D) )
      !   a  = 0.5 + (1/6)*( exp(-GSI/15) - exp(-20/3) )
      ! ----------------------------------------------------------------
      XMB  = XMI  * exp( (GSI - 100.0d0) / (28.0d0 - 14.0d0*DISTD) )
      XS   = exp( (GSI - 100.0d0) / (9.0d0 - 3.0d0*DISTD) )
      XA   = 0.5d0 + (one/6.0d0) * ( exp(-GSI/15.0d0) - exp(-20.0d0/3.0d0) )

      ! ----------------------------------------------------------------
      ! Auto tensile strength: sigma_t = s * sigci / mb   (compression-positive sign: positive = tension)
      ! If user supplies TENS > 0, it overrides the auto value.
      ! ----------------------------------------------------------------
      if (TENS <= 0.0d0) then
        TENS = XS * SIGCI / XMB
      end if

      ! ----------------------------------------------------------------
      ! Elastic stiffness constants
      ! ----------------------------------------------------------------
      FAC = two * G / (one - two * ENU)
      D1  = FAC * (one - ENU)
      D2  = FAC * ENU

      ! ----------------------------------------------------------------
      ! Elastic stress increment
      ! ----------------------------------------------------------------
      DSTRANVOL = DSTRAN(1) + DSTRAN(2) + DSTRAN(3)
      dSig(1) = (D1 - D2) * DSTRAN(1) + D2 * DSTRANVOL
      dSig(2) = (D1 - D2) * DSTRAN(2) + D2 * DSTRANVOL
      dSig(3) = (D1 - D2) * DSTRAN(3) + D2 * DSTRANVOL
      dSig(4) = G * DSTRAN(4)
      if (NTENS == 6) then
        dSig(5) = G * DSTRAN(5)
        dSig(6) = G * DSTRAN(6)
      end if

      ! Elastic trial stress (tension-positive, Anura3D convention)
      SigE  = STRESS + dSig
      STRESS = SigE
      SigEQ  = SigE
      SigC   = SigE

      ! Elastic tangent stiffness matrix
      DDSDDE = 0.0d0
      DDSDDE(1:3,1:3) = D2
      DDSDDE(1,1) = D1
      DDSDDE(2,2) = D1
      DDSDDE(3,3) = D1
      DDSDDE(4,4) = G
      if (NTENS == 6) then
        DDSDDE(5,5) = G
        DDSDDE(6,6) = G
      end if

      ! ----------------------------------------------------------------
      ! Principal stress decomposition  (tension-positive sign)
      ! Note: PrnSig called WITHOUT NTENS (built-in version)
      !       Ordering: Sig1 <= Sig2 <= Sig3 (most negative = most compressive)
      ! ----------------------------------------------------------------
      call PrnSig(1, SigE, xN1, xN2, xN3, Sig1, Sig2, Sig3, P, Q)

      ! ----------------------------------------------------------------
      ! Convert to compression-positive for HB yield surface evaluation
      !   S1C = -Sig3  (most compressive principal stress, positive value)
      !   S3C = -Sig1  (least compressive principal stress, positive value)
      ! ----------------------------------------------------------------
      S1C = -Sig3
      S3C = -Sig1

      ! ----------------------------------------------------------------
      ! HB yield function (compression-positive):
      !   F_HB = S1C - S3C - SIGCI * (mb * S3C / SIGCI + s)^a
      ! Tension cutoff (tension-positive convention):
      !   FT3  = Sig1 - TENS   (Sig1 is most tensile = largest value)
      ! ----------------------------------------------------------------
      TENS_ARG = XMB * S3C / SIGCI + XS
      if (TENS_ARG < 0.0d0) TENS_ARG = 0.0d0
      F_HB = S1C - S3C - SIGCI * TENS_ARG**XA

      FT1 = SIG1 - TENS
      FT2 = SIG2 - TENS
      FT3 = SIG3 - TENS

      ! ----------------------------------------------------------------
      ! Check for yielding
      ! ----------------------------------------------------------------
      IF (F_HB > 1.0d-12 .OR. FT3 > 1.0d-12) GOTO 240

      ! Elastic: no correction needed
      IF (IPL == 0) GOTO 360
      IPL = 0
      GOTO 360

  240 IPL = 1

      ! ----------------------------------------------------------------
      ! Compute instantaneous Mohr-Coulomb equivalent parameters from
      ! the HB tangent at the current confining stress S3C.
      !
      ! From Hoek & Brown (2002) / Hoek, Carranza-Torres & Corkum (2002):
      !
      !   dS1C/dS3C = 1 + a * mb * (mb*S3C/sigci + s)^(a-1)
      !
      ! Instantaneous friction angle phi_eq:
      !   sin(phi_eq) = (dS1C/dS3C - 1) / (dS1C/dS3C + 1)
      !
      ! Instantaneous cohesion c_eq:
      !   c_eq = (S1C - (1+sin(phi_eq))/(1-sin(phi_eq)) * S3C) *
      !          cos(phi_eq) / 2
      !   or equivalently from the MC tangent intercept.
      ! ----------------------------------------------------------------
      TENS_ARG = XMB * S3C / SIGCI + XS
      if (TENS_ARG < 0.0d0) TENS_ARG = 0.0d0

      if (TENS_ARG > 0.0d0) then
        DERIV = one + XA * XMB * TENS_ARG**(XA - one)
      else
        DERIV = one
      end if

      ! Clamp derivative to avoid singular friction angle
      if (DERIV < 1.0d0 + 1.0d-10) DERIV = 1.0d0 + 1.0d-10

      SPHI = (DERIV - one) / (DERIV + one)
      if (SPHI > 0.999d0) SPHI = 0.999d0
      CPHI = sqrt(one - SPHI*SPHI)

      ! S1C_yield = S3C + SIGCI * TENS_ARG^a  (point on HB envelope)
      S1C_YIELD = S3C + SIGCI * TENS_ARG**XA

      ! MC equivalent cohesion from the tangent line:
      !   S1C = [(1+SPHI)/(1-SPHI)] * S3C + 2*c/(1-SPHI)
      ! => 2*c = (S1C_yield - DERIV * S3C) * (1-SPHI)
      COHS = 0.5d0 * (S1C_YIELD - DERIV * S3C) * (one - SPHI)
      if (COHS < 0.0d0) COHS = 0.0d0

      ! ----------------------------------------------------------------
      ! With equivalent MC parameters (SPHI, COHS, SPSI, TENS),
      ! re-evaluate MC yield functions in tension-positive convention.
      ! (MC: F = 0.5*(S3-S1) + 0.5*(S3+S1)*SPHI - c = 0,
      !       where S1 <= S2 <= S3 in tension-positive ordering)
      ! ----------------------------------------------------------------
      F21 = 0.5d0*(SIG2-SIG1) + 0.5d0*(SIG2+SIG1)*SPHI - COHS
      F32 = 0.5d0*(SIG3-SIG2) + 0.5d0*(SIG3+SIG2)*SPHI - COHS
      F31 = 0.5d0*(SIG3-SIG1) + 0.5d0*(SIG3+SIG1)*SPHI - COHS

      call PrnSig(0, SigEQ, Tmp1, Tmp2, Tmp3, SigV1, SigV2, SigV3, Dum1, Dum2)

      FF1 = 0.5d0*    (SIGV2-SIGV1) + 0.5d0*(SIGV2+SIGV1)*SPHI - COHS
      FF2 = 0.5d0*DABS(SIGV2-SIGV3) + 0.5d0*(SIGV2+SIGV3)*SPHI - COHS
      FF3 = 0.5d0*DABS(SIGV3-SIGV1) + 0.5d0*(SIGV3+SIGV1)*SPHI - COHS

      IF (FF1 > 1.0d-12 .OR. FF2 > 1.0d-12 .OR. FF3 > 1.0d-12) GOTO 250

      ! Plastic negative points
      IPL = -1

  250 CONTINUE
      VNUI0 = VNU
      VNUI1 = one - two*VNUI0
      VNUI2 = one + two*VNUI0
      VNUIQ = VNUI2 / VNUI1

      DUM    = SPSI / VNUI1
      PSIMIN = G * (-one + DUM)
      PSIMET = G * ( one + DUM)
      PSINU  = two * G * VNUI0 * DUM

      HA  = ( one - SPHI + SPSI - SPHI*SPSI) * SIG1 + &
            (-two - two/VNUI1*SPHI*SPSI) * SIG2 + &
            ( one - SPHI - SPSI + VNUIQ*SPHI*SPSI) * SIG3 + &
            two*(one + SPSI) * COHS

      HB  = ( one + SPHI + SPSI + VNUIQ*SPHI*SPSI) * SIG1 - &
            ( two + two/VNUI1*SPHI*SPSI) * SIG2 + &
            ( one + SPHI - SPSI - SPHI*SPSI) * SIG3 - &
            two*(one - SPSI) * COHS

      HAB = (VNUI1+SPSI)*SIG1 + &
            (VNUI1-SPSI)*(SIG3-TENS) - &
            (VNUI1+SPSI)*(TENS*(one+SPHI)-two*COHS)/(one-SPHI)

      HBA = (one-VNUI0)*SIG1 - VNUI0*SIG3 - &
            (one-VNUI0)*(TENS*(one+SPHI)-two*COHS)/(one-SPHI) + &
            VNUI0*TENS

      HAO = (one-VNUI0)*SIG2 - VNUI0*SIG3 - VNUI1*TENS

      HOC = (one-VNUI0)*SIG1 - VNUI0*SIG3 - VNUI1*TENS

      HAA = (VNUI1+VNUI2*SPSI)*(SIG1-(TENS*(one+SPHI)-two*COHS)/(one-SPHI)) + &
            (VNUI1-SPSI)*(SIG2-TENS) + (VNUI1-SPSI)*(SIG3-TENS)

      HAAB = VNUI0*(SIG1-(TENS*(one+SPHI)-two*COHS)/(one-SPHI)) - &
             (SIG2-TENS) + VNUI0*(SIG3-TENS)

      HAAO = (SIG1-(TENS*(one+SPHI)-two*COHS)/(one-SPHI)) - &
             VNUI0*(SIG2-TENS) - VNUI0*(SIG3-TENS)

      HBB  = (VNUI1+SPSI)*(SIG1-(TENS*(one+SPHI)-two*COHS)/(one-SPHI)) + &
             (VNUI1+SPSI)*(SIG2-(TENS*(one+SPHI)-two*COHS)/(one-SPHI)) + &
             (VNUI1-VNUI2*SPSI)*(SIG3-TENS)

      ! ----------------------------------------------------------------
      ! Determine return area (same logic as MC)
      ! ----------------------------------------------------------------
      IAREA  = 0
      IASIGN = 0

      IF (F31 > 0.0d0 .AND. HA >= 0.0d0 .AND. HB < 0.0d0 .AND. HAB < 0.0d0) THEN
        IAREA  = 2
        IASIGN = IASIGN + 1
      END IF

      IF (F31 > 0.0d0 .AND. HB >= 0.0d0) THEN
        IASIGN = IASIGN + 1
        IF (HBB < 0.0d0) THEN
          IAREA = 1
          GOTO 260
        ELSE
          IAREA = 8
        END IF
      END IF

      IF (F31 > 0.0d0 .AND. HA < 0.0d0) THEN
        IASIGN = IASIGN + 1
        IF (HAA < 0.0d0) THEN
          IAREA = 3
        ELSE
          IAREA = 7
        END IF
      END IF

      IF (FT3 > 1.0d-12 .AND. &
          HBA >= 0.0d0   .AND. &
          HAO <  0.0d0   .AND. &
          HOC <  0.0d0 ) THEN
        IAREA  = 4
        IASIGN = IASIGN + 1
      END IF

      IF (FT3 > 1.0d-12 .AND. HAB >= 0.0d0 .AND. HBA < 0.0d0) THEN
        IASIGN = IASIGN + 1
        IF (HAAB > 0.0d0) THEN
          IAREA = 5
        ELSE
          IF (IAREA == 0) IAREA = 7
        END IF
      END IF

      IF (FT3 > 1.0d-12 .AND. HAO >= 0.0d0) THEN
        IASIGN = IASIGN + 1
        IF (HAAO >= 0.0d0) THEN
          IAREA = 6
        ELSE
          IF (IAREA == 0) IAREA = 7
        END IF
      END IF

  260 CONTINUE

      DSP1 = 0.0d0
      DSP2 = 0.0d0
      DSP3 = 0.0d0

      ! ----------------------------------------------------------------
      ! Return mapping (identical to MC, using equivalent parameters)
      ! ----------------------------------------------------------------

      ! Area 1: Extension point
      IF (IAREA == 1) THEN
        A11   = G * (one + SPHI*SPSI/VNUI1)
        A12   = 0.5d0*G * (one + SPHI + SPSI + VNUIQ*SPHI*SPSI)
        DETER = A11*A11 - A12*A12
        RLAM31 = (F31*A11 - F32*A12) / DETER
        RLAM32 = (F32*A11 - F31*A12) / DETER
        RLAM21 = 0.0d0
        DSP1 = RLAM31*PSIMIN + RLAM32*PSINU
        DSP2 = RLAM31*PSINU  + RLAM32*PSIMIN
        DSP3 = RLAM31*PSIMET + RLAM32*PSIMET
      END IF

      ! Area 2: Regular yield surface
      IF (IAREA == 2) THEN
        A11    = G * (one + SPHI*SPSI/VNUI1)
        RLAM31 = F31 / A11
        RLAM21 = 0.0d0
        RLAM32 = 0.0d0
        DSP1 = RLAM31*PSIMIN
        DSP2 = RLAM31*PSINU
        DSP3 = RLAM31*PSIMET
      END IF

      ! Area 3: Compression point
      IF (IAREA == 3) THEN
        A11   = G * (one + SPHI*SPSI/VNUI1)
        A12   = 0.5d0*G * (one - SPHI - SPSI + VNUIQ*SPHI*SPSI)
        DETER = A11*A11 - A12*A12
        RLAM31 = (F31*A11 - F21*A12) / DETER
        RLAM21 = (F21*A11 - F31*A12) / DETER
        RLAM32 = 0.0d0
        DSP1 = RLAM31*PSIMIN + RLAM21*PSIMIN
        DSP2 = RLAM31*PSINU  + RLAM21*PSIMET
        DSP3 = RLAM31*PSIMET + RLAM21*PSINU
      END IF

      ! Area 4: Pure tension cutoff
      IF (IAREA == 4) THEN
        DUM    = two*G / VNUI1
        RLAMT3 = FT3 / (DUM*(one - VNUI0))
        DSP1 = DUM*RLAMT3*VNUI0
        DSP2 = DSP1
        DSP3 = DUM*RLAMT3*(one - VNUI0)
        ITENS = 1
        IPL   = 2
      END IF

      ! Area 5: Shear + tension cutoff combined
      IF (IAREA == 5) THEN
        A11   = VNUI1 + SPHI*SPSI
        A12   = VNUI1 + SPHI
        A21   = VNUI1 + SPSI
        A22   = two*(one - VNUI0)
        DUM   = G / VNUI1
        DETER = A11*A22 - A12*A21
        RLAM31 = (F31*A22 - FT3*A12) / DETER / DUM
        RLAMT3 = (FT3*A11 - F31*A21) / DETER / DUM
        DSP1 = DUM*(RLAM31*(-VNUI1+SPSI) + RLAMT3*two*VNUI0)
        DSP2 = DUM*(RLAM31*two*VNUI0*SPSI + RLAMT3*two*VNUI0)
        DSP3 = DUM*(RLAM31*(VNUI1+SPSI) + RLAMT3*two*(one-VNUI0))
        ITENS = 1
        IPL   = 2
      END IF

      ! Area 6: Biaxial tension cutoff
      IF (IAREA == 6) THEN
        RLAMT2 = (FT2*(one-VNUI0) - FT3*VNUI0) / (two*G)
        RLAMT3 = (FT3*(one-VNUI0) - FT2*VNUI0) / (two*G)
        DUM  = two*G / VNUI1
        DSP1 = DUM*VNUI0*(RLAMT2 + RLAMT3)
        DSP2 = DUM*(RLAMT2*(one-VNUI0) + RLAMT3*VNUI0)
        DSP3 = DUM*(RLAMT2*VNUI0       + RLAMT3*(one-VNUI0))
        ITENS = 1
        IPL   = 2
      END IF

      ! Area 7: Apex – triaxial tension cutoff
      IF (IAREA == 7) THEN
        DSP1 = SIG1 - (TENS*(one+SPHI)-two*COHS)/(one-SPHI)
        DSP2 = SIG2 - TENS
        DSP3 = SIG3 - TENS
        ITENS = 1
        IPL   = 2
      END IF

      ! Area 8: Extension apex
      IF (IAREA == 8) THEN
        DSP1 = SIG1 - (TENS*(one+SPHI)-two*COHS)/(one-SPHI)
        DSP2 = SIG2 - (SIG1-DSP1)
        DSP3 = SIG3 - TENS
        ITENS = 1
        IPL   = 2
      END IF

      ! Area 9: Full apex
      IF (IAREA == 9) THEN
        DSP1 = SIG1 - TENS
        DSP2 = SIG2 - TENS
        DSP3 = SIG3 - TENS
        IAPEX = 1
        IPL   = 2
      END IF

      ! ----------------------------------------------------------------
      ! First check: verify corrected stresses satisfy yield surface
      ! ----------------------------------------------------------------
      SIG1R = SIG1 - DSP1
      SIG2R = SIG2 - DSP2
      SIG3R = SIG3 - DSP3

      F21R = 0.5d0*(SIG2R-SIG1R) + 0.5d0*(SIG2R+SIG1R)*SPHI - COHS
      F32R = 0.5d0*(SIG3R-SIG2R) + 0.5d0*(SIG3R+SIG2R)*SPHI - COHS
      F31R = 0.5d0*(SIG3R-SIG1R) + 0.5d0*(SIG3R+SIG1R)*SPHI - COHS
      FT1R = SIG1R - TENS
      FT2R = SIG2R - TENS
      FT3R = SIG3R - TENS

      IF ( F21R > 1.0d-6 .OR. &
           F32R > 1.0d-6 .OR. &
           F31R > 1.0d-6 .OR. &
           FT1R > 1.0d-6 .OR. &
           FT2R > 1.0d-6 .OR. &
           FT3R > 1.0d-6 ) THEN
        ICREC = 0

        IF (IAREA == 5 .AND. F32R > 1.0d-6) THEN
          ! Correct to zone 8
          ICREC = ICREC + 1
          DSP1 = SIG1 - (TENS*(one+SPHI)-two*COHS)/(one-SPHI)
          DSP2 = SIG2 - (SIG1-DSP1)
          DSP3 = SIG3 - TENS
        END IF

        IF (IAREA == 6 .AND. FT1R > 1.0d-6) THEN
          ! Correct to zone 9
          ICREC = ICREC + 1
          DSP1 = SIG1 - TENS
          DSP2 = SIG2 - TENS
          DSP3 = SIG3 - TENS
          ITENS = 0
          IAPEX = 1
        END IF
      END IF

      ! ----------------------------------------------------------------
      ! Second check (final validation)
      ! ----------------------------------------------------------------
      SIG1R = SIG1 - DSP1
      SIG2R = SIG2 - DSP2
      SIG3R = SIG3 - DSP3

      F21 = 0.5d0*(SIG2R-SIG1R) + 0.5d0*(SIG2R+SIG1R)*SPHI - COHS
      F32 = 0.5d0*(SIG3R-SIG2R) + 0.5d0*(SIG3R+SIG2R)*SPHI - COHS
      F31 = 0.5d0*(SIG3R-SIG1R) + 0.5d0*(SIG3R+SIG1R)*SPHI - COHS
      FT1 = SIG1R - TENS
      FT2 = SIG2R - TENS
      FT3 = SIG3R - TENS

      ! ----------------------------------------------------------------
      ! Reconstruct Cartesian stresses
      ! ----------------------------------------------------------------
      call CarSig(Sig1R, Sig2R, Sig3R, xN1, xN2, xN3, NTENS, SigC)
      STRESS = SigC

      ! TAUMAX — used to assess accuracy at plastic points
      TAUMAX = 0.5d0*(SIG3-DSP3 - SIG1+DSP1)
      TAUMAX = MAX(TAUMAX, COHS, 0.5d0)

      ! ----------------------------------------------------------------
      ! Update state variables
      !   STATEV(1) : accumulated equivalent plastic strain
      !   STATEV(2) : yield flag (0=elastic, 1=shear, 2=tension)
      !   STATEV(3) : equivalent friction angle (degrees)
      !   STATEV(4) : equivalent cohesion (MPa)
      ! ----------------------------------------------------------------
      DEPS_EQ = SQRT( (two/3.0d0) * (DSP1**2 + DSP2**2 + DSP3**2) )
      STATEV(1) = STATEV(1) + DEPS_EQ

      if (IPL == 1) then
        STATEV(2) = 1.0d0  ! shear yield
      else if (IPL == 2) then
        STATEV(2) = 2.0d0  ! tension cutoff
      else
        STATEV(2) = 0.0d0  ! elastic
      end if

      ! Convert sin(phi) to friction angle in degrees
      if (SPHI > 0.0d0 .AND. SPHI < one) then
        STATEV(3) = ASIN(SPHI) * 180.0d0 / 3.14159265358979d0
      else
        STATEV(3) = 0.0d0
      end if

      STATEV(4) = COHS

  360 CONTINUE

      RETURN

      End Subroutine UMAT_HB

!***********************************************************************
      Subroutine MatTranspose(A, Ia, AT, IAt, N1, N2)
!  Transposes matrix A(N1,N2) into AT(N2,N1)
      implicit none
      integer :: Ia, N1, N2, IAt
      double precision :: A(Ia,*), AT(IAt,*)
      integer :: I, J

      Do I = 1, N1
        Do J = 1, N2
          AT(I,J) = A(J,I)
        End Do
      End Do

      End Subroutine MatTranspose

!***********************************************************************
      Subroutine MatMat(A, Ia, B, Ib, N1, N2, N3, C, Ic)
!  Matrix multiplication C(N1,N3) = A(N1,N2) * B(N2,N3)
      implicit none
      integer :: Ia, Ib, Ic, N1, N2, N3
      double precision :: A(Ia,*), B(Ib,*), C(Ic,*)
      integer :: I, J, K
      double precision :: Sum

      Do I = 1, N1
        Do J = 1, N3
          Sum = 0.0d0
          Do K = 1, N2
            Sum = Sum + A(I,K) * B(K,J)
          End Do
          C(I,J) = Sum
        End Do
      End Do

      End Subroutine MatMat

!***********************************************************************
      Subroutine CarSig(S1, S2, S3, xN1, xN2, xN3, ntens, Stress)
!  Returns Cartesian stresses from principal stresses and directions.
!
!  S1, S2, S3      I   R     principal stresses
!  xN1, xN2, xN3   I   R(3)  principal directions
!  ntens           I   I     stress vector size (4 or 6)
!  Stress          O   R(ntens) Cartesian stress vector
!
      implicit none

      double precision, intent(in)  :: S1, S2, S3
      double precision, intent(in),  dimension(3) :: xN1, xN2, xN3
      integer,          intent(in)  :: ntens
      double precision, intent(out), dimension(ntens) :: Stress

      integer :: I, IDim
      double precision, dimension(:,:), allocatable :: SM, T, TT, STT

      IDim = 3
      allocate(SM(IDim,IDim), T(IDim,IDim), TT(IDim,IDim), STT(IDim,IDim))

      ! Assemble transformation matrix T whose columns are the principal directions
      Do I = 1, 3
        T(I,1)  = xN1(I)
        T(I,2)  = xN2(I)
        T(I,3)  = xN3(I)
        TT(1,I) = T(I,1)
        TT(2,I) = T(I,2)
        TT(3,I) = T(I,3)
      End Do

      ! Diagonal principal stress matrix
      SM     = 0.0d0
      SM(1,1) = S1
      SM(2,2) = S2
      SM(3,3) = S3

      ! Cartesian stress: T * S * T^T
      Call MatMat(SM, IDim, TT,  IDim, IDim, IDim, IDim, STT, IDim)
      Call MatMat(T,  IDim, STT, IDim, IDim, IDim, IDim, SM,  IDim)

      Do I = 1, IDim
        Stress(I) = SM(I,I)
      End Do
      Stress(4) = SM(2,1)
      if (ntens == 6) then
        Stress(5) = SM(3,2)
        Stress(6) = SM(3,1)
      end if

      End Subroutine CarSig

    end module ModHoekBrown
