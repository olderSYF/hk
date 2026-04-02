! =============================================================================
!  A3DMohrCoulombStandard.f
!
!  Anura3D 2025 – Hoek-Brown constitutive model
!  embedded via the "Mohr-Coulomb shell" approach.
!
!  Structure:
!    module ModMohrCoulomb
!      ESM_MC       – wrapper subroutine (UNCHANGED from Anura3D original)
!      UMAT         – replaced with 2-D/3-D Hoek-Brown implementation
!      CalcPrincipalHB – principal stress + direction subroutine
!      Solve3x3HB   – 3x3 linear system solver (Gaussian elimination)
!      MatTranspose – matrix transpose utility (unchanged)
!
!  PROPS mapping (借壳 Mohr-Coulomb / "Shell" parameter mapping):
!    PROPS(1)  G        shear modulus [kPa]   (set by ExternalSoilModel.for)
!    PROPS(2)  ENU      Poisson ratio
!    PROPS(3)  SCI      sigma_ci [kPa]        (GiD "Friction Angle" field, raw)
!    PROPS(4)  GSI_val  Geological Strength Index  (GiD "Cohesion" field, raw)
!    PROPS(5)  AMI      m_i intact-rock constant   (GiD "Dilatancy Angle" field, raw)
!    PROPS(6)  DDF      D  disturbance factor      (GiD "Tensile Strength" field)
!
!  Stress sign convention: compression POSITIVE (Anura3D standard)
!  Stress/strain tensor order (NTENS=4, 2-D plane-strain):
!    (1) sigma_xx  (2) sigma_yy  (3) sigma_zz  (4) sigma_xy
!  Stress/strain tensor order (NTENS=6, 3-D):
!    (1) sigma_xx  (2) sigma_yy  (3) sigma_zz
!    (4) sigma_xy  (5) sigma_yz  (6) sigma_xz
!
!  State variables:
!    STATEV(1)  cumulative equivalent plastic strain
!    STATEV(2)  plastic indicator (0=elastic, 1=HB yield, 2=tensile cutoff)
!
!  Hoek-Brown parameters derived from PROPS:
!    mb  = m_i * exp( (GSI-100) / (28 - 14*D) )
!    s   = exp( (GSI-100) / (9 - 3*D) )
!    a   = 0.5 + (1/6)*( exp(-GSI/15) - exp(-20/3) )
!    sigma_t (tensile strength) = s * sigma_ci / mb
!
!  Yield function (compression-positive):
!    F = sigma_1 - sigma_3 - sigma_ci * (mb*sigma_3/sigma_ci + s)^a
!    Tensile cutoff: sigma_3 >= -sigma_t
! =============================================================================

    module ModMohrCoulomb
    contains 
    
! =============================================================================
!  ESM_MC  –  WRAPPER SUBROUTINE (UNCHANGED FROM ANURA3D ORIGINAL)
!
!  This subroutine is kept exactly as supplied by Anura3D so that the
!  module interface remains fully compatible with the solver.
! =============================================================================
    Subroutine ESM_MC(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER, &
     DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)


      implicit double precision (a-h, o-z) 
      integer :: NTENS, NSTATEV, NADDVAR, NPROPS, NPT, NOEL, IDSET, NUMBEROFPHASES
      double precision :: EUNLOADING, PLASTICMULTIPLIER
      CHARACTER*80 CMNAME   
      DIMENSION STRESS(NTENS), DSTRAN(NTENS),STATEV(NSTATEV),ADDITIONALVAR(NADDVAR),PROPS(NPROPS)
      

!---Local variables required in standard UMAT
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
        double precision :: Value
        double precision :: Porosity, WaterPressure, WaterPressure0
        double precision :: GasPressure, GasPressure0, DegreeSaturation  
              
        integer :: ndi, nshr, layer, kspt, kstep, kinc     

        allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
              coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )
    
!---Initialization
        Eunloading = 0.0
        PlasticMultiplier = 0.0
          
!---Rename additional variables
        Porosity         = AdditionalVar(1)
        WaterPressure    = AdditionalVar(2)
        WaterPressure0   = AdditionalVar(3)
        GasPressure      = AdditionalVar(4)
        GasPressure0     = AdditionalVar(5)
        DegreeSaturation = AdditionalVar(6)
        time(1) = AdditionalVar(7)   !TotalRealTime
        time(2) = AdditionalVar(8)   !OverallTotalTime
        dtime   = AdditionalVar(9)   !TimeIncrement
        IStep   = AdditionalVar(10)    
        TimeStep = AdditionalVar(11)  !Note: Very first time and load step: Istep=1 and TimeStep=1

!---Call the UMAT (Hoek-Brown implementation below)
          call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, &
           dfgrd1, noel, npt, layer, kspt, kstep, kinc)

!---Definition of Eunloading -> required to define the max time step
      Eunloading = max(ddsdde(1,1), ddsdde(2,2), ddsdde(3,3))

!---PlasticMultiplier: pass through plastic indicator stored in STATEV(2)
!   (enables plastic point marking via SetIPL in ExternalSoilModel.for)
      PlasticMultiplier = statev(2)

        return

    end subroutine ESM_MC

! =============================================================================
!  UMAT  –  2-D / 3-D HOEK-BROWN CONSTITUTIVE MODEL
!
!  Implements the Hoek-Brown failure criterion with:
!    - Elastic predictor / plastic corrector (return-mapping)
!    - Newton-Raphson iteration for HB surface (3x3 system)
!    - Tensile cutoff with apex return
!    - Elastic tangent stiffness (DDSDDE)
!
!  References:
!    Hoek & Brown (1980), Hoek, Carranza-Torres & Corkum (2002)
!    de Souza Neto, Peric & Owen (2008) – return mapping algorithm
! =============================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
      RPL,DDSDDT,DRPLDE,DRPLDT, &
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      implicit double precision (a-h, o-z) 

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS), STATEV(NSTATEV), &
                DDSDDE(NTENS,NTENS), DDSDDT(NTENS), DRPLDE(NTENS), &
                STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1), &
                PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)

      ! ---- Local arrays ----
      DIMENSION dSig(6), SigE(6), SigC(6)
      DIMENSION xN1(3), xN2(3), xN3(3)

      ! J_A starts with 'J' (i-n range = implicit integer), so declare explicitly
      double precision :: J_A(3,3), RHS_A(3), dX_A(3)

      integer :: IPL, ITER, INFO_3x3
      integer, parameter :: MAXITER = 100

      ! ------------------------------------------------------------------
      ! 1.  Read material parameters
      ! ------------------------------------------------------------------
      ! PROPS(1)  G        shear modulus [kPa]
      ! PROPS(2)  ENU      Poisson ratio
      ! PROPS(3)  SCI      sigma_ci [kPa]  (HB uniaxial compressive strength)
      ! PROPS(4)  GSI_val  Geological Strength Index (dimensionless, 0-100)
      ! PROPS(5)  AMI      m_i (intact-rock Hoek-Brown constant)
      ! PROPS(6)  DDF      D   (disturbance factor, 0-1)

      GG      = PROPS(1)   ! shear modulus G
      ENU     = PROPS(2)   ! Poisson ratio  nu
      SCI     = PROPS(3)   ! sigma_ci
      GSI_val = PROPS(4)   ! GSI
      AMI     = PROPS(5)   ! m_i
      DDF     = PROPS(6)   ! D

      ! Guard against degenerate inputs
      if (GG  <= 0.0d0) GG  = 1.0d0
      if (SCI <= 0.0d0) SCI = 1.0d0
      if (ENU >= 0.5d0) ENU = 0.49d0
      if (ENU <= 0.0d0) ENU = 0.001d0
      if (DDF < 0.0d0)  DDF = 0.0d0
      if (DDF > 1.0d0)  DDF = 1.0d0

      ! ------------------------------------------------------------------
      ! 2.  Hoek-Brown derived parameters
      !     mb = m_i * exp( (GSI-100) / (28 - 14*D) )
      !     s  = exp( (GSI-100) / (9 - 3*D) )
      !     a  = 0.5 + (1/6)*( exp(-GSI/15) - exp(-20/3) )
      ! ------------------------------------------------------------------
      AMB = AMI * dexp((GSI_val - 100.0d0) / (28.0d0 - 14.0d0*DDF))
      SVL = dexp((GSI_val - 100.0d0) / (9.0d0 - 3.0d0*DDF))
      AVL = 0.5d0 + (1.0d0/6.0d0) * &
            (dexp(-GSI_val/15.0d0) - dexp(-20.0d0/3.0d0))

      if (AMB < 1.0d-12) AMB = 1.0d-12

      ! Tensile strength magnitude: sigma_t = s*sigma_ci/mb (>0)
      TENS_HB = SVL * SCI / AMB

      ! ------------------------------------------------------------------
      ! 3.  Elastic moduli
      !     D1 = lambda + 2*mu = 2G(1-nu)/(1-2nu)
      !     D2 = lambda        = 2G*nu/(1-2nu)
      ! ------------------------------------------------------------------
      FAC = 2.0d0 * GG / (1.0d0 - 2.0d0*ENU)
      D1  = FAC * (1.0d0 - ENU)
      D2  = FAC * ENU

      ! ------------------------------------------------------------------
      ! 4.  Elastic stiffness matrix DDSDDE
      ! ------------------------------------------------------------------
      DDSDDE = 0.0d0
      DDSDDE(1,1) = D1;  DDSDDE(2,2) = D1;  DDSDDE(3,3) = D1
      DDSDDE(1,2) = D2;  DDSDDE(2,1) = D2
      DDSDDE(1,3) = D2;  DDSDDE(3,1) = D2
      DDSDDE(2,3) = D2;  DDSDDE(3,2) = D2
      DDSDDE(4,4) = GG
      if (NTENS == 6) then
        DDSDDE(5,5) = GG
        DDSDDE(6,6) = GG
      end if

      ! ------------------------------------------------------------------
      ! 5.  Elastic trial stress predictor
      !     dSig = C_e : dEps
      ! ------------------------------------------------------------------
      DVOL = DSTRAN(1) + DSTRAN(2) + DSTRAN(3)
      dSig(1) = (D1 - D2)*DSTRAN(1) + D2*DVOL
      dSig(2) = (D1 - D2)*DSTRAN(2) + D2*DVOL
      dSig(3) = (D1 - D2)*DSTRAN(3) + D2*DVOL
      dSig(4) = GG * DSTRAN(4)
      if (NTENS == 6) then
        dSig(5) = GG * DSTRAN(5)
        dSig(6) = GG * DSTRAN(6)
      end if

      do I = 1, NTENS
        SigE(I) = STRESS(I) + dSig(I)
      end do

      ! ------------------------------------------------------------------
      ! 6.  Principal stresses and directions of the trial stress
      !     Sp1 >= Sp2 >= Sp3  (compression-positive ordering)
      !     xN1, xN2, xN3 : unit eigenvectors (3-component vectors)
      ! ------------------------------------------------------------------
      call CalcPrincipalHB(SigE, NTENS, Sp1, Sp2, Sp3, xN1, xN2, xN3)

      Sp1tr = Sp1;  Sp2tr = Sp2;  Sp3tr = Sp3

      ! ------------------------------------------------------------------
      ! 7.  Yield function check
      !     F = s1 - s3 - sci * (mb*s3/sci + s)^a
      !     Tensile cutoff: if s3 < -sigma_t
      ! ------------------------------------------------------------------
      PHItr = AMB * Sp3tr / SCI + SVL
      if (PHItr < 0.0d0) PHItr = 0.0d0

      if (PHItr > 0.0d0) then
        F_yield = Sp1tr - Sp3tr - SCI * PHItr**AVL
      else
        F_yield = Sp1tr - Sp3tr   ! degenerate case (phi<=0)
      end if

      ! F_tens > 0 means sigma_3 is more tensile than the HB apex
      F_tens = -Sp3tr - TENS_HB

      IPL  = 0
      DLAM = 0.0d0
      SigC(1:NTENS) = SigE(1:NTENS)    ! default: elastic step

      ! ------------------------------------------------------------------
      ! 8.  Plasticity check and return mapping
      ! ------------------------------------------------------------------
      if (F_yield <= 0.0d0 .and. F_tens <= 0.0d0) then
        ! ---  Purely elastic  ---
        goto 900
      end if

      IPL = 1

      ! ---- Case A: Tensile failure – return to apex ----
      ! When sigma_3 is more tensile than the HB apex (phi<=0), all
      ! principal stresses are returned to -sigma_t (compression-positive).
      if (F_tens > 0.0d0) then
        Sp3 = -TENS_HB
        Sp2 = min(Sp2tr, -TENS_HB)
        Sp1 = min(Sp1tr, -TENS_HB)
        IPL = 2
        DLAM = 0.0d0   ! not used for back-transform; delta already known
        goto 800
      end if

      ! ---- Case B: Hoek-Brown yield surface (sigma_3 >= -sigma_t) ----
      ! Newton-Raphson return mapping in principal stress space.
      ! Variables: x = [Sp1, Sp3, DLAM]
      !
      ! Update equations (associative flow rule):
      !   Sp1 = Sp1tr - DLAM * g1   where g1 = D1 - D2*(1+h)
      !   Sp3 = Sp3tr - DLAM * g3   where g3 = D2 - D1*(1+h)
      !   h   = a*mb*(mb*Sp3/sci+s)^(a-1)   (gradient component)
      !
      ! Residual equations:
      !   R1 = Sp1 - Sp1tr + DLAM*g1 = 0   (Sp1 update consistency)
      !   R3 = Sp3 - Sp3tr + DLAM*g3 = 0   (Sp3 update consistency)
      !   RF = Sp1 - Sp3 - sci*(mb*Sp3/sci+s)^a = 0  (yield condition)

      Sp1_k = Sp1tr
      Sp3_k = Sp3tr
      DLAM  = 0.0d0

      do ITER = 1, MAXITER

        ! Current phi, h (gradient), hp (derivative of h wrt Sp3)
        PHI_k = AMB * Sp3_k / SCI + SVL
        if (PHI_k < 1.0d-12) then
          ! Approaching tensile apex during iteration → switch to apex return
          Sp3 = -TENS_HB
          Sp2 = min(Sp2tr, -TENS_HB)
          Sp1 = min(Sp1tr, -TENS_HB)
          IPL = 2
          goto 800
        end if

        h_k = AVL * AMB * PHI_k**(AVL - 1.0d0)

        if (PHI_k > 1.0d-10) then
          hp_k = AVL * (AVL - 1.0d0) * AMB**2 / SCI * PHI_k**(AVL - 2.0d0)
        else
          hp_k = 0.0d0
        end if

        ! Stiffness-weighted gradient components
        !   n1 = dF/dSp1 = 1,  n3 = dF/dSp3 = -(1+h)
        !   g1 = D1*n1 + D2*n3 = D1 - D2*(1+h)
        !   g3 = D2*n1 + D1*n3 = D2 - D1*(1+h)
        g1_k = D1 - D2*(1.0d0 + h_k)
        g3_k = D2 - D1*(1.0d0 + h_k)

        ! Residuals
        R1_k = Sp1_k - Sp1tr + DLAM * g1_k
        R3_k = Sp3_k - Sp3tr + DLAM * g3_k
        RF_k = Sp1_k - Sp3_k - SCI * PHI_k**AVL

        ! Convergence check (relative + absolute tolerance)
        TOL_ABS = 1.0d-8 * (dabs(Sp1tr) + dabs(Sp3tr) + SCI) + 1.0d-12
        RNORM   = dabs(R1_k) + dabs(R3_k) + dabs(RF_k)
        if (RNORM <= TOL_ABS) exit

        ! 3x3 Jacobian  J = dR/dx
        !
        !   dR1/dSp1 = 1
        !   dR1/dSp3 = DLAM * dg1/dSp3 = -DLAM*D2*hp
        !   dR1/dDLAM = g1
        !   dR3/dSp1 = 0
        !   dR3/dSp3 = 1 + DLAM*dg3/dSp3 = 1 - DLAM*D1*hp
        !   dR3/dDLAM = g3
        !   dRF/dSp1 = 1
        !   dRF/dSp3 = -(1+h)
        !   dRF/dDLAM = 0
        J_A(1,1) = 1.0d0
        J_A(1,2) = -DLAM * D2 * hp_k
        J_A(1,3) = g1_k
        J_A(2,1) = 0.0d0
        J_A(2,2) = 1.0d0 - DLAM * D1 * hp_k
        J_A(2,3) = g3_k
        J_A(3,1) = 1.0d0
        J_A(3,2) = -(1.0d0 + h_k)
        J_A(3,3) = 0.0d0

        ! Right-hand side:  rhs = -[R1, R3, RF]
        RHS_A(1) = -R1_k
        RHS_A(2) = -R3_k
        RHS_A(3) = -RF_k

        ! Solve J * dx = rhs
        call Solve3x3HB(J_A, RHS_A, dX_A, INFO_3x3)

        if (INFO_3x3 /= 0) then
          ! Singular Jacobian – try a scalar update on DLAM
          DENOM_TMP = g1_k - g3_k - (1.0d0 + h_k)*g3_k
          if (dabs(DENOM_TMP) > 1.0d-20) then
            DLAM = DLAM + RF_k / DENOM_TMP
          end if
          exit
        end if

        Sp1_k = Sp1_k + dX_A(1)
        Sp3_k = Sp3_k + dX_A(2)
        DLAM  = DLAM  + dX_A(3)

        ! Prevent DLAM from going negative (physically inadmissible)
        if (DLAM < 0.0d0) DLAM = 0.0d0

      end do  ! Newton-Raphson loop

      ! Compute corrected Sp1 and Sp3 from final DLAM
      PHI_k = AMB * Sp3_k / SCI + SVL
      if (PHI_k < 1.0d-12) PHI_k = 1.0d-12
      h_k   = AVL * AMB * PHI_k**(AVL - 1.0d0)

      ! Corrected Sp2 (middle principal stress) from elastic relation:
      !   dSp2 = -DLAM * (D2*n1 + D2*n3) = DLAM*D2*h
      Sp1 = Sp1_k
      Sp2 = Sp2tr + DLAM * D2 * h_k
      Sp3 = Sp3_k

      ! If sigma_3 ended up below tensile limit, snap to apex
      if (Sp3 < -TENS_HB - 1.0d-10) then
        Sp3 = -TENS_HB
        Sp2 = min(Sp2, -TENS_HB)
        Sp1 = min(Sp1, -TENS_HB)
        IPL = 2
      end if

 800  continue
      ! ------------------------------------------------------------------
      ! 9.  Back-transform corrected principal stresses to Cartesian
      !     sigma_cart = Sp1*(xN1 x xN1) + Sp2*(xN2 x xN2) + Sp3*(xN3 x xN3)
      !     (dyadic / outer-product reconstruction)
      ! ------------------------------------------------------------------
      SigC(1) = Sp1*xN1(1)*xN1(1) + Sp2*xN2(1)*xN2(1) + Sp3*xN3(1)*xN3(1)
      SigC(2) = Sp1*xN1(2)*xN1(2) + Sp2*xN2(2)*xN2(2) + Sp3*xN3(2)*xN3(2)
      SigC(3) = Sp1*xN1(3)*xN1(3) + Sp2*xN2(3)*xN2(3) + Sp3*xN3(3)*xN3(3)
      SigC(4) = Sp1*xN1(1)*xN1(2) + Sp2*xN2(1)*xN2(2) + Sp3*xN3(1)*xN3(2)
      if (NTENS == 6) then
        SigC(5) = Sp1*xN1(2)*xN1(3) + Sp2*xN2(2)*xN2(3) + Sp3*xN3(2)*xN3(3)
        SigC(6) = Sp1*xN1(1)*xN1(3) + Sp2*xN2(1)*xN2(3) + Sp3*xN3(1)*xN3(3)
      end if

 900  continue
      ! ------------------------------------------------------------------
      ! 10.  Update STRESS
      ! ------------------------------------------------------------------
      do I = 1, NTENS
        STRESS(I) = SigC(I)
      end do

      ! ------------------------------------------------------------------
      ! 11.  Update state variables
      !        STATEV(1) : cumulative equivalent plastic strain
      !        STATEV(2) : plastic indicator (0=elastic, 1=HB, 2=tensile)
      ! ------------------------------------------------------------------
      if (IPL == 1 .and. DLAM > 0.0d0) then
        ! Equivalent plastic strain increment (von Mises-like scalar measure)
        ! Only updated for HB yield (IPL=1); tensile apex (IPL=2) skipped
        ! since DLAM is not defined via NR in that path.
        PHI_fin = AMB * Sp3 / SCI + SVL
        if (PHI_fin < 1.0d-12) PHI_fin = 1.0d-12
        h_fin = AVL * AMB * PHI_fin**(AVL - 1.0d0)
        DEpsP_eq = DLAM * dsqrt( 1.0d0 + (1.0d0 + h_fin)**2 )
        STATEV(1) = STATEV(1) + DEpsP_eq
      end if
      STATEV(2) = dble(IPL)

      RETURN

      END SUBROUTINE UMAT


! =============================================================================
!  CalcPrincipalHB
!
!  Compute eigenvalues (principal stresses) and eigenvectors (principal
!  directions) of the symmetric stress tensor.
!
!  For NTENS=4 (2-D plane strain):
!    - sigma_zz is always a principal stress (z-direction decoupled).
!    - The in-plane 2x2 block gives two more eigenvalues analytically.
!
!  For NTENS=6 (3-D):
!    - Uses the trigonometric solution of the characteristic cubic.
!    - Eigenvectors computed via null-space of (S - lambda*I).
!
!  Output:
!    Sp1 >= Sp2 >= Sp3  (compression-positive; largest = most compressive)
!    xN1, xN2, xN3     (corresponding unit eigenvectors in R^3)
! =============================================================================
      subroutine CalcPrincipalHB(Sig, NTENS, Sp1, Sp2, Sp3, xN1, xN2, xN3)

      implicit double precision (a-h, o-z)
      integer, intent(in)  :: NTENS
      double precision, intent(in)  :: Sig(NTENS)
      double precision, intent(out) :: Sp1, Sp2, Sp3
      double precision, intent(out) :: xN1(3), xN2(3), xN3(3)

      double precision :: sxx, syy, szz, sxy, syz, sxz
      double precision :: pmean, qrad, theta_p
      double precision :: SA, SB, SC
      double precision :: NAx, NAy, NBx, NBy
      ! 3-D eigenvalue variables
      double precision :: p1, q1, p2, r_val, phi_val, arg
      double precision :: eig1, eig2, eig3, tmp
      double precision :: vnorm
      integer :: i, j
      double precision :: eigvals(3)

      double precision, parameter :: PICONST = 3.14159265358979323846d0
      double precision, parameter :: SMALL   = 1.0d-14

      sxx = Sig(1); syy = Sig(2); szz = Sig(3); sxy = Sig(4)
      if (NTENS == 6) then
        syz = Sig(5); sxz = Sig(6)
      else
        syz = 0.0d0; sxz = 0.0d0
      end if

      ! ===================================================================
      ! 2-D plane-strain path  (xz and yz shear are both zero)
      ! ===================================================================
      if (sxz == 0.0d0 .and. syz == 0.0d0) then

        pmean = 0.5d0 * (sxx + syy)
        qrad  = dsqrt( (0.5d0*(sxx - syy))**2 + sxy**2 )

        SA = pmean + qrad   ! larger  in-plane principal stress
        SB = pmean - qrad   ! smaller in-plane principal stress
        SC = szz            ! out-of-plane principal stress (z is principal dir)

        ! In-plane principal directions
        if (qrad > SMALL) then
          theta_p = 0.5d0 * datan2(2.0d0*sxy, sxx - syy)
        else
          theta_p = 0.0d0
        end if
        NAx =  dcos(theta_p); NAy = dsin(theta_p)   ! direction of SA
        NBx = -dsin(theta_p); NBy = dcos(theta_p)   ! direction of SB
        ! (z-direction [0,0,1] is direction of SC)

        ! Sort SA, SB, SC in descending order → Sp1 >= Sp2 >= Sp3
        if (SA >= SB .and. SA >= SC) then
          Sp1 = SA
          xN1 = (/ NAx,  NAy,  0.0d0 /)
          if (SB >= SC) then
            Sp2 = SB; xN2 = (/ NBx, NBy, 0.0d0 /)
            Sp3 = SC; xN3 = (/ 0.0d0, 0.0d0, 1.0d0 /)
          else
            Sp2 = SC; xN2 = (/ 0.0d0, 0.0d0, 1.0d0 /)
            Sp3 = SB; xN3 = (/ NBx, NBy, 0.0d0 /)
          end if

        else if (SB >= SA .and. SB >= SC) then
          Sp1 = SB
          xN1 = (/ NBx, NBy, 0.0d0 /)
          if (SA >= SC) then
            Sp2 = SA; xN2 = (/ NAx, NAy, 0.0d0 /)
            Sp3 = SC; xN3 = (/ 0.0d0, 0.0d0, 1.0d0 /)
          else
            Sp2 = SC; xN2 = (/ 0.0d0, 0.0d0, 1.0d0 /)
            Sp3 = SA; xN3 = (/ NAx, NAy, 0.0d0 /)
          end if

        else
          ! SC is the largest
          Sp1 = SC
          xN1 = (/ 0.0d0, 0.0d0, 1.0d0 /)
          if (SA >= SB) then
            Sp2 = SA; xN2 = (/ NAx, NAy, 0.0d0 /)
            Sp3 = SB; xN3 = (/ NBx, NBy, 0.0d0 /)
          else
            Sp2 = SB; xN2 = (/ NBx, NBy, 0.0d0 /)
            Sp3 = SA; xN3 = (/ NAx, NAy, 0.0d0 /)
          end if
        end if

        return
      end if

      ! ===================================================================
      ! 3-D general path:
      !   Characteristic polynomial: lambda^3 - I1*lambda^2 + I2*lambda - I3 = 0
      !   Solved with the trigonometric (Cardano) method for real symmetric
      !   matrices (always three real eigenvalues).
      ! ===================================================================

      ! Stress invariants
      S_I1 = sxx + syy + szz
      S_I2 = sxx*syy + syy*szz + szz*sxx - sxy**2 - syz**2 - sxz**2
      S_I3 = sxx*syy*szz + 2.0d0*sxy*syz*sxz &
           - sxx*syz**2 - syy*sxz**2 - szz*sxy**2

      ! Deviatoric substitution  lambda = mu + I1/3
      p1    = S_I1 / 3.0d0
      q1    = S_I2 - S_I1**2 / 3.0d0
      r_val = 2.0d0*S_I1**3/27.0d0 - S_I1*S_I2/3.0d0 + S_I3

      ! Trigonometric solution (valid when q1 < 0)
      if (q1 < -SMALL) then
        p2      = dsqrt(-q1 / 3.0d0)
        arg     = r_val / (2.0d0 * p2**3)
        arg     = max(-1.0d0, min(1.0d0, arg))
        phi_val = dacos(arg) / 3.0d0
        eig1 = p1 + 2.0d0*p2*dcos(phi_val)
        eig2 = p1 + 2.0d0*p2*dcos(phi_val - 2.0d0*PICONST/3.0d0)
        eig3 = p1 + 2.0d0*p2*dcos(phi_val - 4.0d0*PICONST/3.0d0)
      else
        ! Near hydrostatic (degenerate): all eigenvalues ≈ I1/3
        eig1 = p1
        eig2 = p1
        eig3 = p1
      end if

      ! Sort eigenvalues: eig1 >= eig2 >= eig3
      eigvals(1) = eig1; eigvals(2) = eig2; eigvals(3) = eig3
      do i = 1, 2
        do j = i+1, 3
          if (eigvals(j) > eigvals(i)) then
            tmp = eigvals(i); eigvals(i) = eigvals(j); eigvals(j) = tmp
          end if
        end do
      end do
      eig1 = eigvals(1); eig2 = eigvals(2); eig3 = eigvals(3)

      Sp1 = eig1; Sp2 = eig2; Sp3 = eig3

      ! Eigenvectors: for each eigenvalue lambda_k, find null-space of (S-lambda_k*I)
      ! Use the cross-product method: take the cross product of any two rows.
      ! Build rows of (S - lambda*I):
      !   Row1: [sxx-lam,  sxy,     sxz    ]
      !   Row2: [sxy,      syy-lam, syz    ]
      !   Row3: [sxz,      syz,     szz-lam]

      call EigenVec3x3Sym(sxx,syy,szz,sxy,syz,sxz, eig1, xN1)
      call EigenVec3x3Sym(sxx,syy,szz,sxy,syz,sxz, eig2, xN2)
      call EigenVec3x3Sym(sxx,syy,szz,sxy,syz,sxz, eig3, xN3)

      ! Ensure right-handedness: xN3 = xN1 x xN2 if nearly degenerate
      xN3(1) = xN1(2)*xN2(3) - xN1(3)*xN2(2)
      xN3(2) = xN1(3)*xN2(1) - xN1(1)*xN2(3)
      xN3(3) = xN1(1)*xN2(2) - xN1(2)*xN2(1)
      vnorm = dsqrt(xN3(1)**2 + xN3(2)**2 + xN3(3)**2)
      if (vnorm > SMALL) then
        xN3 = xN3 / vnorm
      end if

      return
      end subroutine CalcPrincipalHB


! =============================================================================
!  EigenVec3x3Sym
!
!  Compute one eigenvector of a 3x3 real symmetric matrix S given
!  eigenvalue lam, using the cross-product of two rows of (S - lam*I).
!  The most linearly independent pair of rows is chosen automatically.
! =============================================================================
      subroutine EigenVec3x3Sym(sxx,syy,szz,sxy,syz,sxz, lam, xNout)

      implicit double precision (a-h, o-z)
      double precision, intent(in)  :: sxx,syy,szz,sxy,syz,sxz,lam
      double precision, intent(out) :: xNout(3)

      double precision :: r1(3), r2(3), r3(3)
      double precision :: c12(3), c13(3), c23(3)
      double precision :: n12, n13, n23, vnorm
      double precision, parameter :: SMALL = 1.0d-14

      ! Rows of (S - lam*I)
      r1 = (/ sxx - lam, sxy,       sxz       /)
      r2 = (/ sxy,       syy - lam, syz       /)
      r3 = (/ sxz,       syz,       szz - lam /)

      ! Cross products of each pair of rows
      c12(1) = r1(2)*r2(3) - r1(3)*r2(2)
      c12(2) = r1(3)*r2(1) - r1(1)*r2(3)
      c12(3) = r1(1)*r2(2) - r1(2)*r2(1)

      c13(1) = r1(2)*r3(3) - r1(3)*r3(2)
      c13(2) = r1(3)*r3(1) - r1(1)*r3(3)
      c13(3) = r1(1)*r3(2) - r1(2)*r3(1)

      c23(1) = r2(2)*r3(3) - r2(3)*r3(2)
      c23(2) = r2(3)*r3(1) - r2(1)*r3(3)
      c23(3) = r2(1)*r3(2) - r2(2)*r3(1)

      n12 = c12(1)**2 + c12(2)**2 + c12(3)**2
      n13 = c13(1)**2 + c13(2)**2 + c13(3)**2
      n23 = c23(1)**2 + c23(2)**2 + c23(3)**2

      ! Choose the pair with the largest cross-product norm
      if (n12 >= n13 .and. n12 >= n23) then
        xNout = c12;  vnorm = dsqrt(n12)
      else if (n13 >= n12 .and. n13 >= n23) then
        xNout = c13;  vnorm = dsqrt(n13)
      else
        xNout = c23;  vnorm = dsqrt(n23)
      end if

      if (vnorm > SMALL) then
        xNout = xNout / vnorm
      else
        ! Degenerate case: eigenvalue has multiplicity >= 2, pick a canonical dir
        xNout = (/ 1.0d0, 0.0d0, 0.0d0 /)
      end if

      return
      end subroutine EigenVec3x3Sym


! =============================================================================
!  Solve3x3HB
!
!  Solve the 3x3 linear system  A * x = b  using partial-pivot Gaussian
!  elimination.  Used by the Newton-Raphson loop in UMAT.
!
!  INFO = 0  on success
!  INFO = 1  if the matrix is (near-)singular
! =============================================================================
      subroutine Solve3x3HB(A, b, x, INFO)

      implicit double precision (a-h, o-z)
      double precision, intent(in)    :: A(3,3), b(3)
      double precision, intent(out)   :: x(3)
      integer,          intent(out)   :: INFO

      double precision :: Aug(3,4), factor, pivot, tmp
      integer :: i, j, k, ipiv
      double precision, parameter :: SMALL = 1.0d-30

      ! Build augmented matrix [A | b]
      do i = 1, 3
        do j = 1, 3
          Aug(i,j) = A(i,j)
        end do
        Aug(i,4) = b(i)
      end do

      INFO = 0

      ! Forward elimination with partial pivoting
      do k = 1, 2
        ! Find pivot row
        pivot = dabs(Aug(k,k)); ipiv = k
        do i = k+1, 3
          if (dabs(Aug(i,k)) > pivot) then
            pivot = dabs(Aug(i,k)); ipiv = i
          end if
        end do

        ! Swap rows k and ipiv
        if (ipiv /= k) then
          do j = k, 4
            tmp = Aug(k,j); Aug(k,j) = Aug(ipiv,j); Aug(ipiv,j) = tmp
          end do
        end if

        if (dabs(Aug(k,k)) < SMALL) then
          INFO = 1; return
        end if

        ! Eliminate column k below diagonal
        do i = k+1, 3
          factor = Aug(i,k) / Aug(k,k)
          do j = k, 4
            Aug(i,j) = Aug(i,j) - factor * Aug(k,j)
          end do
        end do
      end do

      ! Check last pivot
      if (dabs(Aug(3,3)) < SMALL) then
        INFO = 1; return
      end if

      ! Back substitution
      x(3) = Aug(3,4) / Aug(3,3)
      x(2) = (Aug(2,4) - Aug(2,3)*x(3)) / Aug(2,2)
      x(1) = (Aug(1,4) - Aug(1,2)*x(2) - Aug(1,3)*x(3)) / Aug(1,1)

      return
      end subroutine Solve3x3HB


! =============================================================================
!  MatTranspose  –  utility (unchanged from Anura3D original)
! =============================================================================
      Subroutine MatTranspose(A, Ia, AT, IAt, N1, N2)
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


    end module ModMohrCoulomb
