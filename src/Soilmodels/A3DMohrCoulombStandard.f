!**********************************************************************
!  A3DMohrCoulombStandard.f
!
!  Module ModMohrCoulomb for Anura3D ESM framework.
!  Contains:
!    ESM_MC  — entry point called by the ESM procedure pointer
!    UMAT    — Hoek-Brown constitutive model (借壳 / shell borrowing)
!    CarSig  — principal-to-Cartesian stress back-transformation
!    PrnSig  — declared EXTERNAL (provided by Anura3D GetPrinStress.FOR)
!
!  PROPS mapping (passed from ExternalSoilModel.for via MC channel):
!    PROPS(1) = G        — shear modulus
!    PROPS(2) = nu       — Poisson's ratio
!    PROPS(3) = sigma_ci — uniaxial compressive strength of intact rock
!    PROPS(4) = GSI      — Geological Strength Index
!    PROPS(5) = mi       — intact rock parameter
!    PROPS(6) = D        — disturbance factor (0–1)
!
!  STATEV mapping:
!    STATEV(1) = accumulated equivalent plastic strain
!    STATEV(2) = plastic flag (0=elastic, 1=shear, 2=tensile)
!    STATEV(3) = mb  (derived, for post-processing)
!    STATEV(4) = s   (derived, for post-processing)
!    STATEV(5) = a   (derived, for post-processing)
!
!  Sign convention: compression positive (consistent with Anura3D).
!**********************************************************************

    module ModMohrCoulomb
    contains

    Subroutine ESM_MC(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER,&
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
        double precision :: Porosity, WaterPressure, WaterPressure0, GasPressure, GasPressure0, DegreeSaturation

        integer :: ndi, nshr, layer, kspt, kstep, kinc

        allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
              coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )

!Initialization
        Eunloading = 0.0
        PlasticMultiplier = 0.0

!Rename additional variables
        Porosity = AdditionalVar(1)
        WaterPressure = AdditionalVar(2)
        WaterPressure0 = AdditionalVar(3)
        GasPressure = AdditionalVar(4)
        GasPressure0 = AdditionalVar(5)
        DegreeSaturation = AdditionalVar(6)
        time(1) = AdditionalVar(7)
        time(2) = AdditionalVar(8)
        dtime = AdditionalVar(9)
        IStep = AdditionalVar(10)
        TimeStep = AdditionalVar(11)

!Call the UMAT (module procedure — no EXTERNAL needed)
        call UMAT(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
         dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, &
         dfgrd1, noel, npt, layer, kspt, kstep, kinc)


!---Definition of Eunloading -> required to define the max time step
      Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))

! PlasticMultiplier: output for plotting plastic points
      PlasticMultiplier = STATEV(2)

        return

    end subroutine ESM_MC

!**********************************************************************
!  UMAT: Hoek-Brown constitutive model
!  Called from ESM_MC via standard UMAT interface.
!  Uses PrnSig from Anura3D framework (src/GetPrinStress.FOR).
!**********************************************************************
    SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
      RPL,DDSDDT,DRPLDE,DRPLDT, &
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      implicit double precision (a-h, o-z)
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV), &
        DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS), &
        STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1), &
        PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      EXTERNAL PrnSig

!---Local variables
      integer :: i, j, iter
      integer, parameter :: MAXITER = 30
      double precision, parameter :: FTOL = 1.0d-10

      double precision :: G, nu, sigci, GSI, xmi, D
      double precision :: mb, s, a
      double precision :: K, lam, twoG
      double precision :: SigE(6)
      double precision :: xN1(3), xN2(3), xN3(3)
      double precision :: Sig1, Sig2, Sig3, P, Q
      double precision :: hb_term, F, dFdSig3, dFdSig1
      double precision :: dLambda, Sig1_new, Sig3_new
      double precision :: sig_t, Sig1_tc, Sig3_tc
      double precision :: epsp_new
      integer :: plastic_flag
      double precision :: SigC(6)
      double precision :: Sig1_iter, Sig3_iter, F_iter, dF_iter
      double precision :: dSig1, dSig3
      logical :: tension_cut

!---Read material properties
      G     = PROPS(1)
      nu    = PROPS(2)
      sigci = PROPS(3)
      GSI   = PROPS(4)
      xmi   = PROPS(5)
      D     = PROPS(6)

!---Compute Hoek-Brown derived parameters (2002 generalized HB criterion)
!   Reference: Hoek, Carranza-Torres & Corkum (2002), Proc. NARMS-TAC, 267-273.
      mb = xmi * exp((GSI - 100.0d0) / (28.0d0 - 14.0d0*D))
      s  = exp((GSI - 100.0d0) / (9.0d0 - 3.0d0*D))
      a  = 0.5d0 + (1.0d0/6.0d0)*(exp(-GSI/15.0d0) - exp(-20.0d0/3.0d0))

!---Store derived parameters for post-processing
      STATEV(3) = mb
      STATEV(4) = s
      STATEV(5) = a

!---Elastic stiffness constants
      K    = G * 2.0d0*(1.0d0 + nu) / (3.0d0*(1.0d0 - 2.0d0*nu))
      lam  = K - 2.0d0*G/3.0d0
      twoG = 2.0d0*G

!---Build elastic stiffness matrix DDSDDE
      do i = 1, NTENS
        do j = 1, NTENS
          DDSDDE(i,j) = 0.0d0
        end do
      end do
      do i = 1, 3
        do j = 1, 3
          DDSDDE(i,j) = lam
        end do
        DDSDDE(i,i) = lam + twoG
      end do
      do i = 4, NTENS
        DDSDDE(i,i) = G
      end do

!---Elastic predictor: trial stress = current stress + D:dStrain
      do i = 1, NTENS
        SigE(i) = STRESS(i)
      end do
      do i = 1, 3
        do j = 1, 3
          SigE(i) = SigE(i) + DDSDDE(i,j)*DSTRAN(j)
        end do
      end do
      do i = 4, NTENS
        SigE(i) = SigE(i) + DDSDDE(i,i)*DSTRAN(i)
      end do

!---Compute principal stresses of trial stress
!   PrnSig(IOpt, S, xN1, xN2, xN3, S1, S2, S3, P, Q)
!   from Anura3D framework (src/GetPrinStress.FOR)
!   Returns Sig1 <= Sig2 <= Sig3 (Sig3 = most compressive)
      call PrnSig(1, SigE, xN1, xN2, xN3, Sig1, Sig2, Sig3, P, Q)

!---Evaluate Hoek-Brown yield function
!   F = Sig3 - Sig1 - sigci*(mb*Sig1/sigci + s)^a
      plastic_flag = 0
      tension_cut  = .false.

      hb_term = mb*Sig1/sigci + s

!---Tension cutoff: HB is not defined when hb_term < 0
      if (hb_term < 0.0d0) then
        tension_cut = .true.
        sig_t = -s*sigci/mb
        Sig1_new = sig_t
        Sig3_new = Sig3
        if (Sig3_new < Sig1_new) Sig3_new = Sig1_new
        plastic_flag = 2

        ! Re-check shear yield with corrected sig1
        hb_term = mb*Sig1_new/sigci + s
        if (hb_term >= 0.0d0) then
          F = Sig3_new - Sig1_new - sigci*(hb_term**a)
          if (F > FTOL) then
            Sig1_iter = Sig1_new
            Sig3_iter = Sig3_new
            do iter = 1, MAXITER
              hb_term  = mb*Sig1_iter/sigci + s
              if (hb_term <= 0.0d0) exit
              F_iter   = Sig3_iter - Sig1_iter - sigci*(hb_term**a)
              if (abs(F_iter) <= FTOL) exit
              dFdSig3  = 1.0d0
              dFdSig1  = -1.0d0 - a*mb*(hb_term**(a - 1.0d0))
              dF_iter  = dFdSig1*(-twoG*dFdSig1) + dFdSig3*(twoG*dFdSig3)
              if (abs(dF_iter) < 1.0d-20) exit
              dLambda  = F_iter / dF_iter
              Sig1_iter = Sig1_iter - dLambda*(-twoG*dFdSig1)
              Sig3_iter = Sig3_iter - dLambda*(twoG*dFdSig3)
            end do
            Sig1_new = Sig1_iter
            Sig3_new = Sig3_iter
            plastic_flag = 1
          end if
        end if

        goto 700

      end if

!---Evaluate yield function
      F = Sig3 - Sig1 - sigci*(hb_term**a)

      if (F <= FTOL) then
        ! Elastic step
        Sig1_new = Sig1
        Sig3_new = Sig3
        plastic_flag = 0
        goto 700
      end if

!---Plastic return mapping (Newton-Raphson)
      plastic_flag = 1
      Sig1_iter = Sig1
      Sig3_iter = Sig3

      do iter = 1, MAXITER
        hb_term = mb*Sig1_iter/sigci + s
        if (hb_term < 0.0d0) then
          Sig1_iter = -s*sigci/mb
          plastic_flag = 2
          exit
        end if
        F_iter  = Sig3_iter - Sig1_iter - sigci*(hb_term**a)
        if (abs(F_iter) <= FTOL) exit
        dFdSig3 = 1.0d0
        dFdSig1 = -1.0d0 - a*mb*(hb_term**(a - 1.0d0))
        dF_iter = dFdSig1*(-twoG*dFdSig1) + dFdSig3*(twoG*dFdSig3)
        if (abs(dF_iter) < 1.0d-20) exit
        dLambda   = F_iter / dF_iter
        Sig1_iter = Sig1_iter - dLambda*(-twoG*dFdSig1)
        Sig3_iter = Sig3_iter - dLambda*(twoG*dFdSig3)
      end do

      Sig1_new = Sig1_iter
      Sig3_new = Sig3_iter

!---Back-transformation to Cartesian stress tensor
  700 continue

      epsp_new = STATEV(1) + sqrt((abs(Sig1_new - Sig1)**2 + &
                 abs(Sig3_new - Sig3)**2) / 2.0d0) / max(twoG, 1.0d-20)
      STATEV(1) = epsp_new
      STATEV(2) = dble(plastic_flag)

      ! Sig2 unchanged (no plastic correction on intermediate principal stress)
      call CarSig(Sig1_new, Sig2, Sig3_new, xN1, xN2, xN3, NTENS, SigC)

      do i = 1, NTENS
        STRESS(i) = SigC(i)
      end do

      return
    END SUBROUTINE UMAT

!**********************************************************************
!  CarSig: Transform principal stresses back to Cartesian (Voigt) form
!**********************************************************************
    Subroutine CarSig(Sig1, Sig2, Sig3, xN1, xN2, xN3, ntens, SigC)
      implicit double precision (a-h, o-z)
      integer, intent(in) :: ntens
      double precision, intent(in)  :: Sig1, Sig2, Sig3
      double precision, intent(in)  :: xN1(3), xN2(3), xN3(3)
      double precision, intent(out) :: SigC(ntens)
      double precision :: xP1(3,3), xP2(3,3), xP3(3,3), SigCart(3,3)
      integer :: i, j

      do i = 1, 3
        do j = 1, 3
          xP1(i,j) = xN1(i) * xN1(j)
          xP2(i,j) = xN2(i) * xN2(j)
          xP3(i,j) = xN3(i) * xN3(j)
        end do
      end do

      do i = 1, 3
        do j = 1, 3
          SigCart(i,j) = Sig1*xP1(i,j) + Sig2*xP2(i,j) + Sig3*xP3(i,j)
        end do
      end do

      SigC(1) = SigCart(1,1)
      SigC(2) = SigCart(2,2)
      SigC(3) = SigCart(3,3)
      if (ntens >= 4) SigC(4) = SigCart(1,2)
      if (ntens >= 5) SigC(5) = SigCart(2,3)
      if (ntens >= 6) SigC(6) = SigCart(1,3)

      return
    end subroutine CarSig

    end module ModMohrCoulomb
