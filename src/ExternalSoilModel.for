!**********************************************************************
!  ExternalSoilModel.for
!  Anura3D external soil model interface — modified for Hoek-Brown
!  "借壳法": routes Hoek-Brown parameters through the MC channel.
!
!  Key change from original Anura3D:
!    props(3)-props(6) in the MC branch are assigned WITHOUT SIN() transforms,
!    passing raw material parameters directly to the UMAT (Hoek-Brown) routine.
!
!  PROPS mapping for Hoek-Brown via MC channel:
!    PROPS(1) = G        (from Particles%ShearModulus)
!    PROPS(2) = nu       (from MatParams%PoissonRatio)
!    PROPS(3) = sigma_ci (from MatParams%FrictionAngle  -- NO SIN transform)
!    PROPS(4) = GSI      (from MatParams%DilatancyAngle -- NO SIN transform)
!    PROPS(5) = mi       (from MatParams%TensileStrength -- direct)
!    PROPS(6) = D        (from MatParams%ESM_Solid(6))
!**********************************************************************

    subroutine ESM(IDpt, IDset, IDTask, IsUndrainedCalculation, &
                   Porosity, WaterPressure, WaterPressure0, &
                   GasPressure, GasPressure0, DegreeSaturation, &
                   dt, h, &
                   X, Y, Z, &
                   iAbort)

    use ModGlobalConstants
    use ModMohrCoulomb
    use ModLinearElasticity
    use ModBingham
    use ModMesh
    use ModMPMData

    implicit double precision (a-h, o-z)

    integer, intent(in)    :: IDpt, IDset, IDTask
    integer, intent(in)    :: IsUndrainedCalculation
    double precision, intent(in) :: Porosity, WaterPressure, WaterPressure0
    double precision, intent(in) :: GasPressure, GasPressure0, DegreeSaturation
    double precision, intent(in) :: dt, h
    double precision, intent(in) :: X, Y, Z
    integer, intent(out)   :: iAbort

    integer :: NOEL, NPT, IDSET_ESM, NSTATEV, NADDVAR, NPROPS, NUMBEROFPHASES, NTENS
    double precision :: EUNLOADING, PLASTICMULTIPLIER
    CHARACTER*80 CMNAME
    double precision, dimension(:), allocatable :: STRESS, DSTRAN, STATEV, ADDITIONALVAR, PROPS

    integer :: i, MaterialID
    double precision :: pi

    pi = acos(-1.0d0)
    iAbort = 0

    ! Map Anura3D particle/material data to local variables
    NOEL          = IDpt
    NPT           = IDpt
    IDSET_ESM     = IDset
    NUMBEROFPHASES = 1
    NTENS         = NTENSOR   ! from ModGlobalConstants (typically 6 for 3D)
    MaterialID    = IDset

    ! Determine NSTATEV and NPROPS based on material model
    select case (MatParams(IDset)%MaterialModel)

    !------------------------------------------------------------------
    ! Mohr-Coulomb branch — used here for Hoek-Brown (借壳法)
    !------------------------------------------------------------------
    case (ESM_MOHRCOULOMB)

      NPROPS  = 8
      NSTATEV = 10
      NADDVAR = 11
      CMNAME  = 'MOHR-COULOMB'

      allocate( STRESS(NTENS), DSTRAN(NTENS), STATEV(NSTATEV), &
                ADDITIONALVAR(NADDVAR), PROPS(NPROPS) )

      ! Copy current stress from particle
      do i = 1, NTENS
        STRESS(i) = Particles(IDpt)%Stress(i)
      end do

      ! Copy strain increment
      do i = 1, NTENS
        DSTRAN(i) = Particles(IDpt)%StrainIncrement(i)
      end do

      ! Copy state variables
      do i = 1, NSTATEV
        STATEV(i) = Particles(IDpt)%StateVariables(i)
      end do

      ! Additional variables (pore pressures, time, etc.)
      ADDITIONALVAR(1)  = Porosity
      ADDITIONALVAR(2)  = WaterPressure
      ADDITIONALVAR(3)  = WaterPressure0
      ADDITIONALVAR(4)  = GasPressure
      ADDITIONALVAR(5)  = GasPressure0
      ADDITIONALVAR(6)  = DegreeSaturation
      ADDITIONALVAR(7)  = TimeInformation%CurrentTime
      ADDITIONALVAR(8)  = TimeInformation%CurrentTime
      ADDITIONALVAR(9)  = dt
      ADDITIONALVAR(10) = TimeInformation%IStep
      ADDITIONALVAR(11) = TimeInformation%TimeStep

      !----------------------------------------------------------------
      ! PROPS assignment for Hoek-Brown (NO SIN transforms):
      !   PROPS(1) = G        -- shear modulus
      !   PROPS(2) = nu       -- Poisson's ratio
      !   PROPS(3) = sigma_ci -- from FrictionAngle field (raw, no SIN)
      !   PROPS(4) = GSI      -- from DilatancyAngle field (raw, no SIN)
      !   PROPS(5) = mi       -- from TensileStrength field (direct)
      !   PROPS(6) = D        -- disturbance factor from ESM_Solid(6)
      !----------------------------------------------------------------
      props(1) = Particles(IDpt)%ShearModulus
      props(2) = MatParams(IDSet)%PoissonRatio
      props(3) = MatParams(IDSet)%FrictionAngle          ! sigma_ci (NO SIN transform)
      props(4) = MatParams(IDSet)%DilatancyAngle          ! GSI      (NO SIN transform)
      props(5) = MatParams(IDSet)%TensileStrength          ! mi       (direct)
      props(6) = MatParams(IDSet)%ESM_Solid(6)             ! D        disturbance factor
      props(7) = 0.0d0
      props(8) = 0.0d0

      ! Call ESM_MC (from module ModMohrCoulomb) — dispatches to UMAT (Hoek-Brown)
      call ESM_MC(NPT, NOEL, IDSET_ESM, STRESS, EUNLOADING, PLASTICMULTIPLIER, &
                  DSTRAN, NSTATEV, STATEV, NADDVAR, ADDITIONALVAR, CMNAME, &
                  NPROPS, PROPS, NUMBEROFPHASES, NTENS)

      ! Write back updated stress and state variables to particle
      do i = 1, NTENS
        Particles(IDpt)%Stress(i) = STRESS(i)
      end do
      do i = 1, NSTATEV
        Particles(IDpt)%StateVariables(i) = STATEV(i)
      end do

      Particles(IDpt)%Eunloading      = EUNLOADING
      Particles(IDpt)%PlasticMultiplier = PLASTICMULTIPLIER

    !------------------------------------------------------------------
    ! Linear Elasticity branch (unchanged)
    !------------------------------------------------------------------
    case (ESM_LINEAR)

      NPROPS  = 2
      NSTATEV = 0
      NADDVAR = 11
      CMNAME  = 'LINEAR-ELASTIC'

      allocate( STRESS(NTENS), DSTRAN(NTENS), STATEV(1), &
                ADDITIONALVAR(NADDVAR), PROPS(NPROPS) )

      do i = 1, NTENS
        STRESS(i) = Particles(IDpt)%Stress(i)
        DSTRAN(i) = Particles(IDpt)%StrainIncrement(i)
      end do

      ADDITIONALVAR(1)  = Porosity
      ADDITIONALVAR(2)  = WaterPressure
      ADDITIONALVAR(3)  = WaterPressure0
      ADDITIONALVAR(4)  = GasPressure
      ADDITIONALVAR(5)  = GasPressure0
      ADDITIONALVAR(6)  = DegreeSaturation
      ADDITIONALVAR(7)  = TimeInformation%CurrentTime
      ADDITIONALVAR(8)  = TimeInformation%CurrentTime
      ADDITIONALVAR(9)  = dt
      ADDITIONALVAR(10) = TimeInformation%IStep
      ADDITIONALVAR(11) = TimeInformation%TimeStep

      props(1) = Particles(IDpt)%ShearModulus
      props(2) = MatParams(IDSet)%PoissonRatio

      call ESM_LINEAR(NPT, NOEL, IDSET_ESM, STRESS, EUNLOADING, PLASTICMULTIPLIER, &
                      DSTRAN, NSTATEV, STATEV, NADDVAR, ADDITIONALVAR, CMNAME, &
                      NPROPS, PROPS, NUMBEROFPHASES, NTENS)

      do i = 1, NTENS
        Particles(IDpt)%Stress(i) = STRESS(i)
      end do

      Particles(IDpt)%Eunloading = EUNLOADING

    !------------------------------------------------------------------
    ! Bingham branch (unchanged)
    !------------------------------------------------------------------
    case (ESM_BINGHAM)

      NPROPS  = 4
      NSTATEV = 0
      NADDVAR = 11
      CMNAME  = 'BINGHAM'

      allocate( STRESS(NTENS), DSTRAN(NTENS), STATEV(1), &
                ADDITIONALVAR(NADDVAR), PROPS(NPROPS) )

      do i = 1, NTENS
        STRESS(i) = Particles(IDpt)%Stress(i)
        DSTRAN(i) = Particles(IDpt)%StrainIncrement(i)
      end do

      ADDITIONALVAR(1)  = Porosity
      ADDITIONALVAR(2)  = WaterPressure
      ADDITIONALVAR(3)  = WaterPressure0
      ADDITIONALVAR(4)  = GasPressure
      ADDITIONALVAR(5)  = GasPressure0
      ADDITIONALVAR(6)  = DegreeSaturation
      ADDITIONALVAR(7)  = TimeInformation%CurrentTime
      ADDITIONALVAR(8)  = TimeInformation%CurrentTime
      ADDITIONALVAR(9)  = dt
      ADDITIONALVAR(10) = TimeInformation%IStep
      ADDITIONALVAR(11) = TimeInformation%TimeStep

      props(1) = Particles(IDpt)%ShearModulus
      props(2) = MatParams(IDSet)%PoissonRatio
      props(3) = MatParams(IDSet)%YieldStress
      props(4) = MatParams(IDSet)%ViscosityCoefficient

      call ESM_BINGHAM(NPT, NOEL, IDSET_ESM, STRESS, EUNLOADING, PLASTICMULTIPLIER, &
                       DSTRAN, NSTATEV, STATEV, NADDVAR, ADDITIONALVAR, CMNAME, &
                       NPROPS, PROPS, NUMBEROFPHASES, NTENS)

      do i = 1, NTENS
        Particles(IDpt)%Stress(i) = STRESS(i)
      end do

      Particles(IDpt)%Eunloading = EUNLOADING

    case default
      iAbort = 1

    end select

    ! Cleanup
    if (allocated(STRESS))        deallocate(STRESS)
    if (allocated(DSTRAN))        deallocate(DSTRAN)
    if (allocated(STATEV))        deallocate(STATEV)
    if (allocated(ADDITIONALVAR)) deallocate(ADDITIONALVAR)
    if (allocated(PROPS))         deallocate(PROPS)

    return
    end subroutine ESM
