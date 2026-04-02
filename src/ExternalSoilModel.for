!**********************************************************************
!  ExternalSoilModel.for
!
!  Anura3D external soil model module — modified for Hoek-Brown 借壳法.
!
!  Compatible with Anura3D v2024 / v2025.
!  Based on the Anura3D MPM framework (LGPL v3).
!
!  Modification from stock Anura3D ExternalSoilModel.for:
!    In the Mohr-Coulomb branch of StressSolid, props(3)-(6) are
!    assigned as RAW Hoek-Brown parameters instead of the standard
!    transformed Mohr-Coulomb parameters.
!
!  ORIGINAL MC props assignment (commented out for reference):
!    props(3) = SIN(MatParams(IDSet)%FrictionAngle*(Pi/180.0))
!    props(4) = Particles(IDpt)%CohesionCosPhi
!    props(5) = SIN(MatParams(IDSet)%DilatancyAngle*(Pi/180.0))
!    props(6) = MatParams(IDSet)%TensileStrength
!
!  MODIFIED HB props assignment (借壳法):
!    props(3) = MatParams(IDSet)%FrictionAngle     -- sigma_ci (NO SIN)
!    props(4) = MatParams(IDSet)%DilatancyAngle    -- GSI      (NO SIN)
!    props(5) = MatParams(IDSet)%TensileStrength   -- mi       (direct)
!    props(6) = MatParams(IDSet)%ESM_Solid(6)      -- D        (disturbance)
!
!  Integration:
!    1. Replace src/ExternalSoilModel.for in your Anura3D source tree
!    2. Replace src/Soilmodels/A3DMohrCoulombStandard.f
!    3. Add src/Soilmodels/UMAT_MohrCoulombStandard.f to the project
!    4. Rebuild the Anura3D solution
!**********************************************************************

module ModExternalSoilModel

use ModMPMData
use ModGlobalConstants
use ModReadCalculationData
use ModReadMaterialData
use ModMPMInit
use user32
use kernel32
use ModMeshInfo
use ModLinearElasticity
use ModMohrCoulomb
use ModBingham

contains


subroutine StressSolid(IDpt, IDel, BMatrix, IEntityID)
!**********************************************************************
!
!  Function: calculate stresses at material point using external soil
!            models. Modified for Hoek-Brown 借壳法 in the MC branch.
!
!**********************************************************************

implicit none

    integer(INTEGER_TYPE), intent(in) :: IDpt
    integer(INTEGER_TYPE), intent(in) :: IDel
    real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(in) :: BMatrix
    integer(INTEGER_TYPE), intent(in) :: IEntityID

    ! local variables
    character(len=80) :: cmname
    integer(INTEGER_TYPE) :: I
    integer(INTEGER_TYPE) :: IDset
    integer(INTEGER_TYPE) :: ntens
    integer(INTEGER_TYPE), parameter :: nAddVar = 12
    real(REAL_TYPE), dimension(NPROPERTIES) :: props
    real(REAL_TYPE), dimension(nAddVar) :: AdditionalVar
    real(REAL_TYPE), dimension(MatParams(MaterialIDArray(IDpt))%UMATDimension) :: Stress, StrainIncr
    real(REAL_TYPE), dimension(NTENSOR) :: Sig0, StressIncr, StressPrinc, TempStrainIncr, TempStrainIncrPrevious
    real(REAL_TYPE), dimension(NSTATEVAR) :: StateVar
    real(REAL_TYPE) :: Eunloading, PlasticMultiplier
    character(len=64) :: NameModel
    logical :: IsUndrEffectiveStress
    real(REAL_TYPE) :: DSigWP
    real(REAL_TYPE) :: DSigGP
    real(REAL_TYPE) :: Bulk
    real(REAL_TYPE) :: DEpsVol
    procedure(DUMMYESM), pointer :: ESM

    ! get constitutive model in integration/material point
    IDset = MaterialIDArray(IDpt)
    NameModel = MatParams(IDset)%MaterialModel
    ntens = MatParams(IDset)%UMATDimension

    ! get strain increments in integration/material point
    TempStrainIncr = GetEpsStep(Particles(IDpt))

    StrainIncr = 0.0
    do I = 1, NTENSOR
        StrainIncr(I) = StrainIncr(I) + TempStrainIncr(I)
    end do

    DEpsVol = StrainIncr(1) + StrainIncr(2) + StrainIncr(3)

    IsUndrEffectiveStress = &
        ((CalParams%ApplyEffectiveStressAnalysis .and. &
          (trim(MatParams(IDSet)%MaterialType) == '2-phase')) .or. &
         (trim(MatParams(IDSet)%MaterialType) == SATURATED_SOIL_UNDRAINED_EFFECTIVE))

    ! initialise water pressure (only needed for undrained analyses)
    DSigWP = 0.0d0
    DSigGP = 0.0d0
    if (IsUndrEffectiveStress) then
        if (Particles(IDpt)%Porosity > 0.0) then
            Bulk = Particles(IDpt)%BulkWater / Particles(IDpt)%Porosity
            DSigWP = Bulk * DEpsVol
        else
            DSigWP = 0.0
        end if
        call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP)
    end if

    ! get stresses in integration/material point
    do I = 1, NTENSOR
        Sig0(I) = SigmaEff0Array(IDpt, I)
    end do
    Stress = 0.0
    do I = 1, NTENSOR
        Stress(I) = Stress(I) + Sig0(I)
    end do

    ! initialise state variables
    StateVar = ESMstatevArray(IDpt, :)

    ! Undrained effective stress pore pressure
    if (IsUndrEffectiveStress) then
        Particles(IDPt)%WaterPressure = Particles(IDPt)%WaterPressure + DSigWP
    end if

    ! get values of variables of interest for UMAT model
    AdditionalVar(1)  = Particles(IDPt)%Porosity
    AdditionalVar(2)  = Particles(IDPt)%WaterPressure
    AdditionalVar(3)  = Particles(IDPt)%WaterPressure0
    AdditionalVar(4)  = Particles(IDPt)%GasPressure
    AdditionalVar(5)  = Particles(IDPt)%GasPressure0
    AdditionalVar(6)  = Particles(IDPt)%DegreeSaturation
    AdditionalVar(7)  = CalParams%TotalRealTime
    AdditionalVar(8)  = CalParams%OverallRealTime
    AdditionalVar(9)  = CalParams%TimeIncrement
    AdditionalVar(10) = CalParams%IStep
    AdditionalVar(11) = CalParams%TimeStep
    AdditionalVar(12) = Particles(IDpt)%BulkWater

    ! get name of DLL
    cmname = MatParams(IDSet)%SoilModelDLL
    ! get material properties
    props = MatParams(IDSet)%ESM_Solid

    if (trim(NameModel)//char(0) == trim('linear_elasticity')//char(0)) then
        props(1) = Particles(IDpt)%ShearModulus
        cmname = UMAT_LINEAR_ELASTICITY
    elseif (trim(NameModel)//char(0) == trim(ESM_MOHR_COULOMB_STANDARD)//char(0)) then
        !---------------------------------------------------------------
        ! Hoek-Brown 借壳法: pass raw HB parameters through MC channel
        !
        ! ORIGINAL Mohr-Coulomb assignment (commented out):
        !   props(3) = SIN(MatParams(IDSet)%FrictionAngle*(Pi/180.0))
        !   props(4) = Particles(IDpt)%CohesionCosPhi
        !   props(5) = SIN(MatParams(IDSet)%DilatancyAngle*(Pi/180.0))
        !   props(6) = MatParams(IDSet)%TensileStrength
        !
        ! Modified for Hoek-Brown (NO SIN transforms):
        !   props(3) = sigma_ci  (stored in FrictionAngle field)
        !   props(4) = GSI       (stored in DilatancyAngle field)
        !   props(5) = mi        (stored in TensileStrength field)
        !   props(6) = D         (stored in ESM_Solid(6))
        !---------------------------------------------------------------
        props(1) = Particles(IDpt)%ShearModulus
        props(2) = MatParams(IDSet)%PoissonRatio
        props(3) = MatParams(IDSet)%FrictionAngle          ! sigma_ci (NO SIN)
        props(4) = MatParams(IDSet)%DilatancyAngle         ! GSI      (NO SIN)
        props(5) = MatParams(IDSet)%TensileStrength        ! mi       (direct)
        props(6) = MatParams(IDSet)%ESM_Solid(6)           ! D        disturbance
        cmname = UMAT_MOHR_COULOMB_STANDARD
    end if

    ! initialise ESM
    if ((CalParams%CPSversion == Anura3D_v2024) .or. &
        (CalParams%CPSversion == Anura3D_v2025)) then
        ESM => MatParams(IDSet)%ESM_POINTER
    end if

    call ESM(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, &
             StrainIncr, NSTATEVAR, StateVar, nAddVar, AdditionalVar, &
             cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)

    ! save unloading stiffness in Particles array
    Particles(IDpt)%ESM_UnloadingStiffness = Eunloading

    if (IsUndrEffectiveStress) then
        Particles(IDpt)%BulkWater = AdditionalVar(12)
    end if

    call SetIPL(IDpt, IDel, int(PlasticMultiplier))

    ! to use objective stress definition
    if (CalParams%ApplyObjectiveStress) then
        call Hill(IdEl, ELEMENTNODES, &
                  IncrementalDisplacementSoil(1:Counters%N, IEntityID), &
                  ReducedDof, ElementConnectivities, BMatrix, &
                  Sig0(1:NTENSOR), Stress(1:NTENSOR), DEpsVol)
    end if

    ! write new stresses to global array
    do I = 1, NTENSOR
        StressIncr(I) = Stress(I) - Sig0(I)
    end do

    ! save updated state variables
    ESMstatevArray(IDpt, :) = StateVar

    call CalculatePrincipalStresses(IDpt, Stress(1:NTENSOR), StressPrinc)
    call AssignStressStrainToGlobalArrayESM(IDpt, NTENSOR, StressIncr, &
                                           StressPrinc, StrainIncr)

    if (CalParams%ApplyBulkViscosityDamping) then
        RateVolStrain(IDEl) = DEpsVol / CalParams%TimeIncrement
        call CalculateViscousDamping_interface(IDpt, IDEl)
    end if

end subroutine StressSolid


subroutine CalculateViscousDamping_interface(ParticleID, IEl)
!**********************************************************************
!
!  Function: compute a pressure term introducing bulk viscosity
!            damping to the equation of motion.
!
!**********************************************************************

implicit none

    integer(INTEGER_TYPE), intent(in) :: ParticleID, IEl

    real(REAL_TYPE) :: ViscousDampingPressure = 0.0
    real(REAL_TYPE) :: Density = 0.0
    real(REAL_TYPE) :: ElementLMinLocal = 0.0
    real(REAL_TYPE) :: RateVolStrainLocal = 0.0
    real(REAL_TYPE) :: MaterialIndex = 0.0
    real(REAL_TYPE) :: DilationalWaveSpeed = 0.0
    logical :: IsUndrEffectiveStress

    if (.not. CalParams%ApplyBulkViscosityDamping) return

    MaterialIndex = MaterialIDArray(ParticleID)

    IsUndrEffectiveStress = &
        ((CalParams%ApplyEffectiveStressAnalysis .and. &
          (trim(MatParams(MaterialIndex)%MaterialType) == '2-phase')) .or. &
         (trim(MatParams(MaterialIndex)%MaterialType) == SATURATED_SOIL_UNDRAINED_EFFECTIVE))

    if (IsUndrEffectiveStress &
        .or. ((CalParams%NumberOfPhases == 2) .or. (CalParams%NumberOfPhases == 3))) then
        Density = MatParams(MaterialIndex)%DensityMixture / 1000.0
    else
        Density = (1 - MatParams(MaterialIndex)%InitialPorosity) * &
                  MatParams(MaterialIndex)%DensitySolid / 1000.0
    end if

    ElementLMinLocal = ElementLMin(IEl)
    RateVolStrainLocal = RateVolStrain(IEl)

    call GetWaveSpeed(ParticleID, DilationalWaveSpeed)

    ViscousDampingPressure = CalParams%BulkViscosityDamping1 * &
        Density * DilationalWaveSpeed * ElementLMinLocal * RateVolStrainLocal

    if ((RateVolStrainLocal < 0.0) .and. (CalParams%BulkViscosityDamping2 > 0.0)) then
        ViscousDampingPressure = ViscousDampingPressure + &
            Density * (CalParams%BulkViscosityDamping2 * ElementLMinLocal * RateVolStrainLocal)**2
    end if

    Particles(ParticleID)%DBulkViscousPressure = ViscousDampingPressure

end subroutine CalculateViscousDamping_interface


end module ModExternalSoilModel
