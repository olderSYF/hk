	!*****************************************************************************
    !                                       ____  _____  
    !           /\                         |___ \|  __ \ 
    !          /  \   _ __  _   _ _ __ __ _  __) | |  | |
    !         / /\ \ | '_ \| | | | '__/ _` ||__ <| |  | |
    !        / ____ \| | | | |_| | | | (_| |___) | |__| |
    !       /_/    \_\_| |_|\__,_|_|  \__,_|____/|_____/ 
    !
    !
	!	Anura3D - Numerical modelling and simulation of large deformations 
    !   and soil-water-structure interaction using the material point method (MPM)
    !
    !	Copyright (C) 2025  Members of the Anura3D MPM Research Community 
    !   (See Contributors file "Contributors.txt")
    !
    !	This program is free software: you can redistribute it and/or modify
    !	it under the terms of the GNU Lesser General Public License as published by
    !	the Free Software Foundation, either version 3 of the License, or
    !	(at your option) any later version.
    !
    !	This program is distributed in the hope that it will be useful,
    !	but WITHOUT ANY WARRANTY; without even the implied warranty of
    !	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !	GNU Lesser General Public License for more details.
    !
    !	You should have received a copy of the GNU Lesser General Public License
    !	along with this program.  If not, see <https://www.gnu.org/licenses/>.
	!
    !*****************************************************************************
    !
    ! ============================================================================
    ! MODIFICATION NOTE (Hoek-Brown "Shell" Approach):
    !
    ! The Mohr-Coulomb routing block below has been modified so that the raw GiD
    ! input values are passed directly into PROPS without any trigonometric
    ! transformation.  This lets the same Mohr-Coulomb data channel carry
    ! Hoek-Brown parameters to the UMAT in A3DMohrCoulombStandard.f:
    !
    !   GiD field           PROPS index   HB parameter
    !   -------------------------------------------------
    !   Young's Modulus     (computed G)  shear modulus G [kPa]
    !   Poisson Ratio       PROPS(2)      nu
    !   Friction Angle      PROPS(3)      sigma_ci [kPa]  <-- raw, no SIN()
    !   Cohesion            PROPS(4)      GSI             <-- raw, no CohesionCosPhi
    !   Dilatancy Angle     PROPS(5)      mi              <-- raw, no SIN()
    !   Tensile Strength    PROPS(6)      D (disturbance)
    !
    ! NOTE on PROPS(4)/Cohesion:
    !   If MatParams(IDSet)%Cohesion does not compile in your Anura3D version,
    !   replace it with:  MatParams(IDSet)%ESM_Solid(4)
    !   (the ESM_Solid array holds the raw GiD values as read from the .GOM file).
    ! ============================================================================


module ModExternalSoilModel
!**********************************************************************
!
!  Function: Contains the routines related to calling user-defined soil
!            models in external DLLs.
!     $Revision: 9064 $
!     $Date: 2021-02-20 10:27:53 +0100 (sab, 20 feb 2021) $
!
!**********************************************************************
      
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
!    Function:  calculate stresses at material point using external
!               soil models
!
!*********************************************************************        
        
implicit none
        
    integer(INTEGER_TYPE), intent(in) :: IDpt        ! global integration/material point number
    integer(INTEGER_TYPE), intent(in) :: IDel        ! global element number
    ! B-matrix at the considered integration point
    ! (here only used if ApplyObjectiveStress=TRUE)
    real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(in) :: BMatrix    
    integer(INTEGER_TYPE), intent(in) :: IEntityID   ! entity ID
    ! local variables
    character(len=80) :: cmname
    integer(INTEGER_TYPE) :: I                       ! counter
    integer(INTEGER_TYPE) :: IDset                   ! ID of material parameter set
    integer(INTEGER_TYPE) :: ntens                   ! Dimension of stress vector to pass to ESM 
    integer(INTEGER_TYPE), parameter :: nAddVar = 12
    real(REAL_TYPE), dimension(NPROPERTIES) :: props ! array of material properties
    real(REAL_TYPE), dimension(nAddVar) :: AdditionalVar
    real(REAL_TYPE), dimension(MatParams(MaterialIDArray(IDpt))%UMATDimension) :: Stress, StrainIncr
    real(REAL_TYPE), dimension(NTENSOR) :: Sig0, StressIncr, StressPrinc, TempStrainIncr, TempStrainIncrPrevious
    real(REAL_TYPE), dimension(NSTATEVAR) :: StateVar
    real(REAL_TYPE) :: Eunloading, PlasticMultiplier
    character(len=64) :: NameModel                   ! name of the constitutive model
    logical :: IsUndrEffectiveStress
    real(REAL_TYPE) :: DSigWP                        ! Change of water pressure at integration point 
    real(REAL_TYPE) :: DSigGP                        ! Change of gas pressure at integration point 
    real(REAL_TYPE) :: Bulk                          ! Bulk modulus
    real(REAL_TYPE) :: DEpsVol                       ! Incremental volumetric strain (water)
    procedure(DUMMYESM), pointer :: ESM
    
    ! get constitutive model in integration/material point
    IDset = MaterialIDArray(IDpt)       ! material number stored in $$MATERIAL_INDEX in the GOM-file
    NameModel = MatParams(IDset)%MaterialModel  ! name of constitutive model as specified in GOM-file
    ntens = MatParams(IDset)%UMATDimension      ! 2D or 3D formulation of the External soil model   
          
    ! get strain increments in integration/material point
    TempStrainIncr = GetEpsStep(Particles(IDpt))     ! incremental strain vector assigned to point
    
        
    StrainIncr = 0.0

    do I=1, NTENSOR
    StrainIncr(I) = StrainIncr(I) + TempStrainIncr(I)
    enddo 
        
    DEpsVol = StrainIncr(1) + StrainIncr(2) + StrainIncr(3) ! volumetric strain, valid for 2D and 3D
          
    IsUndrEffectiveStress = &
    !code version 2016 and previous
    ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(IDSet)%MaterialType)=='2-phase')) .or. &
    !code version 2017.1 and following
    (trim(MatParams(IDSet)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
          
    ! initialise water pressure (only needed for undrained analyses)
    DSigWP = 0.0d0
    DSigGP = 0.0d0
    ! for effective stress analysis
    if (IsUndrEffectiveStress) then
        if (Particles(IDpt)%Porosity > 0.0) then
        Bulk = Particles(IDpt)%BulkWater / Particles(IDpt)%Porosity ! kN/m2
        DSigWP = Bulk * DEpsVol
        else
        DSigWP = 0.0
        end if
        call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP)
    end if ! effective stress analysis
          
    ! get stresses in integration/material point      
    do I = 1, NTENSOR
    Sig0(I) = SigmaEff0Array(IDpt, I) ! get initial stress of current step assigned to point 
    end do
    Stress=0.0
    do I=1, NTENSOR
        Stress(I) = Stress(I) + Sig0(I) !To use 3D UMAT also for 2D formulations
    enddo 
          
    ! initialise state variables
    StateVar = ESMstatevArray(IDpt,:)
          
    ! Undrained effective stress pore pressure
    if (IsUndrEffectiveStress) then
    Particles(IDPt)%WaterPressure = Particles(IDPt)%WaterPressure + DSigWP
    end if 
    
    !get values of variables of interest for UMAT model
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
    ! get material properties (initialise with raw ESM_Solid values)
    props = MatParams(IDSet)%ESM_Solid

    ! -----------------------------------------------------------------------
    ! Model routing: assign PROPS for the selected constitutive model
    ! -----------------------------------------------------------------------
    if (trim(NameModel)//char(0) == trim('linear_elasticity')//char(0)) then
    props(1) = Particles(IDpt)%ShearModulus ! shear modulus, G
    cmname = UMAT_LINEAR_ELASTICITY

    elseif (trim(NameModel)//char(0) == trim(ESM_MOHR_COULOMB_STANDARD)//char(0)) then
    ! *** Hoek-Brown "shell" modification ***
    ! Raw GiD values are passed directly so the UMAT in
    ! A3DMohrCoulombStandard.f can interpret them as HB parameters.
    ! Mapping:  GiD field           -> PROPS index -> HB meaning
    !           Young's Modulus(G)  -> PROPS(1)    -> shear modulus G [kPa]
    !           Poisson Ratio       -> PROPS(2)    -> nu
    !           Friction Angle      -> PROPS(3)    -> sigma_ci [kPa]  (NO SIN transform)
    !           Cohesion            -> PROPS(4)    -> GSI              (NO CohesionCosPhi)
    !           Dilatancy Angle     -> PROPS(5)    -> mi               (NO SIN transform)
    !           Tensile Strength    -> PROPS(6)    -> D (disturbance factor)
    props(1) = Particles(IDpt)%ShearModulus            ! shear modulus, G
    props(2) = MatParams(IDSet)%PoissonRatio           ! Poisson ratio nu
    props(3) = MatParams(IDSet)%FrictionAngle          ! sigma_ci (GiD "Friction Angle" field, raw value in kPa)
    props(4) = MatParams(IDSet)%Cohesion               ! GSI      (GiD "Cohesion" field, raw value)
    !           MatParams%Cohesion is the standard field name in Anura3D 2025.
    !           If your build reports an undefined member, use instead:
    !           props(4) = MatParams(IDSet)%ESM_Solid(4)
    !           (ESM_Solid holds the raw GiD values as-read from the .GOM file.)
    props(5) = MatParams(IDSet)%DilatancyAngle         ! mi       (GiD "Dilatancy Angle" field, raw value)
    props(6) = MatParams(IDSet)%TensileStrength        ! D        (GiD "Tensile Strength" field)
    cmname = UMAT_MOHR_COULOMB_STANDARD

    endif          

    ! initialise UMAT pointer
    if ((CalParams%CPSversion==Anura3D_v2024).or.(CalParams%CPSversion==Anura3D_v2025)) then
    ESM => MatParams(IDSet)%ESM_POINTER
    else
    !ESM = GetProcAddress(MatParams(IDSet)%SoilModelDLLHandle) 
    end if

    call ESM(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, &
             NSTATEVAR, StateVar, nAddVar, AdditionalVar, cmname, NPROPERTIES, props, &
             CalParams%NumberOfPhases, ntens)

    ! save unloading stiffness in Particles array  
    Particles(IDpt)%ESM_UnloadingStiffness = Eunloading
                 
    if (IsUndrEffectiveStress) then
        Particles(IDpt)%BulkWater = AdditionalVar(12)
    end if
    
    call SetIPL(IDpt, IDel, int(PlasticMultiplier))

    ! to use objective stress definition
    if (CalParams%ApplyObjectiveStress) then ! Consider large deformation terms
    call Hill(IdEl, ELEMENTNODES, IncrementalDisplacementSoil(1:Counters%N, IEntityID),  &
                     ReducedDof, ElementConnectivities, BMatrix, Sig0(1:NTENSOR), Stress(1:NTENSOR), DEpsVol)
    end if ! objective stress            
            
    ! write new stresses to global array
    do I=1, NTENSOR
        StressIncr(I) = Stress(I) - Sig0(I)
    enddo             
                               
    ! save updated state variables in Particles array
    ESMstatevArray(IDpt,:) = StateVar
          
    call CalculatePrincipalStresses(IDpt, Stress(1:NTENSOR), StressPrinc)
    call AssignStressStrainToGlobalArrayESM(IDpt, NTENSOR, StressIncr, StressPrinc, StrainIncr)

    if (CalParams%ApplyBulkViscosityDamping) then
    RateVolStrain(IDEl) = DEpsVol / CalParams%TimeIncrement
    call CalculateViscousDamping_interface(IDpt, IDEl)
    end if  

end subroutine StressSolid


        subroutine CalculateViscousDamping_interface(ParticleID, IEl)
        !**********************************************************************
        !
        !> Computes a pressure term introducing bulk viscosity damping to the
        !> equation of motion.
        !>
        !! \param[in] ParticleID ID of considered material point.
        !! \param[in] IEl        ID of element of the considered material point.
        !
        !*********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ParticleID, IEl

          real(REAL_TYPE) :: ViscousDampingPressure = 0.0
          real(REAL_TYPE) :: Density = 0.0
          real(REAL_TYPE) :: ElementLMinLocal = 0.0
          real(REAL_TYPE) :: RateVolStrainLocal = 0.0
          real(REAL_TYPE) :: MaterialIndex = 0.0
          real(REAL_TYPE) :: DilationalWaveSpeed = 0.0
          logical :: IsUndrEffectiveStress

          if (.not.CalParams%ApplyBulkViscosityDamping) return

          MaterialIndex = MaterialIDArray(ParticleID)
          
           IsUndrEffectiveStress = &
              !code version 2016 and previous
              ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialIndex)%MaterialType)=='2-phase')) .or. &
              !code version 2017.1 and following
              (trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
           
          if (IsUndrEffectiveStress &
          .or.((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3))) then
            Density = MatParams(MaterialIndex)%DensityMixture / 1000.0
          else
            Density = (1 - MatParams(MaterialIndex)%InitialPorosity) * MatParams(MaterialIndex)%DensitySolid / 1000.0
          end if

          ElementLMinLocal  = ElementLMin(IEl)
          RateVolStrainLocal = RateVolStrain(IEl)

          call GetWaveSpeed(ParticleID, DilationalWaveSpeed)

          ViscousDampingPressure = CalParams%BulkViscosityDamping1 *  &
            Density * DilationalWaveSpeed * ElementLMinLocal * RateVolStrainLocal

          if ((RateVolStrainLocal < 0.0).and.(CalParams%BulkViscosityDamping2 > 0.0)) then
            ViscousDampingPressure = ViscousDampingPressure + &
              Density * (CalParams%BulkViscosityDamping2 * ElementLMinLocal * RateVolStrainLocal)**2
          end if

          Particles(ParticleID)%DBulkViscousPressure = ViscousDampingPressure

        end subroutine CalculateViscousDamping_interface

end module ModExternalSoilModel
