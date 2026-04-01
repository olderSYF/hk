! =============================================================================
! ReadMaterialData_HB_patch.FOR
!
! REFERENCE / PATCH FILE
! Shows how to add Hoek-Brown material parameter reading in
! src/ReadMaterialData.FOR of the Anura3D_OpenSource repository.
!
! The Hoek-Brown model uses ESM_Solid(1:8).  ESM_Solid(1) and (2) are the
! standard shear modulus G and Poisson ratio already handled by the base reader.
! This patch adds reading of ESM_Solid(3:8) from the GOM file.
!
! Instructions:
!   1. Open src/ReadMaterialData.FOR
!   2. Locate the section that reads material parameters for ESM models
!      (search for "ESM_MOHR_COULOMB_STANDARD" or "ReadMohrCoulombData").
!   3. Add the HB reading block shown below in the equivalent location.
! =============================================================================


! -----------------------------------------------------------------------------
! GOM file keyword block for Hoek-Brown material
! (Example GOM snippet — see README.md for full example)
!
!   $$MATERIAL_INDEX 1
!   $$MATERIAL_MODEL hoek_brown
!   $$DENSITY_SOLID  2700.0
!   $$YOUNG_MODULUS  10000.0   ! E = 2G(1+nu)
!   $$POISSON_RATIO  0.25
!   $$SHEAR_MODULUS  4000.0    ! G  (MPa)
!   $$SIGCI          50.0      ! Intact rock UCS (MPa)
!   $$MI             10.0      ! Intact rock constant mi
!   $$GSI            60.0      ! Geological Strength Index
!   $$DISTURBANCE    0.0       ! Disturbance factor D
!   $$DILATANCY_ANGLE 5.0      ! dilation angle (degrees)
!   $$TENSILE_STRENGTH 0.0     ! 0 = auto-calculate from HB parameters
! -----------------------------------------------------------------------------


! -----------------------------------------------------------------------------
! FORTRAN code to insert into ReadMaterialData.FOR
! (inside the material reading loop, in the elseif chain for model names)
! -----------------------------------------------------------------------------
!
!   elseif (trim(NameModel)//char(0) == trim(ESM_HOEK_BROWN)//char(0)) then
!
!     ! ESM_Solid(1) = G    — read from $$SHEAR_MODULUS (standard field)
!     ! ESM_Solid(2) = ENU  — read from $$POISSON_RATIO (standard field)
!
!     ! SIGCI — intact rock uniaxial compressive strength (MPa, positive)
!     call ReadRealValue(FileUnit, '$$SIGCI', MatParams(ISet)%ESM_Solid(3), IsError)
!     if (IsError) call GiveError('HB: $$SIGCI not found for material ' // trim(String(ISet)))
!
!     ! mi — intact rock material constant
!     call ReadRealValue(FileUnit, '$$MI', MatParams(ISet)%ESM_Solid(4), IsError)
!     if (IsError) call GiveError('HB: $$MI not found for material ' // trim(String(ISet)))
!
!     ! GSI — Geological Strength Index (10 to 100)
!     call ReadRealValue(FileUnit, '$$GSI', MatParams(ISet)%ESM_Solid(5), IsError)
!     if (IsError) call GiveError('HB: $$GSI not found for material ' // trim(String(ISet)))
!
!     ! D — disturbance factor (0 = undisturbed, 1 = fully disturbed)
!     call ReadRealValue(FileUnit, '$$DISTURBANCE', MatParams(ISet)%ESM_Solid(6), IsError)
!     if (IsError) call GiveError('HB: $$DISTURBANCE not found for material ' // trim(String(ISet)))
!
!     ! Dilation angle: read in degrees, convert to sin(psi) and store
!     call ReadRealValue(FileUnit, '$$DILATANCY_ANGLE', TmpVal, IsError)
!     if (IsError) call GiveError('HB: $$DILATANCY_ANGLE not found for material ' // trim(String(ISet)))
!     MatParams(ISet)%ESM_Solid(7) = SIN(TmpVal * Pi / 180.0)
!
!     ! Tensile strength cutoff (MPa, positive; 0 = auto-calculate from HB)
!     call ReadRealValue(FileUnit, '$$TENSILE_STRENGTH', MatParams(ISet)%ESM_Solid(8), IsError)
!     if (IsError) then
!       MatParams(ISet)%ESM_Solid(8) = 0.0  ! auto-calculate
!     end if
!
!     ! Also store in the standard TensileStrength field for consistency
!     MatParams(ISet)%TensileStrength = MatParams(ISet)%ESM_Solid(8)
!
!     ! Register model ID
!     MatParams(ISet)%MaterialModel = ESM_HOEK_BROWN
!     MatParams(ISet)%ModelID = HOEK_BROWN
!     MatParams(ISet)%SoilModelDLL = UMAT_HOEK_BROWN


! =============================================================================
! END OF PATCH FILE
! =============================================================================
