! =============================================================================
! ExternalSoilModel_HB_patch.for
!
! REFERENCE / PATCH FILE
! Shows the modifications needed to src/ExternalSoilModel.for in the
! Anura3D_OpenSource repository to integrate the Hoek-Brown model.
!
! Instructions:
!   1. Open src/ExternalSoilModel.for
!   2. Apply the two changes described below.
! =============================================================================


! -----------------------------------------------------------------------------
! CHANGE 1 — Add USE statement (around line 50, after "use ModBingham")
! -----------------------------------------------------------------------------
!
! BEFORE:
!   use ModBingham
!
! AFTER:
!   use ModBingham
!   use ModHoekBrown
!
! Note: ModHoekBrown is defined in src/Soilmodels/A3DHoekBrown.f


! -----------------------------------------------------------------------------
! CHANGE 2 — Add elseif branch in subroutine StressSolid
!            (after the existing ESM_MOHR_COULOMB_STANDARD branch,
!             around lines 175-183 of ExternalSoilModel.for)
! -----------------------------------------------------------------------------
!
! BEFORE (end of existing MC branch):
!
!   elseif (trim(NameModel)//char(0) == trim(ESM_MOHR_COULOMB_STANDARD)//char(0)) then
!     props(1) = Particles(IDpt)%ShearModulus ! shear modulus, G
!     props(2) = MatParams(IDSet)%PoissonRatio
!     props(3) = SIN(MatParams(IDSet)%FrictionAngle*(Pi/180.0))
!     props(4) = Particles(IDpt)%CohesionCosPhi
!     props(5) = SIN(MatParams(IDSet)%DilatancyAngle*(Pi/180.0))
!     props(6) = MatParams(IDSet)%TensileStrength
!     cmname = UMAT_MOHR_COULOMB_STANDARD
!   endif
!
! AFTER (add the HB branch before "endif"):
!
!   elseif (trim(NameModel)//char(0) == trim(ESM_MOHR_COULOMB_STANDARD)//char(0)) then
!     props(1) = Particles(IDpt)%ShearModulus
!     props(2) = MatParams(IDSet)%PoissonRatio
!     props(3) = SIN(MatParams(IDSet)%FrictionAngle*(Pi/180.0))
!     props(4) = Particles(IDpt)%CohesionCosPhi
!     props(5) = SIN(MatParams(IDSet)%DilatancyAngle*(Pi/180.0))
!     props(6) = MatParams(IDSet)%TensileStrength
!     cmname = UMAT_MOHR_COULOMB_STANDARD
!   elseif (trim(NameModel)//char(0) == trim(ESM_HOEK_BROWN)//char(0)) then
!     props(1) = Particles(IDpt)%ShearModulus          ! G  (MPa)
!     props(2) = MatParams(IDSet)%PoissonRatio          ! nu
!     props(3) = MatParams(IDSet)%ESM_Solid(3)          ! SIGCI (MPa)
!     props(4) = MatParams(IDSet)%ESM_Solid(4)          ! mi
!     props(5) = MatParams(IDSet)%ESM_Solid(5)          ! GSI
!     props(6) = MatParams(IDSet)%ESM_Solid(6)          ! D
!     props(7) = SIN(MatParams(IDSet)%DilatancyAngle*(Pi/180.0))  ! sin(psi)
!     props(8) = MatParams(IDSet)%TensileStrength        ! TENS (MPa)
!     cmname = UMAT_HOEK_BROWN
!   endif


! =============================================================================
! END OF PATCH FILE
! =============================================================================
