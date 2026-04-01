! =============================================================================
! GlobalConstants_HB_patch.FOR
!
! REFERENCE / PATCH FILE
! Shows the additions needed to src/GlobalConstants.FOR in the
! Anura3D_OpenSource repository to register the Hoek-Brown model.
!
! Instructions:
!   1. Open src/GlobalConstants.FOR
!   2. Apply the three additions described below.
! =============================================================================


! -----------------------------------------------------------------------------
! ADDITION 1 — Model ID constant
!              Insert after line 181 (MOHR_COULOMB_TEUNISSEN = 109)
! -----------------------------------------------------------------------------
!
!   integer(INTEGER_TYPE), parameter :: HOEK_BROWN = 110
!
! Context (surrounding lines for reference):
!
!   integer(INTEGER_TYPE), parameter :: MOHR_COULOMB_TEUNISSEN = 109
!   integer(INTEGER_TYPE), parameter :: HOEK_BROWN = 110          ! <-- ADD THIS
!   integer(INTEGER_TYPE), parameter :: External = 149


! -----------------------------------------------------------------------------
! ADDITION 2 — Model name string constant (ESM name used in GOM file)
!              Insert after line 187 (ESM_MOHR_COULOMB_STRAIN_SOFTENING)
! -----------------------------------------------------------------------------
!
!   character(len=64), parameter :: ESM_HOEK_BROWN = 'hoek_brown'
!
! Context:
!
!   character(len=64), parameter :: ESM_MOHR_COULOMB_STRAIN_SOFTENING = 'strain_softening_mohr_coulomb'
!   character(len=64), parameter :: ESM_HOEK_BROWN = 'hoek_brown'    ! <-- ADD THIS


! -----------------------------------------------------------------------------
! ADDITION 3 — DLL name string constant
!              Insert after line 211 (UMAT_MOHR_COULOMB_STRAIN_SOFTENING)
! -----------------------------------------------------------------------------
!
!   character(len=64), parameter :: UMAT_HOEK_BROWN = 'A3DHoekBrown.dll'
!
! Context:
!
!   character(len=64), parameter :: UMAT_MOHR_COULOMB_STRAIN_SOFTENING = 'A3DMohrCoulombStrainSoftening.dll'
!   character(len=64), parameter :: UMAT_HOEK_BROWN = 'A3DHoekBrown.dll'  ! <-- ADD THIS


! =============================================================================
! END OF PATCH FILE
! =============================================================================
