!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
      INCLUDE 'ABA_PARAM.INC'

!**********************************************************************
!  UMAT_MohrCoulombStandard.f
!  Hoek-Brown constitutive model for Anura3D (MPM) via ESM "借壳" approach
!
!  PROPS mapping (passed from ExternalSoilModel.for via MC channel):
!    PROPS(1) = G        -- shear modulus
!    PROPS(2) = nu       -- Poisson's ratio
!    PROPS(3) = sigma_ci -- uniaxial compressive strength of intact rock
!    PROPS(4) = GSI      -- Geological Strength Index
!    PROPS(5) = mi       -- intact rock parameter
!    PROPS(6) = D        -- disturbance factor (0-1)
!
!  STATEV mapping:
!    STATEV(1) = accumulated equivalent plastic strain
!    STATEV(2) = plastic flag (0=elastic, 1=shear, 2=tensile)
!    STATEV(3) = mb  (derived, for post-processing)
!    STATEV(4) = s   (derived, for post-processing)
!    STATEV(5) = a   (derived, for post-processing)
!
!  Sign convention: compression positive (consistent with Anura3D)
!**********************************************************************

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      implicit double precision (a-h, o-z)
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

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

!---Compute Hoek-Brown derived parameters (2002 generalized HB criterion, valid for GSI > 25)
!   Reference: Hoek, E., Carranza-Torres, C. & Corkum, B. (2002). Proc. NARMS-TAC, 267-273.
!   mb = mi * exp((GSI - 100) / (28 - 14*D))
!   s  = exp((GSI - 100) / (9 - 3*D))
!   a  = 0.5 + (1/6)*(exp(-GSI/15) - exp(-20/3))
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
!   PrnSig(mode, NTENS, SigVoigt, xN1, xN2, xN3, Sig1, Sig2, Sig3, P, Q)
!   Returns Sig1 <= Sig2 <= Sig3 (Sig3 = most compressive, compression positive)
      call PrnSig(1, NTENS, SigE, xN1, xN2, xN3, Sig1, Sig2, Sig3, P, Q)

!---Evaluate Hoek-Brown yield function
!   F = Sig3 - Sig1 - sigci*(mb*Sig1/sigci + s)^a
!   (Sig1 = minor principal = most tensile; Sig3 = major principal = most compressive)
      plastic_flag = 0
      tension_cut  = .false.

      hb_term = mb*Sig1/sigci + s

!---Tension cutoff: HB is not defined when hb_term < 0
      if (hb_term < 0.0d0) then
        tension_cut = .true.
        ! Tension cutoff stress (Sig1 = sig_t where F=0 and Sig1=Sig3 -> deviatoric=0)
        ! Simple tension cutoff: Sig1 = -s*sigci/mb
        sig_t = -s*sigci/mb
        Sig1_new = sig_t
        Sig3_new = Sig3
        ! Re-check if Sig3 also needs correction
        if (Sig3_new < Sig1_new) Sig3_new = Sig1_new
        plastic_flag = 2

        ! After tension cutoff, re-check shear yield with corrected sig1
        hb_term = mb*Sig1_new/sigci + s
        if (hb_term >= 0.0d0) then
          F = Sig3_new - Sig1_new - sigci*(hb_term**a)
          if (F > FTOL) then
            ! Shear return also needed after tension cutoff
            ! Use Newton-Raphson from corrected state
            Sig1_iter = Sig1_new
            Sig3_iter = Sig3_new
            do iter = 1, MAXITER
              hb_term  = mb*Sig1_iter/sigci + s
              if (hb_term <= 0.0d0) exit
              F_iter   = Sig3_iter - Sig1_iter - sigci*(hb_term**a)
              if (abs(F_iter) <= FTOL) exit
              ! dF/dSig3 = 1,  dF/dSig1 = -1 - sigci*a*(hb_term^(a-1))*(mb/sigci)
              dFdSig3  = 1.0d0
              dFdSig1  = -1.0d0 - a*mb*(hb_term**(a - 1.0d0))
              ! Associated flow: dg/dSig1 = dF/dSig1, dg/dSig3 = dF/dSig3
              ! Sig3_iter = Sig3_new - dLambda*dFdSig3*twoG (no change on elastic predictor for Sig3 in tension cut branch)
              ! Increment using consistent tangent approximation
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
          ! Hit tension regime during iteration -> switch to tension cutoff
          Sig1_iter = -s*sigci/mb
          plastic_flag = 2
          exit
        end if
        F_iter  = Sig3_iter - Sig1_iter - sigci*(hb_term**a)
        if (abs(F_iter) <= FTOL) exit
        ! Partial derivatives of yield function
        dFdSig3 = 1.0d0
        dFdSig1 = -1.0d0 - a*mb*(hb_term**(a - 1.0d0))
        ! Consistent tangent denominator (associated flow rule)
        dF_iter = dFdSig1*(-twoG*dFdSig1) + dFdSig3*(twoG*dFdSig3)
        if (abs(dF_iter) < 1.0d-20) exit
        dLambda   = F_iter / dF_iter
        Sig1_iter = Sig1_iter - dLambda*(-twoG*dFdSig1)
        Sig3_iter = Sig3_iter - dLambda*(twoG*dFdSig3)
      end do

      Sig1_new = Sig1_iter
      Sig3_new = Sig3_iter

!---Back-transformation to Cartesian stress tensor (label 700)
  700 continue

      ! Update equivalent plastic strain (deviatoric measure combining both principal plastic increments)
      epsp_new = STATEV(1) + sqrt((abs(Sig1_new - Sig1)**2 + abs(Sig3_new - Sig3)**2) / 2.0d0) &
                             / max(twoG, 1.0d-20)
      STATEV(1) = epsp_new
      STATEV(2) = dble(plastic_flag)

      ! Sig2 stays the same (no plastic correction on intermediate principal stress
      ! in standard HB return mapping)
      call CarSig(Sig1_new, Sig2, Sig3_new, xN1, xN2, xN3, NTENS, SigC)

      do i = 1, NTENS
        STRESS(i) = SigC(i)
      end do

      return
      END SUBROUTINE UMAT

!**********************************************************************
!  CarSig: Transform principal stresses back to Cartesian (Voigt) form
!  (standalone version, outside module, for use by UMAT)
!  NOTE: A module-internal version also exists in A3DMohrCoulombStandard.f.
!  Both are needed because UMAT is a standalone subroutine (outside any module)
!  and cannot access module-private procedures.
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

!**********************************************************************
!  MatTranspose: Transpose a 3x3 matrix (standalone version)
!  Retained for potential use in future extensions (e.g., rotating the
!  elastic stiffness tensor into the material frame). Not called in the
!  current HB implementation.
!**********************************************************************
      Subroutine MatTranspose(A, AT)
      implicit double precision (a-h, o-z)
      double precision, intent(in)  :: A(3,3)
      double precision, intent(out) :: AT(3,3)
      integer :: i, j
      do i = 1, 3
        do j = 1, 3
          AT(i,j) = A(j,i)
        end do
      end do
      return
      end subroutine MatTranspose
