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

!Call the UMAT
        call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
         dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, &
         dfgrd1, noel, npt, layer, kspt, kstep, kinc)

      
!---Definition of Eunloading -> required to define the max time step
      Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))

! PlasticMultiplier: output for plotting plastic points
      PlasticMultiplier = STATEV(2)

        return

    end subroutine ESM_MC

!----------------------------------------------------------
!  CarSig: Transform principal stresses back to Cartesian
!----------------------------------------------------------
    Subroutine CarSig(Sig1, Sig2, Sig3, xN1, xN2, xN3, ntens, SigC)
      implicit double precision (a-h, o-z)
      integer, intent(in) :: ntens
      double precision, intent(in)  :: Sig1, Sig2, Sig3
      double precision, intent(in)  :: xN1(3), xN2(3), xN3(3)
      double precision, intent(out) :: SigC(ntens)
      double precision :: xP1(3,3), xP2(3,3), xP3(3,3), SigCart(3,3)
      integer :: i, j

      ! Construct spectral decomposition: Sig = Sig1*n1 x n1 + Sig2*n2 x n2 + Sig3*n3 x n3
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

      ! Pack back to Voigt notation
      SigC(1) = SigCart(1,1)
      SigC(2) = SigCart(2,2)
      SigC(3) = SigCart(3,3)
      if (ntens >= 4) SigC(4) = SigCart(1,2)
      if (ntens >= 5) SigC(5) = SigCart(2,3)
      if (ntens >= 6) SigC(6) = SigCart(1,3)

      return
    end subroutine CarSig

!----------------------------------------------------------
!  MatTranspose: Transpose a 3x3 matrix
!  Retained for future extensions (e.g., rotating stiffness tensor).
!----------------------------------------------------------
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

    end module ModMohrCoulomb 
