! This program is just a test, it will be deprecated when testing phase is finished
program yarko_compare
   use used_const
   use yarko_force
   implicit none
   real(kind=dkind) :: semiaxm, ecc, R
   real(kind=dkind) :: rho, K, C, yarko_ell2
   real(kind=dkind) :: gam, rotPer
   real(kind=dkind) :: alpha, epsi, kep(6)
   real(kind=dkind) :: yarko, yarko_ell, expo, yarko_quad, start, finish, yarko_v
   call readData(C, kep,  alpha, epsi)
   ! Read input data
   write(*,*) "======================================"
   write(*,*) "========== INPUT PARAMETERS =========="
   write(*,*) "======================================"
   write(*,*) "Semimajor axis:           ", kep(1) 
   write(*,*) "Eccentricity:             ", kep(2) 
   write(*,*) "Heat capacity:            ", C 
   write(*,*) "Rotation period:          ", rotPer
   write(*,*) "Absorption coefficient:   ", alpha
   write(*,*) "Emissivity:               ", epsi
   write(*,*) "======================================"
   expo   = 0.d0
   rho    = 2400.d0
   K      = 0.001d0
   R      = 1000.d0
   gam    = 160.d0
   rotPer = 4.d0
   write(*,*) "Par: ", kep(1), kep(2), rho, K, C, R, gam, rotPer, alpha, epsi, expo
   call computeYarko_circular(rho, K, C, R, kep(1), gam, rotPer, alpha, epsi, yarko_v)
   call yarko_eccentric(kep(1), kep(2), rho, K, C, R, gam, rotPer, alpha, epsi, expo, yarko_ell)
!   call yarko_eccentric2(kep(1), kep(2), rho, K, C, R, gam, rotPer, alpha, epsi, expo, yarko_ell2)
   write(*,*) "Yarko circular       ", yarko_v
   write(*,*) "Yarko eccentric      ", yarko_ell
   write(*,*) "Yarko eccentric2     ", yarko_ell2
end program yarko_compare

!========================================================
!============== INPUT/OUTPUT SUBROUTINES ================
!========================================================

subroutine readData(C, kep, absCoeff, emissiv)
   use used_const
   implicit none
   real(kind=dkind), intent(out) :: C
   real(kind=dkind), intent(out) :: kep(6)
   real(kind=dkind), intent(out) :: absCoeff
   real(kind=dkind), intent(out) :: emissiv
   ! end interface
   real(kind=dkind) :: semiaxm, ecc, inc, OmNod, omega
   real(kind=dkind) :: thermalCondMin, thermalCondMax
   real(kind=dkind) :: expo
   integer          :: method 
   character(80)    :: filename
   integer          :: max_iter
   namelist /asteroid/ C, thermalCondMin, thermalCondMax, &
      & semiaxm, ecc,  absCoeff, emissiv, method, filename, max_iter, expo
   ! read the input namelist
   open(unit=1,file="input/gamma_est_mc.nml",status="old",action="read")
   read(1,asteroid)
   close(1)
   kep(1) = semiaxm
   kep(2) = ecc 
   kep(3) = 0.d0 
   kep(4) = 0.d0 
   kep(5) = 0.d0 
   kep(6) = 0.d0 
end subroutine readData


