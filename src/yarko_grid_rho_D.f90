program yarko_grid_rho_D
   use used_const
   use yarko_force
   implicit none
   real(kind=dkind) :: rho, Kmin, Kmax, deltaK, K, C
   real(kind=dkind) :: rhomin, rhomax, deltarho, thermal_inertia
   real(kind=dkind) :: radius
   real(kind=dkind) :: semiaxm, ecc
   real(kind=dkind) :: rotPer
   real(kind=dkind) :: epsi, alpha
   real(kind=dkind) :: yarko_90 
   real(kind=dkind) :: gam, expo
   integer          :: expK
   real(kind=dkind) :: time 
   ! Read input data
   call readData(rhomin, rhomax, deltarho, Kmin, Kmax, C, radius, & 
               & semiaxm, ecc, rotPer, epsi, alpha)
   write(*,*) "=========================="
   write(*,*) "Density min:              ", rhomin
   write(*,*) "Density max:              ", rhomax
   write(*,*) "Density delta:            ", deltarho
   write(*,*) "Thermal conductivity min: ", Kmin
   write(*,*) "Thermal conductivity max: ", Kmax
   write(*,*) "Specific heat capacity:   ", C
   write(*,*) "Radius:                   ", radius
   write(*,*) "Semimajor axis:           ", semiaxm
   write(*,*) "Eccentricity:             ", ecc
   write(*,*) "Rotation period:          ", rotPer
   write(*,*) "Emissivity:               ", epsi
   write(*,*) "Absorption coefficient:   ", alpha
   write(*,*) "=========================="
   open(unit=4, file='output/yarko_grid_rho_D.txt', action='write')    
   ! Set the obliquity to 90 deg
   gam  = 0.13377161828308E+03 
   expo = 0.d0
   ! Take the initial value of density
   rho    = rhomin
   ! Do loop on density
   do while(rho.le.rhomax)
      K = Kmin 
      ! Do loop on thermal conductivity
      do while(K.le.Kmax)
         ! Variable deltaK for a good plot in log scale
         expK = floor(log10(K)) - 1
         deltaK = 0.5*10**real(expK)
         call yarko_eccentric(semiaxm, ecc, rho, K, C, radius, gam, rotPer, alpha, epsi, expo, yarko_90)
         thermal_inertia = sqrt(rho*K*C)
         ! Write results on a file
         write(4,'(4e20.4,1x)') rho, K, thermal_inertia, yarko_90 
         ! Increase the thermal conductivity
         K = K + deltaK
      enddo 
      rho = rho + deltarho
   enddo
   close(4)
end program 

subroutine readData(rhomin, rhomax, deltarho, thermalCondMin, thermalCondMax, heatCap,&
      & radius, semiaxm, ecc, rotPer, absCoeff, emissiv)
   ! This subroutine reads the input data and compute the multiplying factor
   ! of Equation (8), independent from r and r_p
   use used_const
   implicit none
   real(kind=dkind), intent(out) :: rhomin, rhomax, deltarho
   real(kind=dkind), intent(out) :: thermalCondMin, thermalCondMax
   real(kind=dkind), intent(out) ::     heatCap
   real(kind=dkind), intent(out) ::      radius
   real(kind=dkind), intent(out) ::     semiaxm
   real(kind=dkind), intent(out) ::         ecc
   real(kind=dkind), intent(out) ::      rotPer
   real(kind=dkind), intent(out) ::    absCoeff
   real(kind=dkind), intent(out) ::     emissiv
   ! end interface
   namelist /asteroid/ rhomin, rhomax, deltarho, thermalCondMin, thermalCondMax, &
      & heatCap, radius, semiaxm, ecc, rotPer, absCoeff, emissiv
   ! read the input namelist
   open(unit=1,file="input/yarko_grid_rho_D.nml",status="old",action="read")
   read(1,asteroid)
   close(1)
end subroutine readData
