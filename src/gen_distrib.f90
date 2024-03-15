! ################################################################
! #     Authors: Marco Fenucci, Bojan Novaković, Dušan Marceta   #
! # Institution: University of Belgrade                          #
! #        Date: July 2022                                       #
! ################################################################
!
! PURPOSE: 
!
program gen_distrib 
   use iso_fortran_env, only: output_unit
   use used_const
   use population_model 
   implicit none
   ! Orbital elements of the NEO
   real(kind=dkind) :: sma, ecc, inc
   ! Absolute magnitude
   real(kind=dkind) :: mean_H,    std_H   
   ! Albedo
   real(kind=dkind) :: mean_pv,   std_pv  
   ! Bulk density     [kg m^-3]
   real(kind=dkind) :: mean_rho,  std_rho 
   ! Diameter         [m]
   real(kind=dkind) :: mean_D,    std_D   
   ! Surface density  [kg m^-3]
   real(kind=dkind) :: mean_rhos, std_rhos 
   ! Obliquity        [deg]
   real(kind=dkind) :: mean_gam,  std_gam  
   ! Rotation period  [h]
   real(kind=dkind) :: mean_P,    std_P    
   ! Yarkovsky drift  [au My^-1]
   real(kind=dkind) :: mean_dadt, std_dadt 
   ! Absorption coefficient
   real(kind=dkind) :: alpha, aa, bond_albedo
   ! Phase integral, use G = 0.15
   real(kind=dkind), parameter :: qq = 0.29d0 + 0.685d0*0.15d0
   ! Probabilites for source routes
   real(kind=dkind) :: proute(7)
   ! 
   real(kind=dkind) :: rho, D, pv, rhos, gam, P, dadt
   ! Define the units for output files
   integer, parameter :: unit_pv    = 100
   integer, parameter :: unit_rho   = 101
   integer, parameter :: unit_D     = 102
   integer, parameter :: unit_rhos  = 103
   integer, parameter :: unit_gam   = 104
   integer, parameter :: unit_P     = 105
   integer, parameter :: unit_dadt  = 106
   integer, parameter :: unit_alpha = 107
   ! Loop variables
   integer            :: iter, max_iter
   ! Initialize seed for random number generator
   call init_random_seed()
   ! Load the Granvik et al. 2018 NEO population model
   call gmb_load()
   ! Read the input file gen_distri.nml
   call read_data(max_iter, sma, ecc, inc, mean_H, std_H, mean_pv, std_pv, & 
      & mean_rho, std_rho, mean_D, std_D, mean_rhos, std_rhos, mean_gam, std_gam, &
      & mean_P, std_P, mean_dadt, std_dadt, alpha)
   call print_data(sma, ecc, inc, mean_H, std_H, mean_pv, std_pv, & 
      & mean_rho, std_rho, mean_D, std_D, mean_rhos, std_rhos, mean_gam, std_gam, &
      & mean_P, std_P, mean_dadt, std_dadt, alpha)
   ! Compute the probabilities from Granvik et al 2018 model
   call gmb_search(sma, ecc, inc, mean_H, proute)
   ! Open the files for saving the distributions
   open(unit=unit_pv,    file='input/pv_mc.txt',       action='write')
   open(unit=unit_rho,   file='input/rho_mc.txt',      action='write')
   open(unit=unit_D,     file='input/diam_mc.txt',     action='write')
   open(unit=unit_rhos,  file='input/rho_surf_mc.txt', action='write')
   open(unit=unit_gam,   file='input/gamma_mc.txt',    action='write')
   open(unit=unit_P,     file='input/period_mc.txt',   action='write')
   open(unit=unit_dadt,  file='input/dadt_mc.txt',     action='write')
   open(unit=unit_alpha, file='input/alpha_mc.txt',    action='write')
   do iter=1, max_iter
      call progress_bar(iter, max_iter)
      ! There are 4 cases
      ! pV  D  rho
      ! 0   0    0
      ! 0   0    1
      ! 1   1    0 
      ! 1   1    1
      if(mean_pv.lt.0.d0 .and. mean_D.lt.0.d0 .and. mean_rho.lt.0.d0)then
         ! 0 0 0: all population based
         call gen_density_diameter(mean_H, std_H, proute, rho, D, pv)
      elseif(mean_pv.lt.0.d0 .and. mean_D.lt.0.d0 .and. mean_rho.gt.0.d0)then
         ! 0 0 1: pv and D population based, rho Gaussian
         call gen_diameter(mean_H, std_H, proute, D, pv)
         rho = gen_normal(mean_rho, std_rho)
      elseif(mean_pv.gt.0.d0 .and. mean_D.gt.0.d0 .and. mean_rho.lt.0.d0)then
         ! 1 1 0: pv and D Gaussian
         call gen_density(mean_H, std_H, proute, rho, pv)
         pv = gen_normal(mean_pv, std_pv)
         D  = gen_normal(mean_D,  std_D)
      elseif(mean_pv.gt.0.d0 .and. mean_D.gt.0.d0 .and. mean_rho.gt.0.d0)then
         ! 1 1 1: all gaussian
         pv  = gen_normal(mean_pv,  std_pv)
         D   = gen_normal(mean_D,   std_D)
         rho = gen_normal(mean_rho, std_rho)
      else
         ! This case should never occur
         write(output_unit, *) "ERROR. Wrong selection of parameters combination."
         write(output_unit, *) " STOPPING PROGRAM"
         stop
      endif
      ! rho_surface is always Gaussian
      rhos = gen_normal(mean_rhos, std_rhos)
      if(mean_gam.lt.0.d0)then
         ! Population based
         call gen_obliquity(gam) 
      else
         ! Gaussian distribution
         gam = gen_normal(mean_gam, std_gam)
      endif
      ! Rotation period is always Gaussian
      P = gen_normal(mean_P, std_P)
      ! dadt is always Gaussian
      dadt = gen_normal(mean_dadt, std_dadt)
      ! Write the values in the corresponding files
      write(unit_pv,   *) pv
      write(unit_rho,  *) rho
      write(unit_D,    *) D
      write(unit_rhos, *) rhos
      write(unit_gam,  *) gam
      write(unit_P,    *) P
      write(unit_dadt, *) dadt
      ! If negative, use population model
      if(alpha.lt.0.d0)then
         bond_albedo = qq*pv
         aa          = 1.d0 - bond_albedo
         write(unit_alpha, *) aa
      endif
   enddo
   ! If alpha is positive, use single value
   if(alpha.gt.0.d0)then
      write(unit_alpha, *) alpha
   endif
   ! Close the files
   close(unit_pv)
   close(unit_rho)
   close(unit_D)
   close(unit_rhos)
   close(unit_gam)
   close(unit_P)
   close(unit_dadt)
   close(unit_alpha)
   ! ending program
   write(output_unit,*)
   write(output_unit,*) "Done"
end program gen_distrib

! PURPOSE: 
subroutine read_data(max_iter, sma, ecc, inc, mean_H, std_H, mean_pv, std_pv, & 
      & mean_rho, std_rho, mean_D, std_D, mean_rhos, std_rhos, mean_gam, std_gam, &
      & mean_P, std_P, mean_dadt, std_dadt, alpha)
   use used_const
   implicit none
   integer,          intent(out) :: max_iter
   real(kind=dkind), intent(out) :: sma, ecc, inc
   real(kind=dkind), intent(out) :: mean_H,    std_H   
   real(kind=dkind), intent(out) :: mean_pv,   std_pv  
   real(kind=dkind), intent(out) :: mean_rho,  std_rho 
   real(kind=dkind), intent(out) :: mean_D,    std_D   
   real(kind=dkind), intent(out) :: mean_rhos, std_rhos 
   real(kind=dkind), intent(out) :: mean_gam,  std_gam  
   real(kind=dkind), intent(out) :: mean_P,    std_P    
   real(kind=dkind), intent(out) :: mean_dadt, std_dadt 
   real(kind=dkind), intent(out) :: alpha  
   ! end interface
   namelist /param_model/ max_iter, sma, ecc, inc, mean_H, std_H, mean_pv, std_pv, & 
      & mean_rho, std_rho, mean_D, std_D, mean_rhos, std_rhos, mean_gam, std_gam, &
      & mean_P, std_P, mean_dadt, std_dadt, alpha
   ! read the input namelist
   open(unit=1,file="input/gen_distrib.nml",status="old",action="read")
   read(1,param_model)
   close(1)
   ! Check for errors in the input file
   if(mean_H.lt.0.d0)then
      write(*,*) "ERROR. Absolute magnitude must be positive." 
      write(*,*) " Selected H: ", mean_H
      write(*,*) " STOPPING PROGRAM"
      stop
   elseif(mean_P.lt.0.d0)then
      write(*,*) "ERROR. The rotation period must be positive."
      write(*,*) " Selected P: ", mean_P 
      write(*,*) " STOPPING PROGRAM"
      stop
   elseif(sma.lt.0.d0)then
      write(*,*) "ERROR. The semi-major axis must be positive. " 
      write(*,*) " Selected a: ", sma 
      write(*,*) " STOPPING PROGRAM"
      stop
   elseif(ecc.gt.1.d0 .or. ecc.le.0)then
      write(*,*) "ERROR. Eccentricity must be between 0 and 1."
      write(*,*) " Selected ecc: ", ecc
      write(*,*) " STOPPING PROGRAM"
      stop
   elseif(mean_rhos.lt.0.d0)then
      mean_rhos = 1200.d0
      std_rhos  = 200.d0
   endif
end subroutine

! PURPOSE: 
subroutine print_data(sma, ecc, inc, mean_H, std_H, mean_pv, std_pv, & 
      & mean_rho, std_rho, mean_D, std_D, mean_rhos, std_rhos, mean_gam, std_gam, &
      & mean_P, std_P, mean_dadt, std_dadt, alpha)
   use used_const
   use iso_fortran_env, only: output_unit
   implicit none
   real(kind=dkind), intent(in) :: sma, ecc, inc
   real(kind=dkind), intent(in) :: mean_H,    std_H   
   real(kind=dkind), intent(in) :: mean_pv,   std_pv  
   real(kind=dkind), intent(in) :: mean_rho,  std_rho 
   real(kind=dkind), intent(in) :: mean_D,    std_D   
   real(kind=dkind), intent(in) :: mean_rhos, std_rhos 
   real(kind=dkind), intent(in) :: mean_gam,  std_gam  
   real(kind=dkind), intent(in) :: mean_P,    std_P    
   real(kind=dkind), intent(in) :: mean_dadt, std_dadt 
   real(kind=dkind), intent(in) :: alpha  
   ! end interface
   ! Write the input parameters on screen
   write(output_unit, *) " ====== INPUT SETTINGS ====== "
   write(output_unit, *) "                              "
   write(output_unit, *) " ORBITAL PARAMETERS:          "
   write(output_unit, *) "   Semimajor axis (au) = ", sma
   write(output_unit, *) "   Eccentricity        = ", ecc
   write(output_unit, *) "   Inclination (deg)   = ", inc
   write(output_unit, *) "                              "
   write(output_unit, *) " PHYSICAL PARAMETERS:         "
   write(output_unit, *) "     H = ", mean_H,  "±", std_H
   if(mean_pv.lt.0.d0)then
      write(output_unit, *) "   p_V =    POPULATION MODEL     "
   else
      write(output_unit, *) "   p_V = ", mean_pv, "±", std_pv
   endif
   if(mean_rho.lt.0.d0)then
      write(output_unit, *) "   rho =    POPULATION MODEL     "
   else
      write(output_unit, *) "   rho = ", mean_rho, "±", std_rho
   endif
   if(mean_D.lt.0.d0)then
      write(output_unit, *) "     D =    POPULATION MODEL     "
   else
      write(output_unit, *) "     D = ", mean_D, "±", std_D
   endif
   write(output_unit, *) " rho_s = ", mean_rhos, "±", std_rhos
   if(mean_gam.lt.0.d0)then
      write(output_unit, *) " gamma =    POPULATION MODEL     "
   else
      write(output_unit, *) " gamma = ", mean_gam, "±", std_gam
   endif
   if(alpha.lt.0.d0)then
      write(output_unit, *) " alpha =    POPULATION MODEL     "
   else
      write(output_unit, *) " alpha = ", alpha 
   endif
   write(output_unit, *) "     P = ", mean_P, "±", std_P
   write(output_unit, *) " da/dt = ", mean_dadt, "±", std_dadt
   write(output_unit, *) "                              "
   write(output_unit, *) " ============================ "
   write(output_unit, *) "                              "
end subroutine

! PURPOSE: Initializes the seed for the random number generator
subroutine init_random_seed()
      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed
      call random_seed(size = n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)
      deallocate(seed)
end subroutine init_random_seed

! PURPOSE: Print a progressbar on screen 
!
! INPUT:
!     iter : current iteration number
! max_iter : maximum number of iterations
! 
! OUTPUT:
subroutine progress_bar(iter, max_iter)
   use iso_fortran_env, only: output_unit
   use used_const
   integer, intent(in) :: iter, max_iter
   ! end interface
   integer          :: uu
   character(len=1) :: bar, back, dot
   back = char(8)
   bar  = '='
   dot  = ' '
   write(output_unit,'(256a1)', advance='no') (back, uu =1,30+10)
   flush(output_unit)
   write(output_unit,'(1x, 1i3,1a1,2x,1a1,256a1,1a1,256a1,1a1)', advance='no') 100*iter/max_iter,'%','[', &
      (bar, uu =1,30*iter/max_iter), '>', (dot, uu=1,(30-30*iter/max_iter)), ']'
   flush(output_unit)
end 
