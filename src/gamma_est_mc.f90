!      Author: Marco Fenucci
! Institution: University of Belgrade
!        Date: August 2020
program gamma_est_mc
   use used_const
   use yarko_force, only: yarko_eccentric
   implicit none
   real(kind=dkind) :: rho, C, Kmin, Kmax
   real(kind=dkind) :: radius 
   real(kind=dkind) :: semiaxm, ecc
   real(kind=dkind) :: gam
   real(kind=dkind) :: rotPer
   real(kind=dkind) :: epsi, alpha
   real(kind=dkind) :: time 
   real(kind=dkind) :: levelCurve
   real(kind=dkind) :: expo
   integer :: method 
   integer          :: n_D, n_rho, n_gamma, n_dadt, n_P
   real(kind=dkind) :: Kcross(6)
   integer          :: nCross
   real(kind=dkind), allocatable,  dimension(:) :: diam_mc, rho_mc, gamma_mc, dadt_mc, period_mc
   real(kind=dkind) :: k1_mean, k2_mean
   real(kind=dkind) :: k1_std, k2_std
   real(kind=dkind) :: thermalInertia 
   real(kind=dkind) :: ecc_yarko, pct_err
   character(80)    :: filename
   integer          :: comparison
   integer :: hh, kk, jj, ii, ll, uu
   integer :: iter
   integer, parameter :: max_iter = 1000000
   character(len=1) :: bar, back, dot
   ! Read input data
   call readData(C, Kmin, Kmax, semiaxm, ecc,  alpha, epsi, method, filename, comparison, expo)
   ! Read the distributions of diameter, density and obliquity
   call readLengths(n_D, n_rho, n_gamma, n_dadt, n_P)
   allocate(diam_mc(1:n_D), rho_mc(1:n_rho), gamma_mc(1:n_gamma), dadt_mc(1:n_dadt), period_mc(1:n_P))
   call readMCdata(diam_mc, rho_mc, gamma_mc, dadt_mc, period_mc, n_D, n_rho, n_gamma, n_dadt, n_P)
   write(*,*) "======================================"
   write(*,*) "========== INPUT PARAMETERS =========="
   write(*,*) "======================================"
   write(*,*) "Heat capacity:            ", C
   write(*,*) "Semimajor axis:           ", semiaxm
   write(*,*) "Emissivity:               ", epsi
   write(*,*) "Absorption coefficient:   ", alpha
   write(*,*) "Method:                   ", method 
   write(*,*) "Filename:                 ", filename
   write(*,*) "Comparison:               ", comparison
   write(*,*) "G scaling exponent:       ", expo
   write(*,*) "n_D                       ", n_D
   write(*,*) "n_rho                     ", n_rho
   write(*,*) "n_gamma                   ", n_gamma
   write(*,*) "n_dadt                    ", n_dadt
   write(*,*) "n_P                       ", n_P
   write(*,*) "======================================"
   write(*,*) "Running simulation... "
   time = 1.d6*y2s
   ! Initialize the seed for the generation of random numbers
   call init_random_seed()
   open(unit=10, file='output/'//filename(1:len_trim(filename)),action='write')    
   if(comparison.eq.1)then
      open(unit=11, file='output/diffEst.txt',action='write')    
   endif
!   open(unit=12, file='output/peakLow.txt',action='write')    
!   open(unit=13, file='output/peakHigh.txt',action='write')    
   ! Loop on the distributions of dadt, diameter, density and obliquity
   back = char(8)
   bar  = '='
   dot  = ' '
   do iter=1, max_iter
      write(*,'(256a1)', advance='no') (back, uu =1,30+10)
      write(*,'(1x, 1i3,1a1,2x,1a1,256a1,1a1,256a1,1a1)', advance='no') 100*iter/max_iter,'%','[', &
         (bar, uu =1,30*iter/max_iter), '>', (dot, uu=1,(30-30*iter/max_iter)), ']'
!      write(*,*) "iter ", iter 
      call random_combination(n_D, n_rho, n_gamma, n_dadt, n_P, hh, kk, jj, ii, ll)
      ! Take the values of the parameters. Note that diameter and density comes in couples,
      ! since they are correlated by the albedo 
      radius = diam_mc(hh)/2.d0 
!      rho    = rho_mc(kk)
      rho    = rho_mc(hh)
      gam    = gamma_mc(jj) 
      rotPer = period_mc(ll)
      levelCurve = dadt_mc(ii)
      call yarkoInvert(rho, C, radius, semiaxm, ecc,&
         & gam, rotPer, alpha, epsi, time, levelCurve, kMin, kMax, kCross, nCross, method, expo)
      do ii=1, nCross
         thermalInertia = sqrt(rho*KCross(ii)*C)
         write(10,*)  KCross(ii), thermalInertia, rho, gam
         ! ==========================================
         ! Test to extract the parameters at the
         ! maxima, it does not give very significant
         ! results...
         ! ==========================================
         ! if(abs(thermalInertia-11.36d0).lt.0.1d0)then
         ! if(abs(Kcross(ii)-7.98d-5).lt.1.d-6)then
         !    write(12, *) radius*2, rho, Kcross(ii), gam 
         ! elseif(abs(thermalInertia-87.85d0).lt.0.1d0)then
         ! elseif(abs(Kcross(ii)-0.0046d0).lt.1.d-4)then
         !    write(13, *) radius*2, rho, Kcross(ii), gam 
         ! endif
      enddo

      ! Compute the comparison with eccentric orbit, if needed
      if(comparison.eq.1)then
         do ii=1, nCross
            call yarko_eccentric(semiaxm, ecc, rho, Kcross(ii), C, radius, gam, rotPer, alpha, epsi, expo, ecc_yarko)
            pct_err = abs(ecc_yarko-levelCurve)/max(abs(ecc_yarko), abs(levelCurve))*100.d0
            write(11, *) ecc_yarko, levelCurve, pct_err
            ! ==========================================
            ! Test to extract the parameters at the
            ! two peaks in the errors.
            ! ==========================================
            ! if(abs(pct_err-0.3d0) .lt. 0.1d0)then
            !    write(12, *) radius, rho, Kcross(ii), gam
            ! elseif(abs(pct_err-7.6d0) .lt. 1.00d0)then
            !    write(13, *) radius, rho, Kcross(ii), gam
            ! endif
         enddo
      endif
   enddo
   write(*,*)
   write(*,*) "Done"
   close(10)
   if(comparison.eq.1)then
      close(11)
   endif
!   close(12)
!   close(13)
   deallocate(diam_mc, rho_mc, gamma_mc, dadt_mc, period_mc)
end program gamma_est_mc 

!========================================================
!============== INPUT/OUTPUT SUBROUTINES ================
!========================================================

subroutine readData(C, thermalCondMin, thermalCondMax, &
      &  semiaxm, ecc, absCoeff, emissiv, method, filename, comparison, expo)
   ! This subroutine reads the input data and compute the multiplying factor
   ! of Equation (8), independent from r and r_p
   use used_const
   implicit none
   real(kind=dkind), intent(out) :: C
   real(kind=dkind), intent(out) :: thermalCondMin, thermalCondMax
   real(kind=dkind), intent(out) :: semiaxm, ecc
   real(kind=dkind), intent(out) :: absCoeff
   real(kind=dkind), intent(out) :: emissiv
   real(kind=dkind), intent(out) :: expo
   integer, intent(out) :: method 
   character(80), intent(out)    :: filename
   integer, intent(out)          :: comparison
   ! end interface
   namelist /asteroid/ C, thermalCondMin, thermalCondMax, &
      & semiaxm, ecc, absCoeff, emissiv, method, filename, comparison, expo
   ! read the input namelist
   open(unit=1,file="input/gamma_est_mc.nml",status="old",action="read")
   read(1,asteroid)
   close(1)
end subroutine readData

subroutine readLengths(n_D, n_rho, n_gamma, n_dadt, n_P)
   use used_const
   implicit none
   integer, intent(out) :: n_D, n_rho, n_gamma, n_dadt, n_P
   ! end interface
   real(kind=dkind) :: x
   integer :: flag
   open(unit=1, file='input/diam_mc.txt', status='old')
   n_D = 0
   do 
      read(1,*,iostat=flag) x
      if(flag.lt.0)then
         exit
      endif
      n_D = n_D+1
   enddo
   close(1)

   open(unit=2, file='input/rho_mc.txt', status='old')
   n_rho = 0
   do 
      read(2,*,iostat=flag) x 
      if(flag.lt.0)then
         exit
      endif
      n_rho = n_rho+1
   enddo
   close(2)

   open(unit=3, file='input/gamma_mc.txt', status='old')
   n_gamma = 0
   do 
      read(3,*,iostat=flag) x 
      if(flag.lt.0)then
         exit
      endif
      n_gamma = n_gamma+1
   enddo
   close(3)

   open(unit=4, file='input/dadt_mc.txt', status='old')
   n_dadt = 0
   do 
      read(4,*,iostat=flag) x 
      if(flag.lt.0)then
         exit
      endif
      n_dadt = n_dadt+1
   enddo
   close(4)

   open(unit=5, file='input/period_mc.txt', status='old')
   n_P = 0
   do 
      read(5,*,iostat=flag) x 
      if(flag.lt.0)then
         exit
      endif
      n_P = n_P+1
   enddo
   close(5)

end subroutine readLengths

subroutine readMCdata(diam_mc, rho_mc, gamma_mc, dadt_mc, period_mc, n_D, n_rho, n_gamma, n_dadt, n_P)
   ! Read the input distributions of diameter, density and obliquity
   use used_const
   implicit none
   integer, intent(in)            :: n_D, n_rho, n_gamma, n_dadt, n_P
   real(kind=dkind), intent(out)  :: diam_mc(1:n_D)
   real(kind=dkind), intent(out)  :: rho_mc(1:n_rho)
   real(kind=dkind), intent(out)  :: gamma_mc(1:n_gamma)
   real(kind=dkind), intent(out)  :: dadt_mc(1:n_dadt)
   real(kind=dkind), intent(out)  :: period_mc(1:n_P)
   ! end interface
   integer :: i
   open(unit=1, file='input/diam_mc.txt', status='old')
   do i=1,n_D
      read(1,*) diam_mc(i)
   enddo
   close(1)

   open(unit=2, file='input/rho_mc.txt', status='old')
   do i=1,n_rho
      read(2,*) rho_mc(i)
   enddo
   close(2)

   open(unit=3, file='input/gamma_mc.txt', status='old')
   do i=1,n_gamma
      read(3,*) gamma_mc(i)
   enddo
   close(3)

   open(unit=4, file='input/dadt_mc.txt', status='old')
   do i=1,n_dadt
      read(4,*) dadt_mc(i)
   enddo
   close(4)

   open(unit=5, file='input/period_mc.txt', status='old')
   do i=1,n_P
      read(5,*) period_mc(i)
   enddo
   close(5)
end subroutine readMCdata

!========================================================
!=========== INVERSION OF YARKOVSKY DRIFT ===============
!========================================================

subroutine yarkoInvert(rho, C, radius, semiaxm, ecc, gam, rotPer, alpha, epsi, &
      & time, levelCurve, kMin, kMax, kCrosses, nCross, method, expo)
   use used_const
   use yarko_force
   implicit none
   real(kind=dkind), intent(in)  :: rho, C, radius, semiaxm, ecc, gam, rotPer
   real(kind=dkind), intent(in)  :: alpha, epsi, time, levelCurve
   real(kind=dkind), intent(in)  :: Kmin, Kmax
   real(kind=dkind), intent(in)  :: expo
   real(kind=dkind), intent(out) :: kCrosses(6)
   integer, intent(out)          :: nCross
   integer, intent(in)           :: method
   ! end interface
   real(kind=dkind) :: K, deltaK
   real(kind=dkind) :: expK
   real(kind=dkind) :: KCross
   real(kind=dkind) :: yarko
   logical :: flag, flagTmp
   nCross   = 0
   kCrosses = 0.d0
   K      = Kmin 
   ! Take the initial sign to see on which side of
   ! the level curve are we
   if(method.eq.1)then
      call computeYarko_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, yarko)
   elseif(method.eq.2)then
       call yarko_eccentric(semiaxm, ecc, rho, K, C, radius, gam, rotPer, alpha, epsi, expo, yarko)
   else
      call computeYarko_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, time, yarko)
   endif
   if(yarko-levelCurve.gt.0)then
      flag = .true.
   else
      flag = .false.
   endif
   ! Start the loop on the thermal conductivity
   do while(K.le.Kmax)
      ! Variable deltaK for a good logarithmic discretization 
      expK = floor(log10(K))-1
!      deltaK = 0.5*10**real(expK)
      deltaK = 5.d0*10**real(expK)
      ! Compute the Yarkovsky drift
      if(method.eq.1)then
         call computeYarko_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, yarko)
      elseif(method.eq.2)then
         call yarko_eccentric(semiaxm, ecc, rho, K, C, radius, gam, rotPer, alpha, epsi, expo, yarko)
      else
         call computeYarko_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, time, yarko)
      endif
      ! Check if we crossed the level curve 
      if(yarko-levelCurve.gt.0)then
         flagTmp = .true.
      else
         flagTmp = .false.
      endif
      ! If we crossed the level curve, start the bisection method to compute 
      ! the solution more accurately
      if(flag.neqv.flagTmp)then
         ! We have crossed the level curve!
         ! Restore the flag and start a bisection 
         ! method to find the exact point
         flag   = flagTmp
         call bisectionMethod(rho, C, radius, semiaxm, ecc, gam, rotPer, alpha, epsi, time, &
            & K-deltaK, K, levelCurve, KCross, method, expo)
         nCross = nCross + 1
         KCrosses(nCross) = kCross
      endif
      K = K+deltaK
   enddo 
end subroutine yarkoInvert

subroutine bisectionMethod(rho, C, radius, semiaxm, ecc, gam, rotPer, alpha, epsi, time, Ka, Kb, levelCurve, KCross, method, expo)
   use used_const
   use yarko_force
   implicit none
   real(kind=dkind), intent(in)  :: rho, C, radius, semiaxm, ecc, gam, rotPer, alpha, epsi, time
   real(kind=dkind), intent(in)  :: Ka, Kb, levelCurve
   real(kind=dkind), intent(in)  :: expo 
   real(kind=dkind), intent(out) :: KCross
   integer, intent(in)           :: method
   ! end interface
   real(kind=dkind)            :: yarko
   real(kind=dkind)            :: K1, K2, Kmean
   real(kind=dkind)            :: FK1, FK2, FKmean
   integer                     :: n
   integer, parameter          :: nmax = 100 
   real(kind=dkind), parameter :: tol = 1e-11
   ! set an absurde value for KCross
   KCross = 1d9
   ! Initialize the points
   K1 = Ka
   K2 = Kb
   ! F = yarko-levelCurve
   if(method.eq.1)then
     call computeYarko_vokrouhlicky(rho, K1, C, radius, semiaxm, gam, rotPer, alpha, epsi, yarko)
   elseif(method.eq.2)then
     call yarko_eccentric(semiaxm, ecc, rho, K1, C, radius, gam, rotPer, alpha, epsi, expo, yarko)
   else
     call computeYarko_xu2020(rho, K1, C, radius, semiaxm, gam, rotPer, alpha, epsi, time, yarko)
   endif
   FK1 = yarko-levelCurve
   if(method.eq.1)then
      call computeYarko_vokrouhlicky(rho, K2, C, radius, semiaxm, gam, rotPer, alpha, epsi, yarko)
   elseif(method.eq.2)then
      call yarko_eccentric(semiaxm, ecc, rho, K2, C, radius, gam, rotPer, alpha, epsi, expo, yarko)
   else
      call computeYarko_xu2020(rho, K2, C, radius, semiaxm, gam, rotPer, alpha, epsi, time, yarko)
   endif
   FK2 = yarko-levelCurve
   do n=1, nmax
      Kmean = (K1+K2)/2.d0
      if(method.eq.1)then
         call computeYarko_vokrouhlicky(rho, Kmean, C, radius, semiaxm, gam, rotPer, alpha, epsi, yarko)
      elseif(method.eq.2)then
         call yarko_eccentric(semiaxm, ecc, rho, Kmean, C, radius, gam, rotPer, alpha, epsi, expo, yarko)
      else
         call computeYarko_xu2020(rho, Kmean, C, radius, semiaxm, gam, rotPer, alpha, epsi, time, yarko)
      endif
      FKmean = yarko-levelCurve
      if((K2-K1)/2.d0 .lt. tol)then
         KCross = Kmean
         exit
      endif

      if(FK1*FKmean.gt.0.d0)then
         K1 = Kmean
      else
         K2 = Kmean
      endif
   enddo
   if(n.eq.nmax)then
      write(*,*) "Warning!! Bisection method failed!!"
      stop
   endif
end subroutine bisectionMethod

!========================================================
!============= TO GENERATE RANDOM NUMBERS ===============
!========================================================

subroutine random_combination(n_D, n_rho, n_gamma, n_dadt, n_P, rand_D, rand_rho, rand_gamma, rand_dadt, rand_P)
   use used_const
   implicit none
   integer, intent(in)  :: n_D, n_rho, n_gamma, n_dadt, n_P
   integer, intent(out) :: rand_D, rand_rho, rand_gamma, rand_dadt, rand_P
   ! end interface
   real(kind=dkind) :: u
   call random_number(u)
   rand_D = 1 + FLOOR((n_D)*u)
   call random_number(u)
   rand_rho = 1 + FLOOR((n_rho)*u)
   call random_number(u)
   rand_gamma = 1 + FLOOR((n_gamma)*u)
   call random_number(u)
   rand_dadt = 1 + FLOOR((n_dadt)*u)
   call random_number(u)
   rand_P  = 1 + FLOOR((n_P)*u)
end subroutine random_combination

! Initial seed for the random number generator
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
