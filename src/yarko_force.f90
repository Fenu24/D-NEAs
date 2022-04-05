module yarko_force
   use used_const
   implicit none
   private ::  Fnu_eval,  compute_c1s_c1d, compute_c2s_c2d, prvec, kep2car, keplereq, trapezoid_average

   public :: computeYarko_vokrouhlicky, computeYarkoMaxMin_vokrouhlicky, &
      & yarko_seasonal_vokrouhlicky, yarko_diurnal_vokrouhlicky, &
      & computeYarko_xu2020, computeYarkoMaxMin_xu2020, yarko_seasonal_xu2020, yarko_diurnal_xu2020, &
      & yarko_eccentric, yarkovsky_vf 
           
   contains

   !========================================================
   !====== ANALYTICAL ESTIMATION OF YARKOVSKY EFFECT =======
   !========================================================

   subroutine computeYarkoMaxMin_xu2020(rho, K, C, radius, semiaxm, rotPer, alpha, epsi, time, yarkomax, yarkomin, gammin)
      real(kind=dkind), intent(in)  :: rho, K, C, radius, semiaxm, rotPer, alpha, epsi, time
      real(kind=dkind), intent(out) :: yarkomax, yarkomin
      real(kind=dkind), intent(out), optional :: gammin
      ! end interface
      real(kind=dkind) :: yarkoAdd
      real(kind=dkind) :: gam
      real(kind=dkind) :: f0, f90, f180
      real(kind=dkind) :: dad, das
      gam = 0.d0
      call yarko_seasonal_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, das)
      call yarko_diurnal_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, dad)
      yarkomax=(das+dad)*time

      f0 = yarkomax

      gam = 90.d0
      call yarko_seasonal_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, das)
      call yarko_diurnal_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, dad)
      f90 = (das+dad)*time

      gam = 180.d0
      call yarko_seasonal_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, das)
      call yarko_diurnal_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, dad)
      f180 = (das+dad)*time

      if(abs(f0/(2.d0*f90)).lt.1.d0)then
         gam = acos(f0/(2.d0*f90))*rad2deg 
         call yarko_seasonal_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, das)
         call yarko_diurnal_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, dad)
         yarkoadd = (das+dad)*time
         yarkomin = min(yarkoadd, f180)
         if(present(gammin))then
            gammin = gam
         endif
      else
         yarkomin = f180
         if(present(gammin))then
            gammin = 180.0d0
         endif
      endif
   end subroutine computeYarkoMaxMin_xu2020

   subroutine computeYarko_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, time, yarko)
      real(kind=dkind), intent(in)  :: rho, K, C, radius, semiaxm, rotPer, gam, alpha, epsi, time
      real(kind=dkind), intent(out) :: yarko
      ! end interface
      real(kind=dkind) :: dad, das
      call yarko_seasonal_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, das)
      call yarko_diurnal_xu2020(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, dad)
      yarko = (das+dad)*time
   end subroutine computeYarko_xu2020

   ! Unified routines, using the linear interpolation to compute the solution 
   ! when the condition R'<<1 or R'>>1 is not satisfied
   subroutine yarko_seasonal_xu2020(rho, K, C, R, a0, gam, rotPer, alpha, epsi, yarko_s)
      real(kind=dkind), intent(out) :: yarko_s
      real(kind=dkind), intent(in) :: rho, K, C, R, a0, gam, rotPer, alpha, epsi
      ! end interface
      real(kind=dkind) :: mAst, mu
      real(kind=dkind) :: gam_rad, omega_rev
      real(kind=dkind) :: ls 
      real(kind=dkind) :: num, den
      real(kind=dkind) :: xs
      real(kind=dkind) :: c1s, c2s, c1d, c2d
      real(kind=dkind) :: Rps, Rsmall, Rbig
      real(kind=dkind) :: delta_a, delta_a_small, delta_a_big
      real(kind=dkind) :: aa, bb, ff
      real(kind=dkind) :: Rup, Rdown
!      Rup   = 2.5d0
!      Rdown = 0.7d0 
!      Rup   = 2.d0
      Rup   = 10.d0
!      Rdown = 0.1d0 
!      Rup   = 9.d0
      Rdown = 0.1d0 
      ! Convert gamma in radians
      gam_rad = gam*deg2rad
      ! Compute mass of the asteroid
      mAst = 4.d0*pi*rho*R**3.d0/3.d0
      ! Compute omega_rev 
      mu = gmsun+uGc*mAst
      omega_rev = sqrt(mu/(a0*au2m)**3.d0)
      ls = sqrt(K/(rho*C*omega_rev))
      ! Compute R'
      Rps = R/ls
      if(Rps .lt. Rdown)then
         ! R' << 1, use formula (20) to compute the Yarkovsky drift
         num = -(alpha*lumSun)**(1.75d0)*(epsi*sBoltz)**(0.25d0)*C*R**2.d0*sin(gam_rad)**2.d0
         den = 120.d0*pi**(1.75d0)*K**2.d0*clight*(a0*au2m)**(3.5d0)
         delta_a= num*m2au/den
      elseif(Rps .gt. Rup)then
         ! R' >> 1, use formula (22) to compute the Yarkovsky drift
         xs = sqrt(2.d0)*R*sqrt(rho*C*omega_rev/K)
         call compute_c1s_c1d(rho, K, C, R, a0, gam, rotPer, alpha, epsi, c1s, c1d)
         call compute_c2s_c2d(rho, K, C, R, a0, gam, rotPer, alpha, epsi, c2s, c2d)
         num = -4.d0*alpha*c1s*c2s*sin(gam_rad)**2.d0
         den = 9.d0*omega_rev*(2.d0+2.d0*c1s+c1s**2.d0)*xs
         delta_a= num*m2au/den
      else
         ! If it is in between, use a linear interpolation between the two cases
         Rsmall = Rdown*ls 
         Rbig   = Rup*ls
         ! Compute Yarkovsky drift for a small object, using formula (20)
         num = -(alpha*lumSun)**(1.75d0)*(epsi*sBoltz)**(0.25d0)*C*Rsmall**2.d0*sin(gam_rad)**2.d0
         den = 120.d0*pi**(1.75d0)*K**2.d0*clight*(a0*au2m)**(3.5d0)
         delta_a_small = num*m2au/den
         ! Compute Yarkovsky drift for a big object, using formula (22)
         xs = sqrt(2.d0)*Rbig*sqrt(rho*C*omega_rev/K)
         call compute_c1s_c1d(rho, K, C, Rbig, a0, gam, rotPer, alpha, epsi, c1s, c1d)
         call compute_c2s_c2d(rho, K, C, Rbig, a0, gam, rotPer, alpha, epsi, c2s, c2d)
         num = -4.d0*alpha*c1s*c2s*sin(gam_rad)**2.d0
         den = 9.d0*omega_rev*(2.d0+2.d0*c1s+c1s**2.d0)*xs
         delta_a_big = num*m2au/den
         ! Compute the Yarkovsky drift by linear interpolation, where we have used
         ! R_small' = 0.1, R_big' = 2
         aa = Rps-0.1d0
         bb = 2.d0-Rps
         ff = aa/(aa+bb)
         ! Linear interpolation
         delta_a = ff*delta_a_big+ (1.d0-ff)*delta_a_small
      endif
      yarko_s = delta_a
   end subroutine yarko_seasonal_xu2020

   subroutine yarko_diurnal_xu2020(rho, K, C, R, a0, gam, rotPer, alpha, epsi, yarko_d)
      real(kind=dkind), intent(out) :: yarko_d
      real(kind=dkind), intent(in) :: rho, K, C, R, a0, gam, rotPer, alpha, epsi
      ! end interface
      real(kind=dkind) :: mAst, mu
      real(kind=dkind) :: gam_rad, omega_rot, omega_rev
      real(kind=dkind) :: ls, ld
      real(kind=dkind) :: num, den
      real(kind=dkind) :: xd
      real(kind=dkind) :: c1s, c2s, c1d, c2d
      real(kind=dkind) :: Rpd, Rsmall, Rbig
      real(kind=dkind) :: delta_a, delta_a_small, delta_a_big
      real(kind=dkind) :: aa, bb, ff
      real(kind=dkind) :: Rup, Rdown
!      Rup   = 2.5d0
!      Rdown = 0.7d0 
!      Rup   = 2.d0
      Rup   = 10.d0
!      Rdown = 0.1d0 
!      Rup   = 9.d0
      Rdown = 0.1d0 
      ! Convert gamma in radians
      gam_rad = gam*deg2rad
      ! Compute mass of the asteroid
      mAst = 4.d0*pi*rho*R**3.d0/3.d0
      ! Compute omega_rev 
      mu = gmsun+uGc*mAst
      omega_rev = sqrt(mu/(a0*au2m)**3.d0)
      omega_rot = duepi/(rotPer*h2s)
      ls = sqrt(K/(rho*C*omega_rev))
      ld = ls*sqrt(omega_rev/omega_rot)
      ! Compute R'
      Rpd = R/ld
      if(Rpd .lt. Rdown)then
         ! R' << 1, use formula (20) to compute the Yarkovsky drift
         num = (alpha*lumSun)**(1.75d0)*(epsi*sBoltz)**(0.25d0)*C*R**2.d0*omega_rot*cos(gam_rad)
         den = 60.d0*pi**(1.75d0)*sqrt(gmSun)*K**2.d0*clight*(a0*au2m)**2.d0
         delta_a= num*m2au/den
      elseif(Rpd .gt. Rup)then
         ! R' >> 1, use formula (22) to compute the Yarkovsky drift
         xd = sqrt(2.d0)*R*sqrt(rho*C*omega_rot/K)
         call compute_c1s_c1d(rho, K, C, R, a0, gam, rotPer, alpha, epsi, c1s, c1d)
         call compute_c2s_c2d(rho, K, C, R, a0, gam, rotPer, alpha, epsi, c2s, c2d)
         num = 8.d0*alpha*c1d*c2d*cos(gam_rad)
         den = 9.d0*omega_rev*(2.d0+2.d0*c1d+c1d**2.d0)*xd
         delta_a= num*m2au/den
      else
         ! If it is in between, use a linear interpolation between the two cases
         Rsmall = Rdown*ld
         Rbig   = Rup*ld
         ! Compute Yarkovsky drift for a small object, using formula (20)
         num = (alpha*lumSun)**(1.75d0)*(epsi*sBoltz)**(0.25d0)*C*Rsmall**2.d0*omega_rot*cos(gam_rad)
         den = 60.d0*pi**(1.75d0)*sqrt(gmSun)*K**2.d0*clight*(a0*au2m)**2.d0
         delta_a_small = num*m2au/den
         ! Compute Yarkovsky drift for a big object, using formula (22)
         xd = sqrt(2.d0)*Rbig*sqrt(rho*C*omega_rot/K)
         call compute_c1s_c1d(rho, K, C, Rbig, a0, gam, rotPer, alpha, epsi, c1s, c1d)
         call compute_c2s_c2d(rho, K, C, Rbig, a0, gam, rotPer, alpha, epsi, c2s, c2d)
         num = 8.d0*alpha*c1d*c2d*cos(gam_rad)
         den =9.d0*omega_rev*(2.d0+2.d0*c1d+c1d**2.d0)*xd
         delta_a_big = num*m2au/den
         ! Compute the Yarkovsky drift by linear interpolation, where we have used
         ! R_small' = 0.1, R_big' = 2
         aa = Rpd-0.1d0
         bb = 2.d0-Rpd
         ff = aa/(aa+bb)
         ! Linear interpolation
         delta_a = ff*delta_a_big+ (1.d0-ff)*delta_a_small
      endif
      yarko_d = delta_a
   end subroutine yarko_diurnal_xu2020

   ! Auxiliary functions to compute the values of formula (16)
   subroutine compute_c1s_c1d(rho, K, C, R, a0, gam, rotPer, alpha, epsi, c1s, c1d)
      real(kind=dkind), intent(in)  :: rho, K, C, R, a0, gam, rotPer, alpha, epsi
      real(kind=dkind), intent(out) :: c1s, c1d
      ! end interface
      real(kind=dkind) :: c1s_num, c1s_den
   !   real(kind=dkind) :: c1d_num, c1d_den
      real(kind=dkind) :: gam_rad, omega_rot, omega_rev
      real(kind=dkind) :: mAst
      real(kind=dkind) :: mu
      gam_rad = gam*deg2rad
      ! Compute omega_rot and omega_rev
      mAst = 4.d0*pi*rho*R**3.d0/3.d0
      mu = gmsun+uGc*mAst
      omega_rev = sqrt(mu/(a0*au2m)**3.d0)
      omega_rot = duepi/(rotPer*h2s)
      ! Compute the numerator and denominator of c1s
      c1s_num = 2.d0*sqrt(2.d0)*sqrt(rho*K*C)*pi**(0.75d0)*mu**(0.25d0)*(a0*au2m)**(0.75d0)
      c1s_den = (alpha*lumSun)**(0.75d0)*(epsi*sBoltz)**(0.25d0)
      ! Compute c1s
      c1s = c1s_num/c1s_den
      ! Compute c1d using the formula with the periods
      c1d = c1s*sqrt(omega_rot/omega_rev)
   end subroutine 

   subroutine compute_c2s_c2d(rho, K, C, R, a0, gam, rotPer, alpha, epsi, c2s, c2d)
      real(kind=dkind), intent(in)  :: rho, K, C, R, a0, gam, rotPer, alpha, epsi
      real(kind=dkind), intent(out) :: c2s, c2d
      ! end interface
      real(kind=dkind) :: c2s_num, c2s_den
   !   real(kind=dkind) :: c2d_num, c2d_den
      real(kind=dkind) :: gam_rad, omega_rot, omega_rev
      real(kind=dkind) :: mAst
      real(kind=dkind) :: mu
      gam_rad = gam*deg2rad
      ! Compute omega_rot and omega_rev
      mAst = 4.d0*pi*rho*R**3.d0/3.d0
      mu = gmsun+uGc*mAst
      omega_rev = sqrt(mu/(a0*au2m)**3.d0)
      omega_rot = duepi/(rotPer*h2s)
      ! Compute the numerator and denominator of c2s
      c2s_num = 3.d0*sqrt(2.d0)*lumSun*sqrt(C)*mu**(0.25d0)
      c2s_den = 16.d0*pi*clight*sqrt(rho*K)*(a0*au2m)**(2.75d0)
      ! Compute c1s
      c2s = c2s_num/c2s_den
      ! Compute c1d using the formula with the periods
      c2d = c2s*sqrt(omega_rot/omega_rev)
   end subroutine 
   
   !============================================
   !====== FORMULAS BY VOKROUHLICKY 1999 =======
   !============================================
 
   subroutine computeYarkoMaxMin_vokrouhlicky(rho, K, C, radius, semiaxm, rotPer, alpha, epsi, yarkomax, yarkomin, gammin)
      real(kind=dkind), intent(in)  :: rho, K, C, radius, semiaxm, rotPer, alpha, epsi
      real(kind=dkind), intent(out) :: yarkomax, yarkomin
      real(kind=dkind), intent(out), optional :: gammin
      ! end interface
      real(kind=dkind) :: yarkoAdd
      real(kind=dkind) :: gam
      real(kind=dkind) :: f0, f90, f180
      real(kind=dkind) :: dad, das
      gam = 0.d0
      call computeYarko_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, yarkomax)

      f0 = yarkomax

      gam = 90.d0
      call computeYarko_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, f90)

      gam = 180.d0
      call computeYarko_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, f180)

      if(abs(f0/(2.d0*f90)).lt.1.d0)then
         gam = acos(f0/(2.d0*f90))*rad2deg 
         call computeYarko_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, yarkoadd)
         yarkomin = min(yarkoadd, f180)
         if(present(gammin))then
            gammin = gam
         endif
      else
         yarkomin = f180
         if(present(gammin))then
            gammin = 180.0d0
         endif
      endif
   end subroutine computeYarkoMaxMin_vokrouhlicky

   subroutine computeYarko_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, yarko)
      real(kind=dkind), intent(in)  :: rho, K, C, radius, semiaxm, rotPer, gam, alpha, epsi
      real(kind=dkind), intent(out) :: yarko
      ! end interface
      real(kind=dkind) :: dad, das
      call yarko_seasonal_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, das)
      call yarko_diurnal_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, dad)
      yarko = (das+dad)*m2au/s2my
   end subroutine computeYarko_vokrouhlicky

   subroutine yarko_seasonal_vokrouhlicky(rho, K, C, R, a0, gam, rotPer, alpha, epsi, yarko_s)
      real(kind=dkind), intent(out) :: yarko_s
      real(kind=dkind), intent(in) :: rho, K, C, R, a0, gam, rotPer, alpha, epsi
      ! end interface
      real(kind=dkind) :: mAst, mu
      real(kind=dkind) :: gam_rad, omega_rev
      real(kind=dkind) :: ls 
      real(kind=dkind) :: Theta, phi
      real(kind=dkind) :: E_star, T_star
      real(kind=dkind) :: dist
      real(kind=dkind) :: F_omega_rev
      real(kind=dkind) :: Rps
      ! Convert gamma in radians
      gam_rad = gam*deg2rad
      ! Compute mass of the asteroid
      mAst = 4.d0*pi*rho*R**3.d0/3.d0
      ! Compute omega_rev 
      mu = gmsun+uGc*mAst
      omega_rev = sqrt(mu/(a0*au2m)**3.d0)
      ls = sqrt(K/(rho*C*omega_rev))
      ! Compute R'
      Rps = R/ls
      ! Convert distance in meters
      dist   = a0*au2m
      ! Compute the solar radiation flux at distance dist
      E_star = lumSun/(4.d0*pi*dist**2.d0)
      ! Compute the radiation force
      phi    = pi*R**2.d0*E_star/(mAst*cLight)
      ! Compute the subsolar temperature
      T_star = (alpha*E_star/(epsi*sBoltz))**0.25d0
      ! Compute the thermal parameter
      Theta = sqrt(rho*K*C*omega_rev)/(epsi*sBoltz*T_star**3.d0)
      ! Compute F_omega_rev
      call Fnu_eval(Rps, Theta, F_omega_rev)
      ! Compute the Yarkovsky acceleration
      yarko_s = 4.d0*alpha*phi*F_omega_rev*sin(gam_rad)**2.d0/(9.d0*omega_rev)
   end subroutine yarko_seasonal_vokrouhlicky

   subroutine yarko_diurnal_vokrouhlicky(rho, K, C, R, a0, gam, rotPer, alpha, epsi, yarko_d)
      real(kind=dkind), intent(out) :: yarko_d
      real(kind=dkind), intent(in) :: rho, K, C, R, a0, gam, rotPer, alpha, epsi
      ! end interface
      real(kind=dkind) :: mAst, mu
      real(kind=dkind) :: gam_rad, omega_rot, omega_rev
      real(kind=dkind) :: ls, ld
      real(kind=dkind) :: phi
      real(kind=dkind) :: Rpd, Rsmall, Rbig
      real(kind=dkind) :: Theta
      real(kind=dkind) :: E_star, T_star
      real(kind=dkind) :: dist
      real(kind=dkind) :: F_omega_rot
      ! Convert gamma in radians
      gam_rad = gam*deg2rad
      ! Compute mass of the asteroid
      mAst = 4.d0*pi*rho*R**3.d0/3.d0
      ! Compute omega_rev 
      mu = gmsun+uGc*mAst
      omega_rev = sqrt(mu/(a0*au2m)**3.d0)
      omega_rot = duepi/(rotPer*h2s)
      ld = sqrt(K/(rho*C*omega_rot))
      ! Compute R'
      Rpd = R/ld
      ! Convert distance in meters
      dist   = a0*au2m
      ! Compute the solar radiation flux at distance dist
      E_star = lumSun/(4.d0*pi*dist**2.d0)
      ! Compute the radiation force
      phi    = pi*R**2.d0*E_star/(mAst*cLight)
      ! Compute the subsolar temperature
      T_star = (alpha*E_star/(epsi*sBoltz))**0.25d0
      ! Compute the thermal parameter
      Theta = sqrt(rho*K*C*omega_rot)/(epsi*sBoltz*T_star**3.d0)
      ! Compute F_omega_rev
      call Fnu_eval(Rpd, Theta, F_omega_rot)
      ! Compute the Yarkovsky acceleration
      yarko_d = -8.d0*alpha*phi*F_omega_rot*cos(gam_rad)/(9.d0*omega_rev)
   end subroutine yarko_diurnal_vokrouhlicky

   subroutine Fnu_eval(R, Theta, Fnu)
      real(kind=dkind), intent(in)  :: R, Theta
      real(kind=dkind), intent(out) :: Fnu
      ! end interface
      real(kind=dkind) :: k1, k2, k3
      real(kind=dkind) :: Rb, Rs
      real(kind=dkind) :: x, chi
      real(kind=dkind) :: A, B, C, D, U, V, den
      real(kind=dkind) :: xx1, xx2
      integer :: j
!      ! Constant values, valid only for R>>1
!      k1  = 0.5d0
!      k2  = 0.5d0
!      k3  = 0.5d0
!      Fnu = -k1*Theta/(1.d0 + 2.d0*k2*Theta + k3*Theta**2.d0)  
!      ! Formulation by Xu et. al. 2020
!!      x   = sqrt(2.d0)*R
!!      chi = Theta/x
!!      A = -(x+2.d0) - exp(x)*((x-2.d0)*cos(x) - x*sin(x))
!!      B = -x - exp(x)*(x*cos(x) + (x-2.d0)*sin(x))
!!      U = 3.d0*(x+2.d0) + exp(x)*(3.d0*(x-2.d0)*cos(x)+x*(x-3.d0)*sin(x))
!      V = x*(x+3.d0) - exp(x)*(x*(x-3.d0)*cos(x) - 3.d0*(x-2.d0)*sin(x))
!      C = A + chi*U/(1.d0 + chi)
!      D = B + chi*V/(1.d0 + chi)
!      Fnu = (B*C-A*D)/(C**2.d0 + D**2.d0)/(1.d0+chi)
     ! Formulation by Vokrouhlicky 1999 
      if(R.gt.30d0)then
         k1  = 0.5d0
         k2  = 0.5d0
         k3  = 0.5d0
      else
         x   = sqrt(2.d0)*R
         A = -(x+2.d0) - exp(x)*((x-2.d0)*cos(x) - x*sin(x))
         B = -x - exp(x)*(x*cos(x) + (x-2.d0)*sin(x))
         U = 3.d0*(x+2.d0) + exp(x)*(3.d0*(x-2.d0)*cos(x)+x*(x-3.d0)*sin(x))
         V = x*(x+3.d0) - exp(x)*(x*(x-3.d0)*cos(x) - 3.d0*(x-2.d0)*sin(x))
         den = x*(A**2.d0+B**2.d0);
         k1 = (A*V-B*U)/den
         k2 = (A*(A+U) + B*(B+V))/den
         k3 = ((A+U)**2.d0 + (B+V)**2.d0)/(den*x)
      endif
      Fnu = -k1*Theta/(1.d0 + 2.d0*k2*Theta + k3*Theta**2.d0)  
   end subroutine Fnu_eval

   subroutine k1k2k3_eval(R, Theta, k1, k2, k3)
      real(kind=dkind), intent(in)  :: R, Theta
      real(kind=dkind), intent(out) :: k1, k2, k3
      ! end interface
      real(kind=dkind) :: Rb, Rs
      real(kind=dkind) :: x, chi
      real(kind=dkind) :: A, B, C, D, U, V, den
      real(kind=dkind) :: xx1, xx2
      integer :: j
      if(R.gt.30d0)then
         k1  = 0.5d0
         k2  = 0.5d0
         k3  = 0.5d0
      else
         x   = sqrt(2.d0)*R
         A = -(x+2.d0) - exp(x)*((x-2.d0)*cos(x) - x*sin(x))
         B = -x - exp(x)*(x*cos(x) + (x-2.d0)*sin(x))
         U = 3.d0*(x+2.d0) + exp(x)*(3.d0*(x-2.d0)*cos(x)+x*(x-3.d0)*sin(x))
         V = x*(x+3.d0) - exp(x)*(x*(x-3.d0)*cos(x) - 3.d0*(x-2.d0)*sin(x))
         den = x*(A**2.d0+B**2.d0);
         k1 = (A*V-B*U)/den
         k2 = (A*(A+U) + B*(B+V))/den
         k3 = ((A+U)**2.d0 + (B+V)**2.d0)/(den*x)
      endif
   end subroutine k1k2k3_eval

   !===========================================
   !==== FOR YARKOVSKY ON ECCENTRIC ORBITS ====
   !===========================================

   ! Vector field due to Yarkovsky effect
   ! Ref:
   ! Vokrouhlický et. al. 2017
   ! Detailed Analysis of the Asteroid Pair (6070) Rheinland and (54827) 2001 NQ8
   subroutine yarkovsky_vf(a0, ecc, posvel, rho, K0, C, R, gam, rotPer, alpha, epsi, expo, yarko)
      real(kind=dkind), intent(in)  :: a0, ecc, posvel(6)
      real(kind=dkind), intent(in)  :: rho, K0, C 
      real(kind=dkind), intent(in)  :: R
      real(kind=dkind), intent(in)  :: gam, rotPer
      real(kind=dkind), intent(in)  :: alpha, epsi
      real(kind=dkind), intent(in)  :: expo
      real(kind=dkind), intent(out) :: yarko(3)
      ! end interface
      real(kind=dkind) :: K
      real(kind=dkind) :: pos(3), vel(3)
      real(kind=dkind) :: E_star, T_star
      real(kind=dkind) :: mAst, mu
      real(kind=dkind) :: omega_rot, omega_rev
      real(kind=dkind) :: gam_rad
      real(kind=dkind) :: kappa
      real(kind=dkind) :: Theta, gamma1, gamma2
      real(kind=dkind) :: Thetabar, gamma1bar, gamma2bar
      real(kind=dkind) :: ls, ld, Rps, Rpd, k1, k2, k3
      real(kind=dkind) :: verspos(3), dist
      real(kind=dkind) :: aux(3)
      real(kind=dkind) :: normalvec(3), normnormalvec
      real(kind=dkind) :: coeff1
      real(kind=dkind) :: e1(3), e2(3), e3(3)
      ! Take position and velocity
      pos(1:3) = posvel(1:3)
      vel(1:3) = posvel(4:6)
      ! Compute the mass of the asteroid
      mAst = 4.d0*pi*rho*R**3.d0/3.d0
      mu = gmsun+uGc*mAst
      ! Convert the obliquity in radians
      gam_rad = gam*deg2rad
      ! Compute the mean motion and the rotation frequency
      omega_rev = sqrt(mu/(a0*au2m)**3.d0)
      omega_rot = duepi/(rotPer*h2s)
      ! Compute penetration depths
      ls = sqrt(K/(rho*C*omega_rev))
      ld = ls*sqrt(omega_rev/omega_rot)
      ! e1 = s : unit vector of the spin axis
      ! e2 = n x s where n = pos/|pos|
      ! e3 = s x (n x s) = s x e2 = e1 x e2
      dist = norm2(pos)
      verspos = pos/dist
      K = K0*(dist*m2au)**(expo) 
      ! Compute e1, assume zero inclination for the asteroid
      e1(1) = sin(gam_rad)
      e1(2) = 0.d0
      e1(3) = cos(gam_rad)
      ! Compute e2
      call prvec(verspos, e1, e2)
      ! Compute e3
      call prvec(e1, e2, e3)
      ! Compute the solar flux and the subsolar temperature
      E_star = lumSun/(4.d0*pi*dist**2.d0)
      T_star = (alpha*E_star/(epsi*sBoltz))**0.25d0
      ! Compute the k coefficient
      kappa = 4.d0*alpha*pi*R**2.d0*E_star/(9.d0*mAst*cLight)
      ! Compute the coefficients of e2 and e3, representing the diurnal effect
      Rpd = R/ld
      Theta  = sqrt(rho*K*C*omega_rot)/(epsi*sBoltz*T_star**3.d0)
      call k1k2k3_eval(Rpd, Theta, k1, k2, k3) 
!      gamma1 = -0.5d0*Theta/(1.d0 + Theta + 0.5d0*Theta**2.d0)
!      gamma2 = (1.d0 + 0.5d0*Theta)/(1.d0 + Theta + 0.5d0*Theta**2.d0)
      gamma1 = -k1*Theta/(1.d0 + 2.d0*k2*Theta + k3*Theta**2.d0)
      gamma2 = (1.d0 + k2*Theta)/(1.d0 + 2.d0*k2*Theta + k3*Theta**2.d0)
      ! Compute the coefficient of e1, representing the seasonal effect
      Rps = R/ls
      Thetabar  = sqrt(rho*K*C*omega_rev)/(epsi*sBoltz*T_star**3.d0)
!      gamma1bar = -0.5d0*Thetabar/(1.d0 + Thetabar + 0.5d0*Thetabar**2.d0)
!      gamma2bar = (1.d0 + 0.5d0*Thetabar)/(1.d0 + Thetabar + 0.5d0*Thetabar**2.d0)
      call k1k2k3_eval(Rps, Thetabar, k1, k2, k3) 
      gamma1bar = -k1*Thetabar/(1.d0 + 2.d0*k2*Thetabar + k3*Thetabar**2.d0)
      gamma2bar = (1.d0 + k2*Thetabar)/(1.d0 + 2.d0*k2*Thetabar + k3*Thetabar**2.d0)
      ! Compute the normal vector to the orbital plane, N
      call prvec(pos, vel, normalvec)
      normnormalvec = norm2(normalvec)
      normalvec = normalvec/normnormalvec
      ! Compute N x n = N x pos/|pos|
      call prvec(normalvec, verspos, aux)  
      ! Compute the coefficient of e1
      ! \bar{gamma2} n \cdot s + \var{gamma1} (N x n) \cdot s
      coeff1 = gamma2bar*dot_product(verspos, e1) + gamma1bar*dot_product(aux, e1)

      ! Compute the total vector field of the Yarkovksy drift
      yarko = kappa*(coeff1*e1+gamma1*e2+gamma2*e3)
   end subroutine yarkovsky_vf

   subroutine yarko_eccentric(semiaxm, ecc, rho, K, C, R, gam, rotPer, alpha, epsi, expo, ecc_ye)
      real(kind=dkind), intent(in)  :: semiaxm, ecc, R
      real(kind=dkind), intent(in)  :: rho, K, C
      real(kind=dkind), intent(in)  :: gam, rotPer
      real(kind=dkind), intent(in)  :: alpha, epsi
      real(kind=dkind), intent(in)  :: expo
      real(kind=dkind), intent(out) :: ecc_ye
      ! end interface
      real(kind=dkind) :: meanMotion, mAst, mu, period
      real(kind=dkind) :: time, deltaTime
      real(kind=dkind) :: deltaEll
      real(kind=dkind) :: kep(6), car(6)
      real(kind=dkind) :: pos(3), vel(3), yarko(3)
      integer, parameter :: npoints = 500
      real(kind=dkind) :: dadt(npoints), yy
      real(kind=dkind) :: normvel
      integer :: j
      ! Compute the mass of the asteroid
      mAst = 4.d0*pi*rho*R**3.d0/3.d0
      ! Compute the mean motion of the asteroid
      mu = gmsun+uGc*mAst
      meanMotion = sqrt(mu/(semiaxm*au2m)**3.d0)
      ! Define the Keplerian elements of the orbits, in standard units
      kep(1) = semiaxm*au2m
      kep(2) = ecc
      kep(3) = 0.d0
      kep(4) = 0.d0
      kep(5) = 0.d0
      kep(6) = 0.d0
      ! Compute the period
      period = duepi
      deltaEll = period/float(npoints)
      ! Start evolution of the orbit, keeping fixed the elements and
      ! varying only the mean anomaly
      j = 0
      do while(kep(6).le.period)
         ! Compute the mean anomaly 
         kep(6) = float(j)*deltaEll 
         ! Convert from Keplerian elements to Cartesian elements
         ! Cartesian elements are in standard units
         call kep2car(kep, car)
         ! Take position and velocity
         pos(1:3) = car(1:3)
         vel(1:3) = car(4:6)
         normvel = norm2(vel)
         ! Compute the Yarkovsky vector field
         call yarkovsky_vf(semiaxm, ecc, car, rho, K, C, R, gam, rotPer, alpha, epsi, expo, yarko)
         ! Vary the force with the eccentricity of the orbit
         ! Here dadt is in standard units, km/s
         dadt(j) = 2.d0*dot_product(yarko, vel)/(meanMotion**2.d0*semiaxm*au2m)
         ! Convert the rate of change in AU/My and save the result
         dadt(j) = dadt(j)*m2au/s2my 
         j = j+1
      enddo
      ! Average the Yarkovsky force over the period
      call trapezoid_average(period, deltaEll, dadt, npoints, ecc_ye)
   end subroutine yarko_eccentric 

!   subroutine yarko_eccentric_old(semiaxm, ecc, rho, K, C, R, gam, rotPer, alpha, epsi, ecc_ye)
!      real(kind=dkind), intent(in)  :: semiaxm, ecc, R
!      real(kind=dkind), intent(in)  :: rho, K, C
!      real(kind=dkind), intent(in)  :: gam, rotPer
!      real(kind=dkind), intent(in)  :: alpha, epsi
!      real(kind=dkind), intent(out) :: ecc_ye
!      ! end interface
!      real(kind=dkind) :: meanMotion, mAst, mu, period
!      real(kind=dkind) :: time, deltaTime
!      real(kind=dkind) :: deltaEll
!      real(kind=dkind) :: kep(6), car(6)
!      real(kind=dkind) :: pos(3), vel(3), yarko(3)
!      integer, parameter :: npoints = 4000
!      real(kind=dkind) :: dadt(npoints), yy
!      real(kind=dkind) :: normvel
!      integer :: j
!      ! Compute the mass of the asteroid
!      mAst = 4.d0*pi*rho*R**3.d0/3.d0
!      ! Compute the mean motion of the asteroid
!      mu = gmsun+uGc*mAst
!      meanMotion = sqrt(mu/(semiaxm*au2m)**3.d0)
!      ! Define the Keplerian elements of the orbits, in standard units
!      kep(1) = semiaxm*au2m
!      kep(2) = ecc
!      kep(3) = 0.d0
!      kep(4) = 0.d0
!      kep(5) = 0.d0
!      kep(6) = 0.d0
!      ! Compute the period
!      period = duepi/meanMotion
!      deltaTime = period/float(npoints)
!      ! Start evolution of the orbit, keeping fixed the elements and
!      ! varying only the mean anomaly
!      time = 0.d0
!      j = 0
!      do while(time.le.period)
!         ! Compute the mean anomaly 
!         kep(6) = meanMotion*time
!         ! Convert from Keplerian elements to Cartesian elements
!         ! Cartesian elements are in standard units
!         call kep2car(kep, car)
!         ! Take position and velocity
!         pos(1:3) = car(1:3)
!         vel(1:3) = car(4:6)
!         normvel = norm2(vel)
!         ! Compute the Yarkovsky vector field
!         call yarkovsky_vf(semiaxm, ecc, car, rho, K, C, R, gam, rotPer, alpha, epsi, expo, yarko)
!         ! Vary the force with the eccentricity of the orbit
!         ! Here dadt is in standard units, km/s
!         dadt(j) = 2.d0*dot_product(yarko, vel)/(meanMotion**2.d0*semiaxm*au2m)
!         ! Convert the rate of change in AU/My and save the result
!         dadt(j) = dadt(j)*m2au/s2my 
!         time = time + deltaTime
!         j = j+1
!      enddo
!      ! Average the Yarkovsky force over the period
!      call trapezoid_average(period*s2my, deltaTime*s2my, dadt, npoints, ecc_ye)
!   end subroutine yarko_eccentric_old

   !=====================================================
   !============== CHANGE OF COORDINATES ================
   !=====================================================

   subroutine keplereq(ell,ecc,u)
      real(kind=dkind), intent(in) :: ell,ecc
      real(kind=dkind), intent(out) :: u
      ! end interface
      integer :: j
      integer,parameter :: jmax=100
      real(kind=dkind) :: ukp1,uk, elle
      uk=pi
      elle = mod(ell,2.d0*pi)
      do j=1,jmax
         ukp1 = uk + (elle-uk+ecc*sin(uk))/(1.d0-ecc*cos(uk))
         if(abs(ukp1-uk).le.1.d-12)then
            exit
         endif
         uk=ukp1
      enddo
      u = ukp1
   end subroutine keplereq

   subroutine kep2car(kep, car)
      real(kind=dkind),intent(in)  :: kep(6) ! Y are Delaunay elements
      real(kind=dkind),intent(out) :: car(6)
      !end interface
      real(kind=dkind) :: kappa, kappa2
      real(kind=dkind) :: inc, ecc, aa
      real(kind=dkind) :: L,G,Z,ell,omeg,Omnod
      real(kind=dkind) :: ci,si,co,so,con,son, cl, sl
      real(kind=dkind) :: didG,didZ
      real(kind=dkind),dimension(3,3) :: R_om,R_i,R_Omn
      real(kind=dkind),dimension(3,3) :: RiRom,ROmnRi
      real(kind=dkind),dimension(3,3) :: dR_om,dR_i,dR_Omn !,dRomRi,RomdRi
      real(kind=dkind),dimension(3,3) :: Rot,dRotdOmn,dRotdi,dRotdi_tmp,dRotdom
      real(kind=dkind) :: pos(3,1)
      real(kind=dkind) :: dposdell(3,1)
      real(kind=dkind) :: dposdL(3,1),dposdG(3,1),dposdZ(3,1)
      real(kind=dkind) :: u,sinu,cosu,dcosudL,dcosudG,dudell,dudL,dudG,dudecc
      real(kind=dkind) :: deccdL,deccdG
      real(kind=dkind) :: dposdu(3,1),vel(3,1),versvel(3,1),normvel
      real(kind=dkind) :: tmpcar(6,1),norma,normpos
      kappa = sqrt(gmSun)!*au2m**(-1.5d0)*d2s
      kappa2 = kappa**2.d0
      aa    = kep(1)
      ecc   = kep(2)
      inc   = kep(3)*deg2rad
      omeg  = kep(4)
      Omnod = kep(5)
      ell   = kep(6)

      L = kappa*sqrt(aa)
      G = L*sqrt(1.d0-ecc**2)
      Z = G*cos(inc)
       
      ! Rotation matrix of angle omega for the asteroid
      co=cos(omeg)
      so=sin(omeg)
      R_om = 0.d0
      R_om(1,1) = co
      R_om(1,2) =-so
      R_om(2,1) = so
      R_om(2,2) = co
      R_om(3,3) = 1.d0

      ! Rotation matrix of angle Inclination for the asteroid
      ci = cos(inc)
      si = sin(inc)
      R_i = 0.d0
      R_i(1,1) = 1.d0
      R_i(2,2) = ci
      R_i(2,3) =-si
      R_i(3,2) = si
      R_i(3,3) = ci
      
      ! R_i*R_om
      RiRom = MATMUL(R_i,R_om)
      
      ! Rotation matrix of angle omega for the asteroid
      con=cos(Omnod)
      son=sin(Omnod)
      R_Omn = 0.d0
      R_Omn(1,1) = con
      R_Omn(1,2) =-son
      R_Omn(2,1) = son
      R_Omn(2,2) = con
      R_Omn(3,3) = 1.d0
    
      ROmnRi = MATMUL(R_Omn,R_i)
      ! R_Omn*R_i*R_om
      Rot = MATMUL(R_Omn,RiRom)
          
      !  write(*,*)'a,ecc=',aa,ecc
      CALL keplereq(ell,ecc,u)
      !  write(*,*)'u=',u,'u-e*sinu=',u-ecc*sin(u)
      cosu = cos(u); sinu = sin(u)
       
      pos(1,1) = aa*(cosu - ecc)
      pos(2,1) = aa*sqrt(1.d0-ecc**2)*sinu 
      pos(3,1) = 0.d0 

      dposdu(1,1) = L**2/kappa2*(-sinu)
      dposdu(2,1) = L*G*cosu/kappa2
      dposdu(3,1) = 0.d0
      
      norma=sqrt(DOT_PRODUCT(dposdu(1:3,1),dposdu(1:3,1)))
      versvel(1:3,1)=dposdu(1:3,1)/norma
      normpos=sqrt(DOT_PRODUCT(pos(1:3,1),pos(1:3,1)))
      normvel =kappa*sqrt(2.d0/normpos-1.d0/aa) 
      vel(1:3,1)=normvel*versvel(1:3,1)
      
      tmpcar(1:3,1) =MATMUL(Rot,pos(1:3,1))
      tmpcar(4:6,1) =MATMUL(Rot,vel(1:3,1))    
      car(1:6) = tmpcar(1:6,1)
   end subroutine kep2car

   subroutine prvec(a,b,c) 
      use used_const
      implicit none
      real(kind=dkind), intent(in)  :: a(3), b(3)
      real(kind=dkind), intent(out) :: c(3)
      c(1)=a(2)*b(3)-a(3)*b(2) 
      c(2)=a(3)*b(1)-a(1)*b(3) 
      c(3)=a(1)*b(2)-a(2)*b(1) 
   end subroutine prvec

   !=================================
   !========== AVERAGING ============
   !=================================

   ! Implementation of the trapezoid rule for the computation of integrals
   subroutine trapezoid_average(time, deltaT, fun, npoints, average)
      use used_const
      implicit none
      integer, intent(in)           :: npoints
      real(kind=dkind), intent(in)  :: time, deltaT
      real(kind=dkind), intent(in)  :: fun(npoints)
      real(kind=dkind), intent(out) :: average
      ! end interface
      real(kind=dkind) :: integral
      integer :: j
     
      integral = 0.d0
      do j=2, npoints-1
         integral = integral + fun(j)
      enddo
      integral = integral + 0.5d0*(fun(1)+fun(npoints))
      integral = integral*deltaT

      average = integral/time
   end subroutine trapezoid_average

end module yarko_force
