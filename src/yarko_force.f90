!      Authors: Marco Fenucci, Bojan Novakovic, Dusan Marceta
!  Institution: University of Belgrade
!         Date: July 2022
! 
! PURPOSE: This module contains the implementations of the Yarkovsky models
module yarko_force
   use used_const
   implicit none
   ! Private subroutines
   private ::  Fnu_eval, cross_prod, kep2car, keplereq, trapezoid_average
   ! Public subroutines
   public :: computeYarko_vokrouhlicky,   computeYarkoMaxMin_vokrouhlicky, &
           & yarko_seasonal_vokrouhlicky, yarko_diurnal_vokrouhlicky, &
           & yarko_eccentric,             yarkovsky_vf 
   contains
   
   !============================================
   !====== FORMULAS BY VOKROUHLICKY 1999 =======
   !============================================
 
   ! PURPOSE:
   !         
   !
   ! INPUT:
   !
   ! OUTPUT:
   subroutine computeYarkoMaxMin_vokrouhlicky(rho, K, C, radius, semiaxm, rotPer, alpha, epsi, yarkomax, yarkomin, gammin)
      real(kind=dkind), intent(in)            :: rho, K, C, radius, semiaxm, rotPer, alpha, epsi
      real(kind=dkind), intent(out)           :: yarkomax, yarkomin
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

   ! PURPOSE: Compute the Yarkovsky drift with the analytical model by Vokrouhlicky 1998,
   !          that assumes a circular orbit for the asteroid
   !
   ! INPUT:
   !
   ! OUTPUT:
   subroutine computeYarko_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, yarko)
      real(kind=dkind), intent(in)  :: rho, K, C, radius, semiaxm, rotPer, gam, alpha, epsi
      real(kind=dkind), intent(out) :: yarko
      ! end interface
      real(kind=dkind) :: dad, das
      ! Compute the seasonal component
      call yarko_seasonal_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, das)
      ! Compute the diurnal component
      call yarko_diurnal_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, dad)
      ! Add the two contributions and convert in au/My
      yarko = (das+dad)*m2au/s2my
   end subroutine computeYarko_vokrouhlicky

   ! PURPOSE: Computes the seasonal component of the Yarkovsky drift
   !
   ! INPUT:
   !
   ! OUTPUT:
   subroutine yarko_seasonal_vokrouhlicky(rho, K, C, R, a0, gam, rotPer, alpha, epsi, yarko_s)
      real(kind=dkind), intent(in)  :: rho, K, C, R, a0, gam, rotPer, alpha, epsi
      real(kind=dkind), intent(out) :: yarko_s
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

   ! PURPOSE: Computes the diurnal component of the Yarkovsky drift
   !
   ! INPUT:
   !
   ! OUTPUT:
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
      omega_rot = twopi/(rotPer*h2s)
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
      Theta  = sqrt(rho*K*C*omega_rot)/(epsi*sBoltz*T_star**3.d0)
      ! Compute F_omega_rev
      call Fnu_eval(Rpd, Theta, F_omega_rot)
      ! Compute the Yarkovsky acceleration
      yarko_d = -8.d0*alpha*phi*F_omega_rot*cos(gam_rad)/(9.d0*omega_rev)
   end subroutine yarko_diurnal_vokrouhlicky

   ! PURPOSE: Evaluates the function of Eq. (A.7) in Fenucci et al. 2021
   !
   ! INPUT: 
   !     R : radius of the asteroid [m]
   ! Theta : thermal parameter
   !
   ! OUTPUT:
   !   Fnu : value of the function 
   subroutine Fnu_eval(R, Theta, Fnu)
      real(kind=dkind), intent(in)  :: R, Theta
      real(kind=dkind), intent(out) :: Fnu
      ! end interface
      real(kind=dkind) :: k1, k2, k3
      real(kind=dkind) :: Rb, Rs
      real(kind=dkind) :: x, chi
      real(kind=dkind) :: A, B, C, D, U, V, den
      real(kind=dkind) :: xx1, xx2
      ! Formulation by Vokrouhlicky 1999 
      if(R.gt.30d0)then
         ! Constant values, valid only for R>>1
         k1  = 0.5d0
         k2  = 0.5d0
         k3  = 0.5d0
      else
         x   = sqrt(2.d0)*R
         A   = -(x+2.d0) - exp(x)*((x-2.d0)*cos(x) - x*sin(x))
         B   = -x - exp(x)*(x*cos(x) + (x-2.d0)*sin(x))
         U   = 3.d0*(x+2.d0) + exp(x)*(3.d0*(x-2.d0)*cos(x) + x*(x-3.d0)*sin(x))
         V   = x*(x+3.d0) - exp(x)*(x*(x-3.d0)*cos(x) - 3.d0*(x-2.d0)*sin(x))
         den = x*(A**2.d0 + B**2.d0)
         k1  = (A*V - B*U)/den
         k2  = (A*(A + U) + B*(B + V))/den
         k3  = ((A + U)**2.d0 + (B + V)**2.d0)/(den*x)
      endif
      Fnu = -k1*Theta/(1.d0 + 2.d0*k2*Theta + k3*Theta**2.d0)  
   end subroutine Fnu_eval

   ! PURPOSE: Computes the coefficients k1 k2 k3 for the evaluation of Eq. (A.7) of Fenucci et al. 2021
   !
   ! INPUT:
   !          R : radius of the asteroid [m]
   !      Theta : thermal parameter
   !
   ! OUTPUT:
   ! k1, k2, k3 : values of the coefficients
   subroutine k1k2k3_eval(R, Theta, k1, k2, k3)
      real(kind=dkind), intent(in)  :: R, Theta
      real(kind=dkind), intent(out) :: k1, k2, k3
      ! end interface
      real(kind=dkind) :: Rb, Rs
      real(kind=dkind) :: x, chi
      real(kind=dkind) :: A, B, C, D, U, V, den
      real(kind=dkind) :: xx1, xx2
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

   ! PURPOSE: Compute the vector field due to Yarkovsky effect on a given position along the orbit
   !          Ref: Vokrouhlický et. al. 2017, Detailed Analysis of the Asteroid Pair (6070) Rheinland and (54827) 2001 NQ8
   !
   ! INPUT: 
   !    kep : Keplerian elements of the asteroid
   ! posvel : Cartesian elements (in SI units)
   !    rho : density of the asteroid
   !     K0 : thermal conductivity
   !      C : heat capacity
   !      R : radius
   !    gam : obliquity
   ! rotPer : rotation period
   !  alpha : absorption coefficient
   !   epsi : emissivity
   !   expo : exponent for the variation of K along the orbit
   !
   ! OUTPUT:
   !  yarko : Yarkovsky vector field
   subroutine yarkovsky_vf(kep, posvel, rho, K0, C, R, gam, rotPer, alpha, epsi, expo, yarko)
      real(kind=dkind), intent(in)  :: kep(6), posvel(6)
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
      ! These angles are in radians
      real(kind=dkind) :: inc, omega, OmNod
      real(kind=dkind) :: co, so, ci, si, con, son
      real(kind=dkind) :: R_om(3,3), R_i(3,3), R_Omn(3,3), RiRom(3,3), Rot(3,3)
      ! Unit vector of the obliquity referred to the orbital plane
      real(kind=dkind) :: e1_op(3)
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
      mu   = gmsun + uGc*mAst
      ! Convert the obliquity in radians
      gam_rad   = gam*deg2rad
      ! Compute the mean motion and the rotation frequency
      omega_rev = sqrt(mu/(kep(1)*au2m)**3.d0)
      omega_rot = twopi/(rotPer*h2s)
      ! Compute penetration depths
      ls = sqrt(K/(rho*C*omega_rev))
      ld = ls*sqrt(omega_rev/omega_rot)
      ! Compute the unit vectors
      ! e1 = s : unit vector of the spin axis
      ! e2 = n x s where n = pos/|pos|
      ! e3 = s x (n x s) = s x e2 = e1 x e2
      dist    = norm2(pos)
      verspos = pos/dist
      ! Compute the thermal conductivity as a function
      ! of the distance from the Sun. This formulation
      ! assumes that rho and C does not change with the 
      ! distance
      K       = K0*(dist*m2au)**(expo) 
      ! Start computing e1 
      ! Take the angles defining the orbital plane
      inc    = kep(3)*deg2rad
      Omnod  = kep(4)*deg2rad
      omega  = kep(5)*deg2rad
      ! Rotation matrix of angle omega 
      co   = cos(omega)
      so   = sin(omega)
      R_om = 0.d0
      R_om(1,1) = co
      R_om(1,2) =-so
      R_om(2,1) = so
      R_om(2,2) = co
      R_om(3,3) = 1.d0
      ! Rotation matrix of angle inclination
      ci  = cos(inc)
      si  = sin(inc)
      R_i = 0.d0
      R_i(1,1) = 1.d0
      R_i(2,2) = ci
      R_i(2,3) =-si
      R_i(3,2) = si
      R_i(3,3) = ci
      ! R_i*R_om
      RiRom = MATMUL(R_i,R_om)
      ! Rotation matrix of angle Omega for the asteroid
      con   = cos(Omnod)
      son   = sin(Omnod)
      R_Omn = 0.d0
      R_Omn(1,1) = con
      R_Omn(1,2) =-son
      R_Omn(2,1) = son
      R_Omn(2,2) = con
      R_Omn(3,3) = 1.d0
      ! R_Omn*R_i*R_om: Rotation matrix for the transformation of a vector
      ! from the orbital plane to the 3d space
      Rot = MATMUL(R_Omn,RiRom)
      ! Unit vector of the obliquity, with gamma defined wrt the orbital plane
      e1_op(1) = sin(gam_rad)
      e1_op(2) = 0.d0
      e1_op(3) = cos(gam_rad)
      ! Computation of e1 by rotation
      e1 = matmul(Rot, e1_op)
      ! Compute e2
      call cross_prod(verspos, e1, e2)
      ! Compute e3
      call cross_prod(e1, e2, e3)
      ! Compute the solar flux and the subsolar temperature
      E_star = lumSun/(4.d0*pi*dist**2.d0)
      T_star = (alpha*E_star/(epsi*sBoltz))**0.25d0
      ! Compute the k coefficient
      kappa = 4.d0*alpha*pi*R**2.d0*E_star/(9.d0*mAst*cLight)
      ! Compute the coefficients of e2 and e3, representing the diurnal effect
      Rpd    = R/ld
      Theta  = sqrt(rho*K*C*omega_rot)/(epsi*sBoltz*T_star**3.d0)
      call k1k2k3_eval(Rpd, Theta, k1, k2, k3) 
      gamma1 = -k1*Theta/(1.d0 + 2.d0*k2*Theta + k3*Theta**2.d0)
      gamma2 = (1.d0 + k2*Theta)/(1.d0 + 2.d0*k2*Theta + k3*Theta**2.d0)
      ! Compute the coefficient of e1, representing the seasonal effect
      Rps       = R/ls
      Thetabar  = sqrt(rho*K*C*omega_rev)/(epsi*sBoltz*T_star**3.d0)
      call k1k2k3_eval(Rps, Thetabar, k1, k2, k3) 
      gamma1bar = -k1*Thetabar/(1.d0 + 2.d0*k2*Thetabar + k3*Thetabar**2.d0)
      gamma2bar = (1.d0 + k2*Thetabar)/(1.d0 + 2.d0*k2*Thetabar + k3*Thetabar**2.d0)
      ! Compute the normal vector to the orbital plane, N
      call cross_prod(pos, vel, normalvec)
      normnormalvec = norm2(normalvec)
      normalvec = normalvec/normnormalvec
      ! Compute N x n = N x pos/|pos|
      call cross_prod(normalvec, verspos, aux)  
      ! Compute the coefficient of e1
      ! \bar{gamma2} n \cdot s + \var{gamma1} (N x n) \cdot s
      coeff1 = gamma2bar*dot_product(verspos, e1) + gamma1bar*dot_product(aux, e1)
      ! Compute the total vector field of the Yarkovksy drift
      yarko = kappa*(coeff1*e1 + gamma1*e2 + gamma2*e3)
   end subroutine yarkovsky_vf

   ! PURPOSE:
   !         
   !
   ! INPUT:
   !
   ! OUTPUT:
   ! TODO: this is deprecated, it can be deleted
   subroutine yarko_eccentric_ell(kep, rho, K, C, R, gam, rotPer, alpha, epsi, expo, ecc_ye)
      real(kind=dkind), intent(in)  :: kep(6), R
      real(kind=dkind), intent(in)  :: rho, K, C
      real(kind=dkind), intent(in)  :: gam, rotPer
      real(kind=dkind), intent(in)  :: alpha, epsi
      real(kind=dkind), intent(in)  :: expo
      real(kind=dkind), intent(out) :: ecc_ye
      ! end interface
      real(kind=dkind)   :: meanMotion, mAst, mu, period
      real(kind=dkind)   :: time, deltaTime
      real(kind=dkind)   :: deltaEll
      real(kind=dkind)   :: car(6), kep_c(6)
      real(kind=dkind)   :: pos(3), vel(3), yarko(3)
      integer, parameter :: npoints = 500
      real(kind=dkind)   :: dadt(npoints+1)
      integer :: j
      ! Compute the mass of the asteroid
      mAst = 4.d0*pi*rho*R**3.d0/3.d0
      ! Compute the mean motion of the asteroid
      mu = gmsun+uGc*mAst
      meanMotion = sqrt(mu/(kep(1)*au2m)**3.d0)
      ! Compute the period
      period   = twopi 
      deltaEll = period/float(npoints)
      ! Start evolution of the orbit, keeping fixed the elements and
      ! varying only the mean anomaly
      j = 1
      ! TODO: Is everything symmetric wrt half period? If so, we can integrate only for half of the period,
      !       and it should be faster.
      ! TODO: I think it is better to integrate in eccentric anomaly! The relation between the two is
      !        ell = u - e sin u
      !       and use constant u. Maybe we can also change the coordinates:
      !       \int f(l) dl = \int f(l(u)) dl/du du = \int f(l(u)) (1-e cos(u)) du
      kep_c    = kep 
      kep_c(6) = 0.d0
      do while(kep_c(6).lt.period)
         ! Compute the mean anomaly
         kep_c(6) = float(j-1)*deltaEll
         ! Convert from Keplerian elements to Cartesian elements
         ! Cartesian elements are in standard units
         call kep2car(kep_c, car)
         ! Take position and velocity
         pos(1:3) = car(1:3)
         vel(1:3) = car(4:6)
         ! Compute the Yarkovsky vector field
         call yarkovsky_vf(kep_c, car, rho, K, C, R, gam, rotPer, alpha, epsi, expo, yarko)
         ! Vary the force with the eccentricity of the orbit
         ! Here dadt is in standard units, km/s
         dadt(j) = 2.d0*dot_product(yarko, vel)/(meanMotion**2.d0*kep(1)*au2m)
         ! Convert the rate of change in AU/My and save the result
         dadt(j) = dadt(j)*m2au/s2my
         j = j+1
      enddo
      ! Average the Yarkovsky force over the period
      call trapezoid_average(period, deltaEll, dadt, npoints+1, ecc_ye)
   end subroutine yarko_eccentric_ell


   ! PURPOSE: Compute the average of the Yarkovsky drift over an orbital period, by integrating the Yarkovsky vector field 
   !          along the orbit. Note that the integration is done in eccentric anomaly, because it provides a better grid 
   !          near the perihelion.
   !
   ! INPUT:
   !    kep : Keplerian elements of the asteroid
   !    rho : density of the asteroid
   !      K : thermal conductivity
   !      C : heat capacity
   !      R : radius
   !    gam : obliquity
   ! rotPer : rotation period
   !  alpha : absorption coefficient
   !   epsi : emissivity
   !   expo : exponent for the variation of K along the orbit
   !
   ! OUTPUT:
   ! ecc_ye : averaged Yarkovsky effect over an orbital period
   subroutine yarko_eccentric(kep, rho, K, C, R, gam, rotPer, alpha, epsi, expo, ecc_ye)
      real(kind=dkind), intent(in)  :: kep(6), R
      real(kind=dkind), intent(in)  :: rho, K, C
      real(kind=dkind), intent(in)  :: gam, rotPer
      real(kind=dkind), intent(in)  :: alpha, epsi
      real(kind=dkind), intent(in)  :: expo
      real(kind=dkind), intent(out) :: ecc_ye
      ! end interface
      real(kind=dkind)   :: meanMotion, mAst, mu, period
      real(kind=dkind)   :: time, deltaTime
      real(kind=dkind)   :: u, ell, deltau, dldu
      real(kind=dkind)   :: car(6), kep_c(6)
      real(kind=dkind)   :: pos(3), vel(3), yarko(3)
      integer, parameter :: npoints = 500
      real(kind=dkind)   :: dadt(npoints+1)
      integer            :: j
      ! Compute the mass of the asteroid
      mAst = 4.d0*pi*rho*R**3/3.d0
      ! Compute the mean motion of the asteroid
      mu = gmsun + uGc*mAst
      meanMotion = sqrt(mu/(kep(1)*au2m)**3)
      ! The period is always 2pi, compute the stepsize 
      period = twopi 
      deltau = period/float(npoints)
      u      = 0.d0
      ! Start evolution of the orbit, keeping fixed the elements and
      ! varying only the eccentric anomaly
      j = 1
      ! Here we integrate in the eccentric anomaly u. The relation with the mean anomaly ell is 
      !        ell = u - e sin u
      ! The integral is done by changing coordinates:
      !       \int f(l) dl = \int f(l(u)) dl/du du = \int f(l(u)) (1-e cos(u)) du
      kep_c = kep
      do while(u.lt.period)
         u = float(j-1)*deltau
         ! Compute the mean anomaly
         kep_c(6) = (u - kep(2)*sin(u))*rad2deg
         ! Convert from Keplerian elements to Cartesian elements
         ! Cartesian elements are in standard units
         call kep2car(kep_c, car)
         ! Take position and velocity
         pos(1:3) = car(1:3)
         vel(1:3) = car(4:6)
         ! Compute the Yarkovsky vector field
         call yarkovsky_vf(kep_c, car, rho, K, C, R, gam, rotPer, alpha, epsi, expo, yarko)
         ! Vary the force with the eccentricity of the orbit
         ! Here dadt is in standard units, km/s
         dldu    = 1.d0 - kep(2)*cos(u)
         dadt(j) = 2.d0*dot_product(yarko, vel)/(meanMotion**2*kep(1)*au2m)*dldu
         ! Convert the rate of change in AU/My and save the result
         dadt(j) = dadt(j)*m2au/s2my
         j = j+1
      enddo
      ! Average the Yarkovsky force over the period
      call trapezoid_average(period, deltau, dadt, npoints+1, ecc_ye)
   end subroutine yarko_eccentric


   !=====================================================
   !============== CHANGE OF COORDINATES ================
   !=====================================================

   ! PURPOSE: Solve the Kepler equation ell = u - e sin u wrt u, using the Newton method
   !
   ! INPUT:
   !    ell : mean anomaly [rad]
   !    ecc : eccentricity [rad]
   !
   ! OUTPUT:
   !      u : eccentric anomaly
   subroutine keplereq(ell, ecc, u)
      real(kind=dkind), intent(in)  :: ell, ecc
      real(kind=dkind), intent(out) :: u
      ! end interface
      real(kind=dkind)  :: ukp1,uk, elle
      integer           :: j
      integer,parameter :: jmax=100
      ! Set the first guess to pi
      uk   = pi
      ! Start the Newton method
      elle = mod(ell, twopi)
      do j=1,jmax
         ukp1 = uk + (elle-uk+ecc*sin(uk))/(1.d0-ecc*cos(uk))
         ! When the difference is smaller than eps = 1.d-12, stop the iterations
         if(abs(ukp1-uk).le.1.d-12)then
            exit
         endif
         uk=ukp1
      enddo
      u = ukp1
   end subroutine keplereq

   ! PURPOSE: Converts Keplerian elements to Cartesian coordinates
   !         
   ! INPUT:
   !    kep : Keplerian elements, semimajor axis in au, angles in degrees
   !
   ! OUTPUT:
   !    car : Cartesian coordinates, in m for the positions and m/s for the velocities
   subroutine kep2car(kep, car)
      real(kind=dkind),intent(in)  :: kep(6) 
      real(kind=dkind),intent(out) :: car(6)
      !end interface
      real(kind=dkind) :: kappa, kappa2
      real(kind=dkind) :: inc, ecc, aa
      real(kind=dkind) :: L,G,Z,ell,omeg,Omnod
      real(kind=dkind) :: ci,si,co,so,con,son, cl, sl
      real(kind=dkind) :: didG,didZ
      real(kind=dkind),dimension(3,3) :: R_om, R_i,R_Omn
      real(kind=dkind),dimension(3,3) :: RiRom, ROmnRi
      real(kind=dkind),dimension(3,3) :: Rot
      real(kind=dkind) :: pos(3,1)
      real(kind=dkind) :: u,sinu,cosu
      real(kind=dkind) :: dposdu(3,1),vel(3,1),versvel(3,1),normvel
      real(kind=dkind) :: tmpcar(6,1),norma,normpos
      kappa  = sqrt(gmSun)
      kappa2 = kappa**2.d0
      aa     = kep(1)*au2m
      ecc    = kep(2)
      inc    = kep(3)*deg2rad
      omeg   = kep(4)*deg2rad
      Omnod  = kep(5)*deg2rad
      ell    = kep(6)*deg2rad
      ! Delaunay Elements
      L = kappa*sqrt(aa)
      G = L*sqrt(1.d0-ecc**2)
      Z = G*cos(inc)
      ! Rotation matrix of angle omega for the asteroid
      co   = cos(omeg)
      so   = sin(omeg)
      R_om = 0.d0
      R_om(1,1) = co
      R_om(1,2) =-so
      R_om(2,1) = so
      R_om(2,2) = co
      R_om(3,3) = 1.d0
      ! Rotation matrix of angle Inclination for the asteroid
      ci  = cos(inc)
      si  = sin(inc)
      R_i = 0.d0
      R_i(1,1) = 1.d0
      R_i(2,2) = ci
      R_i(2,3) =-si
      R_i(3,2) = si
      R_i(3,3) = ci
      ! R_i*R_om
      RiRom = matmul(R_i,R_om)
      ! Rotation matrix of angle omega for the asteroid
      con   = cos(Omnod)
      son   = sin(Omnod)
      R_Omn = 0.d0
      R_Omn(1,1) = con
      R_Omn(1,2) =-son
      R_Omn(2,1) = son
      R_Omn(2,2) = con
      R_Omn(3,3) = 1.d0
      ROmnRi = matmul(R_Omn,R_i)
      ! R_Omn*R_i*R_om
      Rot = matmul(R_Omn,RiRom)
      ! Solve the Kepler equation
      call keplereq(ell,ecc,u)
      cosu = cos(u)
      sinu = sin(u)
      ! Compute the position
      pos(1,1) = aa*(cosu - ecc)
      pos(2,1) = aa*sqrt(1.d0-ecc**2)*sinu 
      pos(3,1) = 0.d0 
      ! Compute the derivative of the position wrt u
      dposdu(1,1) = L**2/kappa2*(-sinu)
      dposdu(2,1) = L*G*cosu/kappa2
      dposdu(3,1) = 0.d0
      ! Compute the unit vector of the velocity and the norm 
      ! of the velocity
      norma          = sqrt(dot_product(dposdu(1:3,1),dposdu(1:3,1)))
      versvel(1:3,1) = dposdu(1:3,1)/norma
      normpos        = sqrt(dot_product(pos(1:3,1),pos(1:3,1)))
      normvel        = kappa*sqrt(2.d0/normpos-1.d0/aa) 
      ! Compute the velocity vector
      vel(1:3,1)     = normvel*versvel(1:3,1)
      ! Rotate everything in the 3d space
      tmpcar(1:3,1)  = matmul(Rot,pos(1:3,1))
      tmpcar(4:6,1)  = matmul(Rot,vel(1:3,1))    
      ! Store the result
      car(1:6) = tmpcar(1:6,1)
   end subroutine kep2car

   ! PURPOSE: Compute the cross product c = a x b
   !
   ! INPUT: 
   !    a, b : vectors of dimension 3  
   !
   ! OUTPUT:
   !       c : result of the cross product
   subroutine cross_prod(a,b,c) 
      real(kind=dkind), intent(in)  :: a(3), b(3)
      real(kind=dkind), intent(out) :: c(3)
      c(1) = a(2)*b(3) - a(3)*b(2) 
      c(2) = a(3)*b(1) - a(1)*b(3) 
      c(3) = a(1)*b(2) - a(2)*b(1) 
   end subroutine 

   !=================================
   !========== AVERAGING ============
   !=================================

   ! PURPOSE: Computes the integral average of a function using the trapezoid rule 
   !
   ! INPUT: 
   !    time : final integration time
   !  deltaT : stepsize of the trapezoid rule
   !     fun : vector containing the evaluations of the function to integrate
   ! npoints : length of the fun vector
   !
   ! OUTPUT:
   ! average : value of the integral average
   subroutine trapezoid_average(time, deltaT, fun, npoints, average)
      integer,          intent(in)  :: npoints
      real(kind=dkind), intent(in)  :: time, deltaT
      real(kind=dkind), intent(in)  :: fun(npoints)
      real(kind=dkind), intent(out) :: average
      ! end interface
      real(kind=dkind) :: integral
      integer          :: j
      ! Set the result to 0 
      integral = 0.d0
      ! Start adding the function evaluations
      do j=2, npoints-1
         integral = integral + fun(j)
      enddo
      ! Add the first and last points
      integral = integral + 0.5d0*(fun(1) + fun(npoints))
      ! Multiply for the step
      integral = integral*deltaT
      ! Divide by the total integration time
      average = integral/time
   end subroutine trapezoid_average

end module yarko_force
