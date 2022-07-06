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
