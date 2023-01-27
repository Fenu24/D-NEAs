! ################################################################
! #     Authors: Marco Fenucci, Bojan Novaković, Dušan Marceta   #
! # Institution: University of Belgrade                          #
! #        Date: July 2022                                       #
! ################################################################
!
! PURPOSE: 
module population_model
   use used_const
   implicit none
   private
   ! Variables to store the Granvik model
   integer, parameter :: ncells = 730560
   real(kind=dkind)   :: gmb_model(ncells, 11)

   ! Public subroutines
   public :: gen_normal, gen_obliquity, gen_density_diameter, gen_density, gen_diameter, gen_albedo
   public :: gmb_load, gmb_search

   contains 

   ! =================================
   !       YARKOVSKY PARAMETERS       
   !         DISTRIBUTIONS            
   ! =================================

   ! PURPOSE: Generates a joined distribution of density and diameter from population model
   !
   ! INPUT:
   !   Hmean : mean value of the absolute magnitude 
   !    Hstd : error of the absolute magnitude 
   !  proute : source region probabilities from Granvik model 
   ! 
   ! OUTPUT:
   !     rho : density  [kg/m^3]
   !       D : diameter [m]
   !      pv : albedo
   subroutine gen_density_diameter(Hmean, Hstd, proute, rho, D, pv)
      real(kind=dkind), intent(in)  :: Hmean, Hstd
      real(kind=dkind), intent(in)  :: proute(7)
      real(kind=dkind), intent(out) :: rho
      real(kind=dkind), intent(out) :: D
      ! Albedo value
      real(kind=dkind), intent(out) :: pv
      ! end interface
      ! For lognormal distributions
      real(kind=dkind)              :: muC, muS, muX
      real(kind=dkind)              :: xxC, xxS, xxX
      real(kind=dkind)              :: ssC, ssS, ssX
      ! Values for C-, S, and X-complexes
      real(kind=dkind), parameter   ::   rhoC = 1200.d0
      real(kind=dkind), parameter   ::   rhoS = 2720.d0
      real(kind=dkind), parameter   ::   rhoX = 2350.d0
      real(kind=dkind), parameter   :: sigmaC = 300.d0
      real(kind=dkind), parameter   :: sigmaS = 540.d0
      real(kind=dkind), parameter   :: sigmaX = 520.d0
      ! Absolute magnitude value
      real(kind=dkind)              :: H
      ! Create the mu and sigma for lognormal distributions
      muC = log(rhoC)
      xxC = (1.d0 + sqrt(1.d0 + 4.d0 * sigmaC**2/exp(2.d0*muC)))/2.d0
      ssC = sqrt(log(xxC))
      muS = log(rhoS)
      xxS = (1.d0 + sqrt(1.d0 + 4.d0 * sigmaS**2/exp(2.d0*muS)))/2.d0
      ssS = sqrt(log(xxS))
      muX = log(rhoX)
      xxX = (1.d0 + sqrt(1.d0 + 4.d0 * sigmaX**2/exp(2.d0*muX)))/2.d0
      ssX = sqrt(log(xxX))
      ! Generate a value of the albedo, assuming Fenucci et al 2021 distribution model
      pv = gen_albedo(proute) 
      ! Generate a value of the absolute magnitude, assuming Gaussian distribution
      H  = gen_normal(Hmean, Hstd)
      ! Convert (pv, H) in D
      D  = pv2dia(pv, H)
      ! Generate the value of the density according to the complex, using lognormal distribution
      if(pv.le.0.1d0)then
         rho = gen_lognormal(muC, ssC)
      elseif( pv .gt. 0.1 .and. pv .le. 0.3)then
         rho = gen_lognormal(muS, ssS)
      else
         rho = gen_lognormal(muX, ssX)
      endif
   end subroutine

   ! INPUT:
   !   Hmean : mean value of the absolute magnitude 
   !    Hstd : error of the absolute magnitude 
   !  proute : source region probabilities from Granvik model 
   ! 
   ! OUTPUT:
   !     rho : density  [kg/m^3]
   !      pv : albedo
   subroutine gen_density(Hmean, Hstd, proute, rho, pv)
      real(kind=dkind), intent(in)  :: Hmean, Hstd
      real(kind=dkind), intent(in)  :: proute(7)
      real(kind=dkind), intent(out) :: rho
      real(kind=dkind), intent(out) :: pv
      ! end interface
      ! For lognormal distributions
      real(kind=dkind)              :: muC, muS, muX
      real(kind=dkind)              :: xxC, xxS, xxX
      real(kind=dkind)              :: ssC, ssS, ssX
      ! Values for C-, S, and X-complexes
      real(kind=dkind), parameter   ::   rhoC = 1200.d0
      real(kind=dkind), parameter   ::   rhoS = 2720.d0
      real(kind=dkind), parameter   ::   rhoX = 2350.d0
      real(kind=dkind), parameter   :: sigmaC = 300.d0
      real(kind=dkind), parameter   :: sigmaS = 540.d0
      real(kind=dkind), parameter   :: sigmaX = 520.d0
      ! Absolute magnitude value
      real(kind=dkind)              :: H
      ! Create the mu and sigma for lognormal distributions
      muC = log(rhoC)
      xxC = (1.d0 + sqrt(1.d0 + 4.d0 * sigmaC**2/exp(2.d0*muC)))/2.d0
      ssC = sqrt(log(xxC))
      muS = log(rhoS)
      xxS = (1.d0 + sqrt(1.d0 + 4.d0 * sigmaS**2/exp(2.d0*muS)))/2.d0
      ssS = sqrt(log(xxS))
      muX = log(rhoX)
      xxX = (1.d0 + sqrt(1.d0 + 4.d0 * sigmaX**2/exp(2.d0*muX)))/2.d0
      ssX = sqrt(log(xxX))
      ! Generate a value of the albedo, assuming Fenucci et al 2021 distribution model
      pv = gen_albedo(proute) 
      ! Generate a value of the absolute magnitude, assuming Gaussian distribution
      H  = gen_normal(Hmean, Hstd)
      ! Generate the value of the density according to the complex, using lognormal distribution
      if(pv.le.0.1d0)then
         rho = gen_lognormal(muC, ssC)
      elseif( pv .gt. 0.1 .and. pv .le. 0.3)then
         rho = gen_lognormal(muS, ssS)
      else
         rho = gen_lognormal(muX, ssX)
      endif
   end subroutine

   ! INPUT:
   !   Hmean : mean value of the absolute magnitude 
   !    Hstd : error of the absolute magnitude 
   !  proute : source region probabilities from Granvik model 
   ! 
   ! OUTPUT:
   !       D : diameter [m]
   !      pv : albedo
   subroutine gen_diameter(Hmean, Hstd, proute, D, pv)
      real(kind=dkind), intent(in)  :: Hmean, Hstd
      real(kind=dkind), intent(in)  :: proute(7)
      real(kind=dkind), intent(out) :: D
      real(kind=dkind), intent(out) :: pv
      ! end interface
      ! Absolute magnitude value
      real(kind=dkind)              :: H
      ! Generate a value of the albedo, assuming Fenucci et al 2021 distribution model
      pv = gen_albedo(proute) 
      ! Generate a value of the absolute magnitude, assuming Gaussian distribution
      H  = gen_normal(Hmean, Hstd)
      ! Convert (pv, H) in D
      D  = pv2dia(pv, H)
   end subroutine

   ! PURPOSE: Generates a distribution of albedo from NEO population models
   !
   ! INPUT:
   !  proute : source region probabilities from Granvik model 
   function gen_albedo(proute)
      real(kind=dkind) :: proute(7)
      real(kind=dkind) :: gen_albedo
      ! end interface
      real(kind=dkind)            :: p_albedo(3, 7)
      real(kind=dkind), parameter :: integral = 0.104526d0
      real(kind=dkind)            :: prob1, prob2
      real(kind=dkind)  :: R, pv, aux
      real(kind=dkind)  :: const, A, B, C
      integer :: j 
      ! Values from the Table 1 of Morbidelli et al 2020
      p_albedo(1, 1) = 0.120d0
      p_albedo(1, 2) = 0.144d0
      p_albedo(1, 3) = 0.294d0
      p_albedo(1, 4) = 0.021d0
      p_albedo(1, 5) = 0.501d0
      p_albedo(1, 6) = 0.399d0
      p_albedo(1, 7) = 1.0d0

      p_albedo(2, 1) = 0.558d0
      p_albedo(2, 2) = 0.782d0
      p_albedo(2, 3) = 0.557d0
      p_albedo(2, 4) = 0.113d0
      p_albedo(2, 5) = 0.452d0
      p_albedo(2, 6) = 0.200d0
      p_albedo(2, 7) = 0.0d0

      p_albedo(3, 1) = 0.322d0
      p_albedo(3, 2) = 0.074d0
      p_albedo(3, 3) = 0.149d0
      p_albedo(3, 4) = 0.866d0
      p_albedo(3, 5) = 0.047d0
      p_albedo(3, 6) = 0.461d0
      p_albedo(3, 7) = 0.0d0
      ! Compute the probabilities of
      ! pv < 0.1, 0.1 < pv < 0.3
      prob1 = 0.d0
      prob2 = 0.d0
      do j=1, 7
         prob1 = prob1 + p_albedo(1, j)*proute(j)
      enddo
      do j=1, 7
         prob2 = prob2 + p_albedo(2, j)*proute(j)
      enddo
      ! Compute the constant in front of the third part of the distribution
      const = 0.d0
      do j=1, 7
         const = const + p_albedo(3, j)*proute(j)
      enddo
      const = const/integral
      ! Start with the inverse sampling algorithm
      ! Generate a random number between 0 and 1
      call random_number(R)
      ! If               R <= prob1       -> category 1
      ! If prob1       < R <= prob1+prob2 -> category 2
      ! If prob1+prob2 < R <= 1           -> category 3
      ! NOTE: The CDF is defined piecewise and it needs to be inverted 
      !       case by case
      if(R.le.prob1)then
         pv = 0.1d0*R/prob1
      elseif(R.gt.prob1 .and. R.le.prob1+prob2)then
         pv = (R - prob1)*0.2d0/prob2 + 0.1d0
      else
         aux = (1.d0 - (R-prob1-prob2)*log(2.6d0)/(0.1d0*const))
         pv  = -log(aux)*0.1d0/log(2.6d0) + 0.3d0
      endif
      gen_albedo = pv
   end function 

   ! PURPOSE: Generates obliquity equal to 0 or 180 with equal probability
   ! 
   ! OUTPUT:
   !   gam : obliquity [deg]
   ! TODO: I have to fix this
   ! The cdf for cos gamma is
   !
   ! f = a y^3/3 + b y^2/2 + cy +a/3 - b/2 +c
   subroutine gen_obliquity(gam)
      real(kind=dkind), intent(out) :: gam
      ! end interface
      real(kind=dkind), parameter :: a = 1.12d0
      real(kind=dkind), parameter :: b = -0.32d0
      real(kind=dkind), parameter :: c = 0.13d0
      real(kind=dkind)   :: f, fp, cosg, cosgp
      real(kind=dkind)   :: x
      integer            :: j
      integer, parameter :: maxit = 100
      call random_number(x)
      cosg = -1.d0
      do j=1, maxit
         ! f  = a y^3/3 + b y^2/2 + cy +a/3 - b/2 +c
         ! f' = a y^2 + b y + c 
         f     = a*cosg**3/3.d0 + b*cosg**2/2.d0 + c*cosg + a/3.d0 - b/2.d0 + c - x 
         fp    = a*cosg**2 + b*cosg + c
         cosgp = cosg - f/fp
         if(abs(cosg-cosgp).lt.1.d-14)then
            exit
         endif
         cosg = cosgp
      enddo
      gam = acos(cosgp)*rad2deg
   end subroutine

   !TODO: gen_alpha

   ! PURPOSE: Converts (pv, H) in diameter. Diameter in output in meters
   !
   ! INPUT:
   !   pv : albedo
   !    H : absolute magnitude 
   function pv2dia(pv, H)
      real(kind=dkind) :: pv2dia
      real(kind=dkind) :: pv, H
      ! end interface
      pv2dia = 1000.d0*1329.d0*10.d0**(-H/5.d0)/sqrt(pv)
   end function

   ! ==================================
   !       GRANVIK MODEL INTERFACE
   ! ==================================

   ! PURPOSE: Load the Gravik et al. 2018 model into the global variable gmb_model
   subroutine gmb_load()
      ! end interface
      real(kind=dkind) :: x
      integer          :: i
      open(100, file='dat/gmb_model.dat', action='read')
      do i=1, ncells
         read(100, *) gmb_model(i, 1:4), x, x, gmb_model(i, 5:11)
      enddo
      close(100)
   end subroutine

   ! PURPOSE: search for the closest point to our orbital elements in the grid points
   !
   ! INPUT:
   !       a : semi-major axis [au]
   !       e : eccentricity 
   !       i : inclination     [deg] 
   ! 
   ! OUTPUT:
   !  proute : source region probabilities from Granvik model 
   subroutine gmb_search(a, e, i, H, proute)
      real(kind=dkind), intent(in)  :: a, e, i, H
      real(kind=dkind), intent(out) :: proute(7)
      ! end interface
      real(kind=dkind) :: ax, ex, ix, Hx
      real(kind=dkind) :: norm, minNorm
      integer :: j, jmin 
      ! Put absurdly big initial min values   
      minNorm = 1e9
      jmin    = -1
      ! Search for the point of the grid closer to our orbital elements
      do j=1, ncells
         ax = gmb_model(j, 1)
         ex = gmb_model(j, 2)
         ix = gmb_model(j, 3)
         Hx = gmb_model(j, 4)
         norm = sqrt( (a-ax)**2 + (e-ex)**2 + (i-ix)**2 + (H-Hx)**2 )
         if(norm .lt. minNorm)then
            minNorm = norm
            jmin    = j
         endif
      enddo
      ! Take the probabilities of coming from each route
      ! nu6     5/2     2/1     HUN     3/1     PHO     JFC 
      proute(1:7) = gmb_model(jmin, 5:11)
   end subroutine 

   ! ==================================
   !         FOR STATISTICS
   ! ==================================

   ! PURPOSE: Generate normal distribution
   !
   ! INPUT:
   !    mu : mean value
   ! sigma : standard deviation
   function gen_normal(mu, sigma)
      real(kind=dkind) :: gen_normal
      real(kind=dkind) :: mu, sigma
      ! end interface
      real(kind=dkind)  :: x
      call gauss_random(x)
      gen_normal = sigma*x + mu
   end function

   ! PURPOSE: Generate lognormal distribution
   !
   ! INPUT:
   !    mu : mean value
   ! sigma : standard deviation
   function gen_lognormal(mu, sigma)
      real(kind=dkind) :: gen_lognormal
      real(kind=dkind) :: mu, sigma
      ! end interface
      real(kind=dkind) :: x
      x = gen_normal(mu, sigma)
      gen_lognormal = exp(x)
   end function

   ! PURPOSE: Generate a random number in (a,b)
   !
   ! INPUT:
   !  a,b : extrema of the interval
   function gen_random(a, b)
      real(kind=dkind) :: a, b
      real(kind=dkind) :: gen_random
      ! end interface
      real(kind=dkind) :: x
      call random_number(x)
      gen_random = (b-a)*x + a
   end function

   ! PURPOSE: Computes the mean of a vector
   !
   ! INPUT:
   !   x : vector with a distribution
   function mean(x)
      real(kind=dkind) :: x(:)
      real(kind=dkind) :: mean
      ! end interface
      integer :: n
      n    = size(x)
      mean = sum(x(1:n))/n
   end function

   ! PURPOSE: Computes the std of a vector
   !
   ! INPUT:
   !   x : vector with a distribution
   function std(x)
      real(kind=dkind) :: x(:)
      real(kind=dkind) :: std 
      ! end interface
      integer          :: n
      real(kind=dkind) :: mm
      n   = size(x)
      mm  = mean(x) 
      std = sqrt(sum((x(1:n)-mm)**2)/n)
   end function

   ! PURPOSE: Generates a random number with normal distribution and Box-Muller algorithm 
   ! 
   ! OUTPUT:
   !  u : normal distributed random variable
   subroutine gauss_random(u)
      real(kind=dkind), intent(out) :: u
      ! end interface
      real(kind=dkind) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u = sqrt(-2.0*log(u1))*cos(twopi*u2)
   end subroutine

end module
