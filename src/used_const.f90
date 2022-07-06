! ################################################################
! #      Authors: Marco Fenucci, Bojan Novakovic, Dusan Marceta  #
! # Institution: University of Belgrade                          #
! #        Date: July 2022                                       #
! ################################################################
!
! PURPOSE: This module contains mathematical and physical constants
!
module used_const
   implicit none
   private
   ! For FORTRAN programming
   integer, parameter         :: dkind = kind(1.d0)
   ! Mathematical Constants 
   real(kind=dkind),parameter :: pi      = 3.14159265358979d0 
   real(kind=dkind),parameter :: twopi   = 2.d0*pi
   real(kind=dkind),parameter :: deg2rad = pi/180.d0
   real(kind=dkind),parameter :: rad2deg = 180.d0/pi
   ! ===========================================
   ! Physical parameters: standard units [kg, m, s, K]
   ! ===========================================
   ! Universal gravitational constant [m3 kg-1 s-2]
   real(kind=dkind), parameter :: uGc = 6.67408d-11 
   ! Gravitational parameters G*m [m3 s-2]
   real(kind=dkind), parameter :: gmSun = 1.32712440018d20
   ! Stefan-Boltzmann constant [W m-2 K-4] = [kg s-3 K-4]
   real(kind=dkind), parameter :: sBoltz = 5.670374419d-8
   ! Solar luminosity [W] = [kg m2 s-3]
   real(kind=dkind), parameter :: lumSun = 3.828d26
   ! Astronomical Unit [m] and conversion factors from/to meters
   real(kind=dkind), parameter :: AU   = 149597870700.d0
   real(kind=dkind), parameter :: au2m = AU
   real(kind=dkind), parameter :: m2au = 1.d0/AU
   ! Speed of light [m s-1]
   real(kind=dkind), parameter :: clight = 299792458.d0
   ! Time conversions
   real(kind=dkind), parameter :: h2s   = 3600.d0
   real(kind=dkind), parameter :: s2h   = 1.d0/h2s
   real(kind=dkind), parameter :: d2s   = 24.d0*h2s
   real(kind=dkind), parameter :: s2d   = 1.d0/d2s
   real(kind=dkind), parameter :: y2s   = 365.25d0*d2s
   real(kind=dkind), parameter :: s2y   = 1.d0/y2s
   real(kind=dkind), parameter :: yr2my = 1.d-6
   real(kind=dkind), parameter :: s2my  = s2y*yr2my
   ! Set the variables to public
   public :: dkind, pi, twopi, deg2rad, rad2deg
   public :: uGc, gmSun, sBoltz, lumSun, au, au2m, m2au, clight
   public :: h2s, s2h, d2s, s2d, y2s, s2y, yr2my, s2my
end module used_const


