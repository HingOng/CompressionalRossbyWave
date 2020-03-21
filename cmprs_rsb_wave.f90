MODULE cmprs_rsb_wave

!=======================================================================
!
! This module contains a subroutine generating analytical solutions for
! the compressional Rossby waves.
!
! The default test uses the following special treatments:
! 1. Generalized equatorial f-plane, where the traditional and
!    nontraditional Coriolis parameters are 0 and 2*omega.
! 2. Fast planetary rotation rate of omega = 6.973339d-3 s^-1
! 3. Barotropic ideal gas, where the Poisson constant (rd/cp) and heat
!    capacity ratio (cp/cv) are 0 and 1.
!
! If width = 2.d6 m, and depth = 12500 m (anelastic) or 12721 m (fully
! compressible), the wave period will be 86400 s.
!
!=======================================================================

  IMPLICIT NONE

!=======================================================================
! physical constants
!=======================================================================
  REAL(8), PARAMETER ::          &
    g     = 9.80616d0,           & ! acceleration due to gravity (m s^-2)
    rd    = 287.0d0,             & ! gas constant for dry air (J K^-1 kg^-1)
    rd_cp = 0.d0,                & ! Poisson constant (rd/cp), i.e., kappa
    cp_cv = 1.d0,                & ! heat capacity ratio (cp/cv), i.e., gamma
    p0    = 1.d5,                & ! reference pressure (Pa)
    pi    = 3.141592653589793d0, & ! pi
    omega = 6.973339d-3            ! planetary rotation rate (s^-1)

!=======================================================================
!  test case parameters
!=======================================================================
  REAL(8), PARAMETER ::            &
    t0      = 311.0d0                ! basic-state temperature (K)

  CONTAINS

  SUBROUTINE cmprs_rsb_wave_test(elastic, amp, wid, dep, x, z, tim, u, v, w, pb, t, phi, pp) &
             BIND(c, name = "cmprs_rsb_wave_test")

    IMPLICIT NONE

!=======================================================================
!   input variables
!=======================================================================
    REAL(8), INTENT(IN) :: &
      elastic, & ! fully compressible (1.) or anelastic (0.)
      amp,     & ! amplitude scaling factor
      wid,     & ! domain width i.e. zonal wavelength (m)
      dep,     & ! domain depth i.e. half vertical wavelength (m)
      x,       & ! Cartesian x-coordinate (m)
      z,       & ! Cartesian z-coordinate (m)
      tim        ! time (s)

!=======================================================================
!   output variables
!=======================================================================
    REAL(8), INTENT(OUT) :: &
      u,       & ! zonal velocity (m s^-1)
      v,       & ! meridional velocity (m s^-1)
      w,       & ! vertical velocity (m s^-1)
      pb,      & ! basic-state pressure (Pa)
      t,       & ! temperature (K)
      phi,     & ! perturbation Exner function times cp*theta (J kg^-1)
      pp         ! perturbation pressure (Pa)

!=======================================================================
!   local variables
!=======================================================================
    REAL(8) :: rH, k_wave, m_wave, omega_wave, u0, phi0, w0, phsh_phi, phsh_u

!=======================================================================
!   Generate wave parameters
!=======================================================================
    rH         = g/(rd*t0)                                                         ! inverse density scale height (m^-1)
    k_wave     = 2.d0*pi/wid                                                       ! wavenumber in zonal (m^-1)
    m_wave     = pi/dep                                                            ! wavenumber in vertical (m^-1)
    omega_wave = 2.d0*omega*k_wave*rH/( k_wave**2 +m_wave**2 +0.25d0*rH**2 &       ! wave frequency (s^-1)
                                       +elastic*4.d0*omega**2/(rd*t0) )
    u0         = 0.09d0*amp                                                        ! wave amplitude in u (m s^-1)
    phi0       = u0*(4.d0*omega**2-omega_wave**2)*( omega**2*rH**2 &               ! wave amplitude in phi (J kg^-1)
                     +4.d0*omega**2*m_wave**2 -2.d0*omega*omega_wave*k_wave*rH &
                     +omega_wave**2*k_wave**2 )**(-0.5)
    w0         = phi0*k_wave*(4.d0*omega**2-omega_wave**2)**(-0.5)                 ! wave amplitude in w (m s^-1)
    phsh_phi   = atan((0.5d0*omega_wave*rH-2.d0*omega*k_wave)/(m_wave*omega_wave)) ! phase shift for phi
    phsh_u     = phsh_phi + atan(2.d0*omega*m_wave/(omega*rH-omega_wave*k_wave))   ! phase shift for u

!=======================================================================
!   Generate analytical solutions
!=======================================================================
    u   = u0  *exp(0.5d0*z*rH)*cos(m_wave*z+phsh_u  )*cos(k_wave*x-omega_wave*tim)
    v   = 0.d0
    w   = w0  *exp(0.5d0*z*rH)*sin(m_wave*z         )*sin(k_wave*x-omega_wave*tim)
    pb  = p0  *exp(     -z*rH)
    t   = t0
    phi = phi0*exp(0.5d0*z*rH)*cos(m_wave*z+phsh_phi)*cos(k_wave*x-omega_wave*tim)
    pp  = phi *pb/(rd*t0)

  ENDSUBROUTINE cmprs_rsb_wave_test

ENDMODULE cmprs_rsb_wave
