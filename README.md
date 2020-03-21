# CompressionalRossbyWave
This module contains a subroutine generating analytical solutions for the compressional Rossby waves.

The default test uses the following special treatments:
1. Generalized equatorial f-plane, where the traditional and nontraditional Coriolis parameters are 0 and 2*omega.
2. Fast planetary rotation rate of omega = 6.973339d-3 s^-1
3. Barotropic ideal gas, where the Poisson constant (rd/cp) and heat capacity ratio (cp/cv) are 0 and 1.

If width = 2.d6 m, and depth = 12500 m (anelastic) or 12721 m (fully compressible), the wave period will be 86400 s.
