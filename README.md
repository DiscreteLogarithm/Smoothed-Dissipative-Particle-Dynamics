# Smoothed Dissipative Particle Dynamics

An isothermal implementation of SDPD for water based on Español and Revenga, Phys. Rev. E 67, 026705 (2003)

## Installation

* Put the patch file in your top-level LAMMPS directory, where the `LICENSE` and `README` files are.

* Apply the patch by typing the following command in your top-level LAMMPS directory:

        patch -p1 < patch.sdpd

* Then include USER-SPH package by typing the following commands in `src` directory:

        make yes-user-sph

* The RIGID package is a prerequisite for `fix rigid/meso`. If you like to be able to integrate rigid bodies
motion, include this as well:

        make yes-rigid

* Now include the USER-SDPD package

        make yes-user-sdpd

* Then build again as usual.

### Improving performane with Zest

Most of the time in SDPD simulations are spent producing random gaussian numbers. By default LAMMPS uses the
Box-Muller algorithm to produce random gaussian numbers which is inefficient. The Zest library can produce
random gaussian numbers much more efficiently. Using Zest, SDPD simulations are sped up y a factor of about 3.
Obtain Zest from https://github.com/DiscreteLogarithm/Zest , then copy the header files either into the `src`
directory of the LAMMPS or into your system's include path. Then compile LAMMPS with `-DUSE_ZEST` flag by
adding `-DUSE_ZEST` to the `LMP_INC` variable of your preferred Makefile.

LAMMPS uses an old integer random number generator - RANMAR - which has a relatively low period compared to
modern standards. When LAMMPS is compiled with Zest the 64-bit version of the Mersenne twister is used
which is the de-facto gold standard of uniform random number generators.

## Uninstallation

* Use following command in `src` directory:

        make no-user-sdpd

* Then build again

