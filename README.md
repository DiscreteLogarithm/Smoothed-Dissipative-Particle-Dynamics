# Smoothed Dissipative Particle Dynamics

An isothermal implementation of SDPD based on Español and Revenga, Phys. Rev. E 67, 026705 (2003)

## Installation

* Put the patch file in your top-level LAMMPS directory, where the `LICENSE` and `README` files are.

* Apply the patch by typing the following command in your top-level LAMMPS directory:

        patch -p1 < patch.sdpd

* Then include this package by typing the following command in `src` directory:

        make yes-user-sph

* Then build again as usual.

## Uninstallation

* Use following command in `src` directory:

        make no-user-sph

* Then build again

