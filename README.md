# Dynamics of the 2D Bose–Hubbard Model

This repository contains the Fortran codes developed during my MSc Physics project, **“Dynamics of Quantum Phase Transition in 2D Bose–Hubbard Model.”**

The project investigates the dynamics of the quantum phase transition from a **Mott-insulating phase to a superfluid phase** using the mean-field approximation and numerical simulations on a 2D optical lattice.

## Project Overview

The simulations involve:

* Self-consistent mean-field calculations of the ground state
* Numerical diagonalization of the single-site Hamiltonian
* Time evolution across the quantum phase transition using a Runge–Kutta method
* Simulations on a 70 × 70 lattice with periodic boundary conditions
* Calculation of the superfluid order parameter and particle density
* Analysis of spatial domain growth during the transition
* Analysis of freeze-out dynamics and Kibble–Zurek scaling

## Repository Structure

* **`phase transition_static/`** — Static mean-field calculation used to study the phase transition.
* **`dynamics/`** — Main program for the time-dependent simulation and quench dynamics.
* **`Subroutines/`** — Fortran subroutines used by the simulation programs, including Hamiltonian construction and Runge–Kutta time evolution.
* **`common.h`** — Shared definitions used by the Fortran programs.
* **`Thesis report/`** — MSc thesis describing the theoretical background, numerical methods, and results.

## Numerical Methods

The calculations use the mean-field approximation for the 2D Bose–Hubbard model. The time-dependent evolution is performed using a Runge–Kutta integration scheme, with the system driven through the quantum critical region by a linear quench of the hopping parameter.

## Results

The simulations were used to study the development of the superfluid order parameter, spatial domain growth, and freeze-out behaviour during the quantum phase transition. The resulting scaling behaviour was compared with the Kibble–Zurek mechanism.

## Author

**Apeksha Phadte**
MSc Physics (Computational Physics), Goa University
