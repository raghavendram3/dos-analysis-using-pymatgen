# VASP Density of States (DOS) Analysis Using Pymatgen

## Overview

This code analyzes the electronic structure properties of a material using VASP (Vienna Ab initio Simulation Package) output data. It specifically focuses on d-orbital band characteristics for both spin-up and spin-down electrons at a particular atomic site. However, the p-orbital band characteristics could be obtained by making slight modifications.

## Purpose

The script calculates various statistical and electronic properties of d-bands:
- **Band Center**: The energy-weighted average position of the d-band
- **Band Filling**: The fraction of d-orbital states that are occupied
- **Band Skewness**: Asymmetry of the d-band distribution
- **Band Kurtosis**: Measure of the "tailedness" of the d-band distribution
- **Band Width**: Energy span of the d-band in the occupied region

## Dependencies

- **pymatgen**: Materials analysis library for reading VASP files and electronic structure analysis
- **matplotlib.pyplot**: (imported but not used in current version)

## Input Files

- **vasprun.xml**: VASP output file containing electronic structure data, recommended by pymatgen for DOS analysis

## Key Parameters

- **Orbital Type**: p/d-orbitals
- **Energy Range**: 
  - For band center, skewness, and kurtosis: -15 to +15 eV relative to Fermi level
  - For band width: -15 to 0 eV (occupied states only)
- **Spin Channels**: Both spin-up and spin-down are analyzed separately

## Output

The code prints:
1. Complete crystal structure information
2. D-band centers for both spin channels
3. Kurtosis values for both spin channels  
4. Skewness values for both spin channels
5. D-band filling fractions for both spin channels
6. D-band widths for both spin channels

## Physical Significance

These parameters are important for understanding:
- **Catalytic activity**: D-band center correlates with adsorption energies
- **Magnetic properties**: Spin-resolved analysis reveals magnetic moments
- **Electronic conductivity**: Band filling and width affect transport properties
- **Chemical bonding**: Band shape parameters indicate hybridization effects
