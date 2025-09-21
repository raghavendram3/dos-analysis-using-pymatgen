#!/usr/bin/env python3
"""
VASP Density of States (DOS) Analysis for d-Band Electronic Properties

This script analyzes the electronic structure properties of d-orbitals from VASP calculations,
computing various statistical parameters for both spin-up and spin-down channels.

Dependencies: pymatgen, vasprun
"""

import matplotlib.pyplot as plt
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType


def analyze_d_band_properties(vasprun_file_path: str, target_site_index: int = 161):
    """
    Analyze d-band electronic properties for a specific atomic site.
    
    Parameters:
    -----------
    vasprun_file_path : str
        Path to the vasprun.xml file from VASP calculation
    target_site_index : int
        Index of the atomic site to analyze (0-indexed)
    
    Returns:
    --------
    dict: Dictionary containing all calculated d-band properties
    """
    
    # Configuration parameters
    ENERGY_RANGE_FULL = [-15, 15]      # eV range for center, skewness, kurtosis
    ENERGY_RANGE_OCCUPIED = [-15, 0]   # eV range for band width (occupied states)
    ORBITAL_TYPE = OrbitalType.d        # Focus on d-orbitals
    
    print(f"Reading VASP data from: {vasprun_file_path}")
    print(f"Analyzing site index: {target_site_index}")
    print("-" * 50)
    
    # Read VASP output data
    try:
        vasprun_data = Vasprun(vasprun_file_path)
        fermi_energy = vasprun_data.efermi
        complete_dos_data = vasprun_data.complete_dos
        crystal_structure = complete_dos_data.structure
        
        print(f"Fermi energy: {fermi_energy:.4f} eV")
        print(f"Total number of sites: {len(crystal_structure)}")
        
    except Exception as e:
        print(f"Error reading VASP file: {e}")
        return None
    
    # Validate site index
    if target_site_index >= len(crystal_structure):
        print(f"Error: Site index {target_site_index} exceeds number of sites ({len(crystal_structure)})")
        return None
    
    target_site = crystal_structure[target_site_index]
    print(f"Target site element: {target_site.species_string}")
    print(f"Target site coordinates: {target_site.frac_coords}")
    
    # Calculate d-band properties for both spin channels
    print("\nCalculating d-band properties...")
    
    # Band centers (energy-weighted average position)
    d_band_center_spin_up = complete_dos_data.get_band_center(
        sites=[target_site],
        spin=Spin.up,
        erange=ENERGY_RANGE_FULL,
        band=ORBITAL_TYPE
    )
    
    d_band_center_spin_down = complete_dos_data.get_band_center(
        sites=[target_site],
        spin=Spin.down,
        erange=ENERGY_RANGE_FULL,
        band=ORBITAL_TYPE
    )
    
    # Band filling (fraction of occupied states)
    d_band_filling_spin_up = complete_dos_data.get_band_filling(
        sites=[target_site],
        spin=Spin.up,
        band=ORBITAL_TYPE
    )
    
    d_band_filling_spin_down = complete_dos_data.get_band_filling(
        sites=[target_site],
        spin=Spin.down,
        band=ORBITAL_TYPE
    )
    
    # Band skewness (asymmetry of distribution)
    d_band_skewness_spin_up = complete_dos_data.get_band_skewness(
        sites=[target_site],
        spin=Spin.up,
        erange=ENERGY_RANGE_FULL,
        band=ORBITAL_TYPE
    )
    
    d_band_skewness_spin_down = complete_dos_data.get_band_skewness(
        sites=[target_site],
        spin=Spin.down,
        erange=ENERGY_RANGE_FULL,
        band=ORBITAL_TYPE
    )
    
    # Band kurtosis (measure of tail heaviness)
    d_band_kurtosis_spin_up = complete_dos_data.get_band_kurtosis(
        sites=[target_site],
        spin=Spin.up,
        erange=ENERGY_RANGE_FULL,
        band=ORBITAL_TYPE
    )
    
    d_band_kurtosis_spin_down = complete_dos_data.get_band_kurtosis(
        sites=[target_site],
        spin=Spin.down,
        erange=ENERGY_RANGE_FULL,
        band=ORBITAL_TYPE
    )
    
    # Band width (energy span of occupied states)
    d_band_width_spin_up = complete_dos_data.get_band_width(
        sites=[target_site],
        spin=Spin.up,
        erange=ENERGY_RANGE_OCCUPIED,
        band=ORBITAL_TYPE
    )
    
    d_band_width_spin_down = complete_dos_data.get_band_width(
        sites=[target_site],
        spin=Spin.down,
        erange=ENERGY_RANGE_OCCUPIED,
        band=ORBITAL_TYPE
    )
    
    # Store results in dictionary
    results = {
        'fermi_energy': fermi_energy,
        'site_index': target_site_index,
        'site_element': target_site.species_string,
        'band_centers': {
            'spin_up': d_band_center_spin_up,
            'spin_down': d_band_center_spin_down
        },
        'band_fillings': {
            'spin_up': d_band_filling_spin_up,
            'spin_down': d_band_filling_spin_down
        },
        'band_skewness': {
            'spin_up': d_band_skewness_spin_up,
            'spin_down': d_band_skewness_spin_down
        },
        'band_kurtosis': {
            'spin_up': d_band_kurtosis_spin_up,
            'spin_down': d_band_kurtosis_spin_down
        },
        'band_widths': {
            'spin_up': d_band_width_spin_up,
            'spin_down': d_band_width_spin_down
        }
    }
    
    return results, crystal_structure


def print_analysis_results(results: dict, crystal_structure):
    """
    Print formatted analysis results.
    
    Parameters:
    -----------
    results : dict
        Dictionary containing calculated properties
    crystal_structure : Structure
        Pymatgen Structure object
    """
    
    print("\n" + "="*60)
    print("ELECTRONIC STRUCTURE ANALYSIS RESULTS")
    print("="*60)
    
    # Crystal structure information
    print("\nCRYSTAL STRUCTURE:")
    print(crystal_structure)
    
    # Site information
    print(f"\nTARGET SITE ANALYSIS:")
    print(f"Site index: {results['site_index']}")
    print(f"Element: {results['site_element']}")
    print(f"Fermi energy: {results['fermi_energy']:.4f} eV")
    
    # D-band properties
    print(f"\nD-BAND ELECTRONIC PROPERTIES:")
    print(f"{'Property':<20} {'Spin-Up':<15} {'Spin-Down':<15}")
    print("-" * 50)
    
    print(f"{'Band Center (eV)':<20} "
          f"{results['band_centers']['spin_up']:<15.4f} "
          f"{results['band_centers']['spin_down']:<15.4f}")
    
    print(f"{'Band Filling':<20} "
          f"{results['band_fillings']['spin_up']:<15.4f} "
          f"{results['band_fillings']['spin_down']:<15.4f}")
    
    print(f"{'Band Skewness':<20} "
          f"{results['band_skewness']['spin_up']:<15.4f} "
          f"{results['band_skewness']['spin_down']:<15.4f}")
    
    print(f"{'Band Kurtosis':<20} "
          f"{results['band_kurtosis']['spin_up']:<15.4f} "
          f"{results['band_kurtosis']['spin_down']:<15.4f}")
    
    print(f"{'Band Width (eV)':<20} "
          f"{results['band_widths']['spin_up']:<15.4f} "
          f"{results['band_widths']['spin_down']:<15.4f}")
    
    # Magnetic properties
    magnetic_moment = (results['band_fillings']['spin_up'] - 
                      results['band_fillings']['spin_down'])
    print(f"\nMAGNETIC PROPERTIES:")
    print(f"d-electron magnetic moment: {magnetic_moment:.4f} μB")
    
    # Analysis insights
    print(f"\nANALYSIS INSIGHTS:")
    center_up = results['band_centers']['spin_up']
    center_down = results['band_centers']['spin_down']
    
    if abs(center_up - center_down) > 0.1:
        print(f"• Significant spin-splitting detected ({abs(center_up - center_down):.3f} eV)")
    
    if abs(magnetic_moment) > 0.1:
        print(f"• Magnetic material with d-band moment: {magnetic_moment:.3f} μB")
    
    avg_center = (center_up + center_down) / 2
    if avg_center < -2.0:
        print(f"• Deep d-band center ({avg_center:.3f} eV) suggests strong bonding")
    elif avg_center > -1.0:
        print(f"• Shallow d-band center ({avg_center:.3f} eV) suggests weak bonding")


def main():
    """Main analysis function."""
    
    # Configuration
    VASPRUN_FILE_PATH = "./vasprun.xml"
    TARGET_SITE_INDEX = 161  # 0-indexed site number
    
    # Perform analysis
    analysis_results = analyze_d_band_properties(
        vasprun_file_path=VASPRUN_FILE_PATH,
        target_site_index=TARGET_SITE_INDEX
    )
    
    if analysis_results is None:
        print("Analysis failed. Please check input file and parameters.")
        return
    
    results, structure = analysis_results
    
    # Print results
    print_analysis_results(results, structure)


if __name__ == "__main__":
    main()
