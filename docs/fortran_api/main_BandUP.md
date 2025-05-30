# `main_BandUP.f90`

## Overview

This file contains the main program `BandUP_main`, which serves as the primary executable for the BandUP application. BandUP (Band Unfolding code for Plane-wave based calculations) is designed to perform band unfolding calculations, a technique used in solid-state physics to analyze the electronic band structure of supercells by relating them back to the primitive cell's band structure. This program orchestrates the entire workflow, from reading user inputs and wavefunction data to performing the unfolding calculations and generating output files with the results. It is authored by Paulo V. C. Medeiros.

## Key Components

### Program `BandUP_main`

This is the main entry point of the BandUP application. Its execution flow can be summarized as follows:

1.  **Initialization**:
    *   Initializes timing routines (`CALL initialize(times)`).
    *   Sets OpenMP dynamic behavior (`CALL omp_set_dynamic(.FALSE.)`).
    *   Prints welcome messages, including the package version.
    *   Parses command-line arguments (`CALL get_commline_args(args)`).

2.  **Input Data Reading and Setup**:
    *   Reads supercell crystal structure (`CALL get_crystal_from_file(crystal_SC, ...)`). If not provided, it's derived from the wavefunction file.
    *   Reads wavefunction data (metadata initially, `read_coeffs=.FALSE.`) (`CALL read_wavefunction(wf, ...)`).
    *   Determines primitive cell for the supercell (`CALL get_prim_cell(crystal_SC, ...)`).
    *   Reads primitive cell crystal structure from file (`CALL get_crystal_from_file(crystal_pc, ...)`).
    *   Determines primitive cell for the input primitive cell (`CALL get_prim_cell(crystal_pc, ...)`).
    *   Calculates reciprocal lattice for the primitive cell (`CALL get_rec_latt(...)`).
    *   Writes information about the attempted primitive cell association (`CALL write_attempted_pc_assoc_with_input_unit_cell_and_SC(...)`).
    *   Verifies commensurability between primitive cell and supercell (`CALL verify_commens(...)`).
    *   Analyzes symmetry relations between primitive cell and supercell (`CALL analise_symm_pc_SC(...)`).
    *   Reads user-selected k-point paths in the primitive cell Brillouin zone (`CALL read_pckpts_selected_by_user(...)`).
    *   Determines all unique primitive cell Brillouin zone directions to be used for the Electronic Band Structure (EBS) based on symmetry (`CALL get_pcbz_dirs_2b_used_for_EBS(...)`).
    *   Prints symmetry analysis for these selected directions.
    *   Defines the specific pc-kpoints to be checked along these directions (`CALL define_pckpts_to_be_checked(...)`).
    *   Gets the list of SC K-points from the input/wavefunction file (`CALL get_list_of_SCKPTS(...)`).
    *   Calculates the Geometric Unfolding Relations (GUR) (`CALL get_geom_unfolding_relations(...)`). These relations map pc-kpoints to SC K-points.
    *   Prints messages about GUR determination and checks for success.
    *   Prints the detailed GUR.
    *   Reads energy range information for the band search (`CALL read_energy_info_for_band_search(...)`).
    *   Creates an energy grid (`CALL real_seq(...)`).
    *   Prints final messages before starting the main unfolding loop.
    *   Allocates memory for the `delta_N` (unfolded quantities) variable (`CALL allocate_UnfoldedQuantities(delta_N, ...)`).

3.  **Main Unfolding Loop**:
    *   Iterates over each SC K-point (`i_SCKPT`).
    *   Iterates over each selected primitive cell Brillouin zone direction (`i_selec_pcbz_dir`).
    *   Iterates over symmetry-equivalent needed directions (`i_needed_dirs`).
    *   Iterates over each primitive cell k-point (`ipc_kpt`) along the current direction.
    *   Checks if the current `pckpt` folds from the current `i_SCKPT` based on GUR.
    *   If it folds:
        *   Updates GUR indices (`CALL update_GUR_indices(...)`).
        *   If wavefunction coefficients for `i_SCKPT` are not yet loaded, reads them (`CALL read_wavefunction(wf, ..., read_coeffs=.TRUE.)`).
        *   Selects the relevant plane-wave coefficients from the wavefunction that contribute to the spectral weight at the current pc-kpt (`CALL select_coeffs_to_calc_spectral_weights(...)`).
        *   If coefficients are found, performs the unfolding calculation (`CALL perform_unfolding(...)` which calculates `delta_N`, spectral functions, etc.).
        *   If no relevant coefficients are found, prints a message and potentially stops.

4.  **Output Generation**:
    *   Gathers the calculated unfolded quantities (`delta_N`) and prepares them for output. This involves creating versions strictly along user-selected directions and symmetry-averaged versions (`CALL get_delta_Ns_for_output(...)`).
    *   Saves the results to output files and prints goodbye messages (`CALL say_goodbye_and_save_results(...)`).
    *   Prints final timing information (`CALL print_final_times(times)`).

## Important Variables/Constants

The program uses several key variables to store data and manage the calculation flow. Many of these are instances of derived types defined in `constants_and_types_mod`.

*   `wf` (type `pw_wavefunction`): Stores the plane-wave wavefunction data (coefficients, energies, k-points, lattice vectors).
*   `delta_N` (type `UnfoldedQuantities`): Stores the primary results of the unfolding – the unfolded band character (spectral weights) for each pc-kpt and energy.
*   `delta_N_only_selected_dirs` (type `UnfoldedQuantitiesForOutput`): Stores unfolded quantities formatted for output, strictly along user-selected k-point paths.
*   `delta_N_symm_avrgd_for_EBS` (type `UnfoldedQuantitiesForOutput`): Stores symmetry-averaged unfolded quantities for the Electronic Band Structure (EBS).
*   `GUR` (type `geom_unfolding_relations_for_each_SCKPT`): Contains the vital geometric unfolding relations mapping pc-kpoints to SC K-points.
*   `all_dirs_used_for_EBS_along_pcbz_dir` (allocatable array of type `irr_bz_directions`): List of all irreducible Brillouin zone directions used for generating the EBS.
*   `pckpts_to_be_checked` (type `selected_pcbz_directions`): The specific pc-kpoints along selected directions that will be analyzed.
*   `crystal_SC` (type `crystal_3D`): Holds information about the supercell crystal structure.
*   `crystal_pc` (type `crystal_3D`): Holds information about the primitive cell crystal structure.
*   `list_of_SCKPTS` (allocatable array of type `vec3d`): List of k-points read from the supercell calculation.
*   `times` (type `timekeeping`): Stores timing information for different stages of the program.
*   `args` (type `comm_line_args`): Stores command-line arguments (not explicitly declared here but passed to subroutines that `USE cla_wrappers`).
*   `k_starts`, `k_ends` (allocatable real(dp) array): Start and end points of user-defined k-paths.
*   `energy_grid` (allocatable real(dp) array): Grid of energy values for calculating spectral functions.
*   `e_fermi` (real(dp)): Fermi energy.
*   `selected_coeff_indices` (allocatable integer array): Indices of plane-wave coefficients selected for spectral weight calculation.

## Usage Examples

`BandUP_main` is an executable. It is typically run from the command line with various arguments specifying input files (wavefunction, cell definitions, k-paths, energy ranges) and output file names. The specific command-line arguments are parsed by routines from `cla_wrappers_mod.f90`.

A conceptual usage might be:

```bash
./BandUP_main --wf_file WAVECAR --prim_cell_file POSCAR_PC --supercell_file POSCAR_SC --pc_kpts_file KPOINTS_PC --energies_file E_RANGE.in ... [other options]
```

(Note: The exact command-line flags are defined by the `comm_line_args` type and handled by `cla_wrappers_mod.f90`).

## Dependencies and Interactions

`BandUP_main` depends on several other modules within the BandUP project:

*   `constants_and_types`: For fundamental constants, precision kinds, and derived data types (e.g., `pw_wavefunction`, `UnfoldedQuantities`, `crystal_3D`).
*   `cla_wrappers`: For parsing command-line arguments (`CALL get_commline_args(args)`).
*   `crystals`: For crystal structure creation and manipulation (e.g., `CALL create_crystal(...)`, `CALL get_rec_latt(...)`).
*   `symmetry`: For symmetry analysis (e.g., `CALL get_prim_cell(...)`, `CALL analise_symm_pc_SC(...)`, `CALL get_pcbz_dirs_2b_used_for_EBS(...)`).
*   `lists_and_seqs`: For utility functions like generating sequences (e.g., `CALL real_seq(...)` for `energy_grid`).
*   `time`: For initializing and managing timing (`CALL initialize(times)`, `CALL print_final_times(times)`).
*   `general_io`: For general input/output operations (used by sub-modules).
*   `io_routines`: For specific I/O tasks like reading input files and writing results (e.g., `CALL get_crystal_from_file(...)`, `CALL read_wavefunction(...)`, `CALL read_pckpts_selected_by_user(...)`, `CALL read_energy_info_for_band_search(...)`, `CALL print_welcome_messages(...)`, `CALL say_goodbye_and_save_results(...)`, etc.).
*   `read_vasp_files`: Specifically used by `io_routines` if the plane-wave code is VASP. (Implicit dependency via `io_routines` if `args%pw_code` dictates it).
*   `band_unfolding`: For the core unfolding logic (e.g., `CALL define_pckpts_to_be_checked(...)`, `CALL get_geom_unfolding_relations(...)`, `CALL allocate_UnfoldedQuantities(...)`, `CALL update_GUR_indices(...)`, `CALL select_coeffs_to_calc_spectral_weights(...)`, `CALL perform_unfolding(...)`, `CALL get_delta_Ns_for_output(...)`).

It also uses `omp_lib` for OpenMP parallelization directives.
