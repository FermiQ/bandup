# `band_unfolding_routines_mod.f90`

## Overview

This module, `band_unfolding`, authored by Paulo V C Medeiros from Linköping University, is central to the BandUP application. It encapsulates the core algorithms and subroutines required to perform the band unfolding process. This includes determining geometric relationships between primitive cell and supercell k-points, selecting relevant wavefunction coefficients, calculating spectral weights, spectral functions, and the unfolding density operator, as well as handling symmetry averaging and preparing data for output.

## Key Components

This module provides several public subroutines and functions that orchestrate the band unfolding calculations.

### Public Subroutines:

*   **`allocate_UnfoldedQuantities(delta_N, pckpts_to_be_checked)`**
    *   **Purpose**: Allocates memory for the `delta_N` variable of type `UnfoldedQuantities`.
    *   **Arguments**:
        *   `delta_N` (OUT, type `UnfoldedQuantities`): The variable to be allocated.
        *   `pckpts_to_be_checked` (IN, type `selected_pcbz_directions`): Provides the structure for allocation based on the number of selected pc-kpt directions and points.
    *   **Description**: Initializes the nested structure of `delta_N` to hold unfolding results for each selected pc-kpt path, each needed (symmetry-equivalent) direction, and each k-point along those paths.

*   **`allocate_UnfoldedQuantitiesForOutput(delta_N, pckpts_to_be_checked)`**
    *   **Purpose**: Allocates memory for variables of type `UnfoldedQuantitiesForOutput`.
    *   **Arguments**:
        *   `delta_N` (OUT, type `UnfoldedQuantitiesForOutput`): The variable to be allocated for storing output-formatted unfolded quantities.
        *   `pckpts_to_be_checked` (IN, type `selected_pcbz_directions`): Provides the structure for allocation.
    *   **Description**: Initializes the structure to hold unfolded quantities that are specifically prepared for output files.

*   **`update_GUR_indices(GUR, i_SCKPT, i_selec_pcbz_dir, i_needed_dirs, ipc_kpt)`**
    *   **Purpose**: Updates the current indices within the Geometric Unfolding Relations (GUR) structure.
    *   **Arguments**:
        *   `GUR` (INOUT, type `geom_unfolding_relations_for_each_SCKPT`): The GUR data structure.
        *   `i_SCKPT` (IN, integer): Current supercell k-point index.
        *   `i_selec_pcbz_dir` (IN, integer): Current selected primitive cell Brillouin zone direction index.
        *   `i_needed_dirs` (IN, integer): Current needed (symmetry-equivalent) primitive cell Brillouin zone direction index.
        *   `ipc_kpt` (IN, integer): Current primitive cell k-point index along the direction.
    *   **Description**: Sets the `current_index` component of the `GUR` variable to keep track of the specific k-point context during unfolding.

*   **`verify_commens(crystal_pc, crystal_SC, args)`**
    *   **Purpose**: Verifies if the primitive cell (`crystal_pc`) and supercell (`crystal_SC`) are commensurate.
    *   **Arguments**:
        *   `crystal_pc` (IN, type `crystal_3D`): Primitive cell structure.
        *   `crystal_SC` (IN, type `crystal_3D`): Supercell structure.
        *   `args` (IN, type `comm_line_args`): Command-line arguments, used here for `args%stop_if_not_commensurate`.
    *   **Description**: Calls `check_if_pc_and_SC_are_commensurate` (from `crystals_mod`) and prints a message. If they are not commensurate and `args%stop_if_not_commensurate` is true, the program stops.

*   **`get_geom_unfolding_relations(GUR, list_of_SCKPTS, pckpts_to_be_checked, input_crystal_pc, input_crystal_SC, verbose)`**
    *   **Purpose**: Determines the Geometric Unfolding Relations (GUR) between primitive cell and supercell k-points. This is a critical step that identifies which pc-kpoints fold into which SC K-points.
    *   **Arguments**:
        *   `GUR` (OUT, type `geom_unfolding_relations_for_each_SCKPT`): The main output, populated GUR structure.
        *   `list_of_SCKPTS` (IN, array of `vec3d`): List of supercell k-points.
        *   `pckpts_to_be_checked` (IN, type `selected_pcbz_directions`): The pc-kpoints to check for folding.
        *   `input_crystal_pc` (IN, type `crystal_3D`): Primitive cell crystal structure.
        *   `input_crystal_SC` (IN, type `crystal_3D`): Supercell crystal structure.
        *   `verbose` (IN, optional, logical): If true, prints more detailed output.
    *   **Description**: This is a wrapper for an internal routine `get_GUR_not_public`. It iteratively calls the internal routine with increasing tolerance if the GUR determination fails initially, up to `max_tol_for_vec_equality`. It populates the `GUR` structure, indicating for each pc-kpt whether it folds to an SC K-point, the corresponding SC K-point, and symmetrized k-vectors.

*   **`define_pckpts_to_be_checked(pckpts_to_be_checked, all_dirs_used_for_EBS_along_pcbz_dir, nkpts_selected_dirs)`**
    *   **Purpose**: Defines the specific discrete primitive cell k-points along the selected (and symmetry-equivalent) Brillouin zone directions that need to be analyzed.
    *   **Arguments**:
        *   `pckpts_to_be_checked` (OUT, type `selected_pcbz_directions`): Output structure containing the coordinates of all pc-kpoints to be checked.
        *   `all_dirs_used_for_EBS_along_pcbz_dir` (IN, array of `irr_bz_directions`): All irreducible BZ directions to be used for the EBS.
        *   `nkpts_selected_dirs` (IN, integer array): Number of k-points along each selected direction.
    *   **Description**: For each direction in `all_dirs_used_for_EBS_along_pcbz_dir`, it generates a line of `nkpts_selected_dirs` k-points between `kstart` and `kend` and stores their coordinates in `pckpts_to_be_checked`.

*   **`select_coeffs_to_calc_spectral_weights(selected_coeff_indices, wf, crystal_pc, GUR)`**
    *   **Purpose**: Selects the plane-wave coefficients from the supercell wavefunction (`wf`) that are relevant for calculating the spectral weight at the current primitive cell k-point (defined by `GUR%current_index`).
    *   **Arguments**:
        *   `selected_coeff_indices` (OUT, allocatable integer array): Indices of the supercell plane-wave coefficients that satisfy the selection criterion.
        *   `wf` (IN, type `pw_wavefunction`): Supercell wavefunction data.
        *   `crystal_pc` (IN, type `crystal_3D`): Primitive cell crystal structure (for its reciprocal lattice vectors).
        *   `GUR` (IN, type `geom_unfolding_relations_for_each_SCKPT`): Geometric Unfolding Relations, used to get the current pc-kpt and SC K-point.
    *   **Description**: A supercell plane-wave G_SC contributes to the spectral weight of a pc-kpt (k_pc) if (G_SC - (S*k_pc - K_sc)) is a reciprocal lattice vector of the primitive cell. (S*k_pc is the symmetrized pc-kpt, K_sc is the SC k-vector whose coefficients are being used). This routine identifies such G_SC vectors and stores their indices. It tries with increasing tolerance if no coefficients are found initially.

*   **`get_delta_Ns_for_output(delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS, delta_N, all_dirs_used_for_EBS_along_pcbz_dir, pckpts_to_be_checked, times)`**
    *   **Purpose**: Prepares the calculated unfolded quantities (`delta_N`) for output. This involves creating two versions: one strictly along user-selected directions and another that is symmetry-averaged for the EBS. It also handles the symmetry averaging of the unfolding density operator and spin projections if calculated.
    *   **Arguments**:
        *   `delta_N_only_selected_dirs` (OUT, type `UnfoldedQuantitiesForOutput`): Unfolded EBS strictly along user-requested directions.
        *   `delta_N_symm_avrgd_for_EBS` (OUT, type `UnfoldedQuantitiesForOutput`): Symmetry-averaged unfolded EBS.
        *   `delta_N` (IN, type `UnfoldedQuantities`): The raw unfolded quantities.
        *   `all_dirs_used_for_EBS_along_pcbz_dir` (IN, array of `irr_bz_directions`): Information about irreducible BZ directions and their weights for averaging.
        *   `pckpts_to_be_checked` (IN, type `selected_pcbz_directions`): Structure of pc-kpts.
        *   `times` (INOUT, optional, type `timekeeping`): For recording timing information.
    *   **Description**: It populates `delta_N_only_selected_dirs` by taking the data for the first `needed_dir`. It calculates `delta_N_symm_avrgd_for_EBS` by performing a weighted average of `delta_N` over all symmetry-equivalent `needed_dir` for each selected pc-kpt path. Similar averaging is done for the unfolding density operator (`rhos`) and spin components if `args%write_unf_dens_op` or `wf%is_spinor` is true.

*   **`perform_unfolding(delta_N, times, GUR, wf, selected_coeff_indices, energy_grid)`**
    *   **Purpose**: This is the core subroutine that performs the unfolding calculation for a single pc-kpt given the SC wavefunction and selected coefficients.
    *   **Arguments**:
        *   `delta_N` (INOUT, type `UnfoldedQuantities`): The data structure where results are stored (specifically, the sub-element corresponding to the current pc-kpt).
        *   `times` (INOUT, type `timekeeping`): For recording timing of sub-steps.
        *   `GUR` (IN, type `geom_unfolding_relations_for_each_SCKPT`): GUR structure providing context (current k-points).
        *   `wf` (IN, type `pw_wavefunction`): Supercell wavefunction.
        *   `selected_coeff_indices` (IN, integer array): Indices of relevant SC plane-wave coefficients.
        *   `energy_grid` (IN, real(dp) array): Energy grid for calculations.
    *   **Description**:
        1.  Calculates the spectral weight for each SC band using `spectral_weight_for_coeff` from the `selected_coeff_indices`.
        2.  Calculates `delta_N` (unfolded band character) vs. energy using `get_delta_Ns_for_EBS`.
        3.  Optionally (if `calc_spec_func_explicitly` is true), calculates the spectral function (`SF`) using `calc_spectral_function`.
        4.  If `args%write_unf_dens_op` is true or the wavefunction is spinorial (`n_WF_components==2`), it calculates the unfolding density operator (`rho`) for energies where `dN` is significant, using the internal `calc_rho` subroutine.
        5.  If the wavefunction is spinorial, it calculates the expectation values of Pauli matrices (`sigma`) and their projections (`spin_proj_perp`, `spin_proj_para`) using the calculated `rho` and SC Pauli matrix elements (obtained via `get_SC_pauli_matrix_elmnts`).

### Private Functions/Subroutines (Notable examples):

*   **`get_GUR_not_public(...)`**: The internal workhorse for `get_geom_unfolding_relations`. It iterates through SC K-points, pc-kpts, and symmetry operations to find folding conditions.
*   **`spectral_weight_for_coeff(coeff, selected_coeff_indices, add_elapsed_time_to)` (Function)**: Calculates the sum of `|coefficient|^2` for the selected indices for each band.
*   **`calc_spectral_function(SF_at_pckpt, energies, SC_calc_ener, spectral_weight, std_dev, add_elapsed_time_to)`**: Calculates the spectral function A(k,E) using Gaussian/Lorentzian broadening.
*   **`lambda(pc_ener, SC_ener, delta_e, std_dev)` (Function)**: Calculates the integral of a delta function (broadened) over an energy interval `delta_e`.
*   **`calc_delta_N_pckpt(delta_N_pckpt, energies, SC_calc_ener, spectral_weight, std_dev)`**: Calculates delta_N(k,E) by summing spectral weights times the `lambda` function.
*   **`get_delta_Ns_for_EBS(...)`**: A wrapper for `calc_delta_N_pckpt`.
*   **`calc_rho(rho, delta_N, pc_ener, delta_e, wf, selected_coeff_indices, std_dev, add_elapsed_time_to)`**: Calculates the unfolding-density operator.
*   **`get_SC_pauli_matrix_elmnts(...)`**: Calculates matrix elements of Pauli matrices between SC wavefunctions.
*   **`calc_spin_projections(...)`**: Calculates parallel and perpendicular projections of the spin vector.

## Important Variables/Constants

This module primarily contains routines. Module-level parameters are not its main feature, as most constants are imported from `constants_and_types_mod`. The logic relies heavily on the derived types defined in `constants_and_types_mod` passed as arguments.

## Usage Examples

The subroutines in this module are called by `main_BandUP.f90` to carry out the unfolding process. For example, `perform_unfolding` is called within the main loop of `BandUP_main`:

```fortran
! Conceptual call within main_BandUP.f90
IF (pckpt_folds) THEN
    ! ... (update GUR indices, read wavefunction if needed, select_coeffs_to_calc_spectral_weights) ...
    IF (ALLOCATED(selected_coeff_indices)) THEN
        CALL perform_unfolding( &
            delta_N, times, GUR, wf, &
            selected_coeff_indices, energy_grid &
        )
    END IF
END IF
```

## Dependencies and Interactions

This module has several dependencies:

*   **`constants_and_types`**: Extensively used for all derived types (e.g., `UnfoldedQuantities`, `pw_wavefunction`, `geom_unfolding_relations_for_each_SCKPT`, `crystal_3D`), numerical kinds (`dp`, `kind_cplx_coeffs`), and various default parameters/tolerances.
*   **`cla_wrappers`**: For accessing command-line arguments (e.g., `args%perform_unfold`, `args%write_unf_dens_op`, `args%no_symm_sckpts`).
*   **`lists_and_seqs`**: For list operations like `append`, `list_index`, and generating k-point lines (`kpts_line`).
*   **`math`**: For mathematical functions like `integral_delta_x_minus_x0`, `delta` (broadening function), `trace_ab`, `versor`, `cross`.
*   **`symmetry`**: For symmetry operations like `get_symm`, `get_star`, `pt_eqv_by_point_group_symop`.
*   **`crystals`**: For crystal-related checks like `vec_in_latt`, `check_if_pc_and_SC_are_commensurate`.
*   **`time`**: For `time_now()` to record elapsed times for various operations.
*   **`general_io` / `io_routines`**: For printing messages/warnings (e.g., `print_message_commens_test`, `print_message_pckpt_cannot_be_parsed`).
*   **`omp_lib`**: For OpenMP directives for parallelization.

The routines within this module are highly interconnected and also interact with routines from the modules listed above.
