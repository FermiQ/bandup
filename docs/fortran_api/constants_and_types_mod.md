# `constants_and_types_mod.f90`

## Overview

This module provides most constants and derived data types used throughout the BandUP application. It serves as a central repository for fundamental definitions that are critical for the core logic of the band unfolding calculations and data representation. The module is authored by Paulo V C Medeiros, Linköping University.

## Key Components

This module defines a collection of parameters (constants) and derived types.

### Parameters (Constants)

The module makes the following parameters publicly available. These parameters define various physical constants, numerical precision kinds, default values for tolerances, and control flags for the application's behavior.

**Integer Kinds:**
*   `sp`: Integer kind for single precision real numbers.
*   `dp`: Integer kind for double precision real numbers.
*   `long_int_kind`: Integer kind for long integers.
*   `kind_cplx_coeffs`: Integer kind for complex coefficients (typically single precision).
*   `qe_dp`: Integer kind for Quantum ESPRESSO's double precision.
*   `str_len`: Integer parameter defining a standard string length (256).

**Real Constants:**
*   `pi`: Value of Pi (π).
*   `twopi`: Value of 2*Pi (2π).
*   `min_dk`: Minimum difference between k-points.
*   `default_tol_for_vec_equality`: Default tolerance for vector equality checks (1E-5).
*   `max_tol_for_vec_equality`: Maximum tolerance for vector equality checks (1E-3).
*   `default_tol_for_int_commens_test`: Default tolerance for integer commensurability tests (1E-5).
*   `default_symprec`: Default symmetry precision (1E-5).
*   `two_m_over_hbar_sqrd`: Constant 2m/ħ² in eV⁻¹Å⁻² (0.262465831).
*   `identity_3D`: 3x3 Identity matrix.
*   `Hartree`: Hartree energy in eV (27.2113834 eV).
*   `Ry`: Rydberg energy in eV (13.60569172 eV).
*   `bohr`: Bohr radius in Ångströms (0.52917721092 Å).

**Logical Flags:**
*   `calc_spec_func_explicitly`: Flag to control explicit calculation of spectral function (.FALSE.).
*   `stop_if_pckpt_cannot_be_parsed`: Flag to stop execution if a pc-kpt cannot be parsed (.TRUE.).
*   `stop_if_GUR_fails`: Flag to stop execution if Geometric Unfolding Relations (GUR) fail (.TRUE.).
*   `get_all_kpts_needed_for_EBS_averaging`: Flag to get all k-points needed for EBS averaging (.TRUE.).
*   `print_GUR_pre_unfolding_utility`: Flag to print GUR before unfolding utility (.FALSE.).
*   `renormalize_wf`: Flag to renormalize wavefunctions (.TRUE.).

**Character Strings:**
*   `package_version`: Character string holding the package version and commit date (e.g., "v3.0.0-beta.6, (unknown commit date)").

### Derived Types

The module defines several derived types to structure the data used in BandUP.

*   **`comm_line_args`**: Holds command-line arguments passed to the program.
    *   `WF_file` (character): Wavefunction file path.
    *   `input_file_prim_cell` (character): Primitive cell input file.
    *   `input_file_supercell` (character): Supercell input file.
    *   `input_file_pc_kpts` (character): PC k-points input file.
    *   `input_file_energies` (character): Energies input file.
    *   `out_file_SC_kpts` (character): Output file for SC k-points.
    *   `output_file_symm_averaged_EBS` (character): Output file for symmetry-averaged EBS.
    *   `output_file_only_user_selec_direcs` (character): Output file for user-selected directions EBS.
    *   `output_file_symm_averaged_unf_dens_op` (character): Output file for symmetry-averaged unfolding density operator.
    *   `output_file_only_user_selec_direcs_unf_dens_op` (character): Output file for user-selected directions unfolding density operator.
    *   `pw_code` (character): Plane-wave code identifier.
    *   `qe_outdir` (character): Quantum ESPRESSO output directory.
    *   `qe_prefix` (character): Quantum ESPRESSO prefix.
    *   `abinit_files_file` (character): File containing ABINIT file paths.
    *   `castep_seed` (character): CASTEP seed name.
    *   `spin_channel` (integer): Spin channel selection.
    *   `n_sckpts_to_skip` (integer): Number of SC k-points to skip.
    *   `saxis` (real(dp), dim(3)): Spin quantization axis.
    *   `normal_to_proj_plane` (real(dp), dim(3)): Normal vector to the projection plane for spin.
    *   `origin_for_spin_proj_cartesian` (real(dp), dim(3)): Cartesian origin for spin projection.
    *   `origin_for_spin_proj_rec` (real(dp), dim(3)): Reciprocal origin for spin projection.
    *   `stop_if_not_commensurate` (logical): Flag to stop if cells are not commensurate.
    *   `write_attempted_pc_corresp_to_input_pc` (logical): Flag to write attempted PC corresponding to input PC.
    *   `write_attempted_pc_corresp_to_SC` (logical): Flag to write attempted PC corresponding to SC.
    *   `no_symm_avg` (logical): Flag to disable symmetry averaging.
    *   `no_symm_sckpts` (logical): Flag to disable SC k-point symmetrization.
    *   `perform_unfold` (logical): Flag to control whether unfolding is performed.
    *   `origin_for_spin_proj_passed_in_rec` (logical): Flag indicating if spin projection origin is in reciprocal coordinates.
    *   `continue_if_npw_smaller_than_expected` (logical): Flag to continue if number of plane waves is smaller than expected.
    *   `write_unf_dens_op` (logical): Flag to control writing of the unfolding density operator.

*   **`timekeeping`**: Stores timing information for various parts of the calculation.
    *   `start`, `end` (real(dp)): Overall start and end times.
    *   `read_wf` (real(dp)): Time spent reading wavefunctions.
    *   `calc_spec_weights` (real(dp)): Time for spectral weights calculation.
    *   `calc_SF` (real(dp)): Time for spectral function calculation.
    *   `calc_dN` (real(dp)): Time for delta_N calculation.
    *   `calc_rho` (real(dp)): Time for unfolding density operator calculation.
    *   `calc_pauli_vec` (real(dp)): Time for Pauli vector calculation.
    *   `calc_pauli_vec_projs` (real(dp)): Time for Pauli vector projections.
    *   `write_dN_files` (real(dp)): Time for writing delta_N files.
    *   `write_unf_dens_op_files` (real(dp)): Time for writing unfolding density operator files.

*   **`vec3d`**: Represents a 3D vector with double precision real coordinates.
    *   `coord` (real(dp), dim(3)): 3D coordinates.

*   **`vec3d_int`**: Represents a 3D vector with integer coordinates.
    *   `coord` (integer, dim(3)): 3D integer coordinates.

*   **`MatrixIndices`**: Represents indices for a matrix element.
    *   `indices` (integer, dim(2)): Pair of indices (row, col).

*   **`pw_wavefunction`**: Stores information about the plane-wave wavefunction.
    *   `i_spin` (integer): Spin index.
    *   `n_pw` (integer): Number of plane waves.
    *   `n_spin` (integer): Number of spin components.
    *   `n_bands` (integer): Number of bands.
    *   `n_spinor` (integer): Number of spinor components (1 or 2).
    *   `n_bands_up`, `n_bands_down` (integer): Number of up/down bands (for collinear spin).
    *   `encut` (real(dp)): Energy cutoff.
    *   `Vcell` (real(dp)): Cell volume.
    *   `e_fermi`, `e_fermi_up`, `e_fermi_down` (real(dp)): Fermi energy / up/down Fermi energies.
    *   `kpt_frac_coords` (real(dp), dim(3)): K-point fractional coordinates.
    *   `kpt_cart_coords` (real(dp), dim(3)): K-point Cartesian coordinates.
    *   `A_matrix` (real(dp), dim(3,3)): Direct lattice vectors (columns).
    *   `B_matrix` (real(dp), dim(3,3)): Reciprocal lattice vectors (columns).
    *   `band_energies` (real(dp), allocatable, dim(:)): Band energies.
    *   `band_occupations` (real(dp), allocatable, dim(:)): Band occupations.
    *   `G_cart` (type(vec3d), allocatable, dim(:)): Cartesian coordinates of G-vectors for plane waves.
    *   `G_frac` (type(vec3d_int), allocatable, dim(:)): Fractional coordinates of G-vectors.
    *   `pw_coeffs` (complex(kind_cplx_coeffs), allocatable, dim(:,:,:)): Plane-wave coefficients.
    *   `is_spinor` (logical): True if it's a spinor wavefunction.
    *   `two_efs` (logical): True if there are separate up/down Fermi energies.

*   **`symmetry_operation`**: Defines a symmetry operation.
    *   `translation_fractional_coords` (integer, dim(3)): Fractional translation vector.
    *   `rotation_fractional_coords` (integer, dim(3,3)): Rotation matrix in fractional coordinates.
    *   `translation_cartesian_coords` (real(dp), dim(3)): Cartesian translation vector.
    *   `rotation_cartesian_coords` (real, dim(3,3)): Rotation matrix in Cartesian coordinates.
    *   `basis` (real(dp), dim(3,3)): Basis vectors for the operation.

*   **`crystal_3D`**: Represents a 3D crystal structure.
    *   `description` (character): Description of the crystal.
    *   `latt_vecs` (real(dp), dim(3,3)): Lattice vectors (columns).
    *   `rec_latt_vecs` (real(dp), dim(3,3)): Reciprocal lattice vectors (columns).
    *   `vol` (real(dp)): Cell volume.
    *   `rec_latt_vol` (real(dp)): Reciprocal cell volume.
    *   `coords_basis_atoms`, `fractional_coords_basis_atoms` (real(dp), allocatable, dim(:,:)): Cartesian/Fractional coordinates of basis atoms.
    *   `unconstrained_dof_basis_atoms` (logical, allocatable, dim(:,:)): Flags for unconstrained degrees of freedom of basis atoms.
    *   `atomic_symbols_basis_atoms` (character(len=3), allocatable, dim(:)): Atomic symbols.
    *   `integer_types_basis_atoms` (integer, allocatable, dim(:)): Integer types for atoms.
    *   `schoenflies` (character(len=10)): Schoenflies symbol.
    *   `international_symb` (character(len=11)): International symbol.
    *   `space_group_num` (integer): Space group number.
    *   `nsym` (integer): Number of symmetry operations.
    *   `symops` (type(symmetry_operation), allocatable, dim(:)): Array of symmetry operations.
    *   `is_prim_cell` (logical): True if this represents a primitive cell.
    *   `corresp_pc` (type(crystal_3D), pointer): Pointer to the corresponding primitive cell.

*   **`star_point_properties`**: Properties of a point in a star (set of symmetry-equivalent k-points).
    *   `coord` (real(dp), dim(3)): Coordinates of the k-point.
    *   `symop` (integer): Index of the symmetry operation that generates this point.

*   **`star`**: Represents a star of k-points.
    *   `neqv` (integer): Number of equivalent points in the star.
    *   `eqv_pt` (type(star_point_properties), allocatable, dim(:)): Array of equivalent points.

*   **`bz_direction`**: Defines a direction in the Brillouin Zone for band structure paths.
    *   `kstart` (real(dp), dim(3)): Starting k-point of the direction.
    *   `kend` (real(dp), dim(3)): Ending k-point of the direction.
    *   `weight` (real(dp)): Weight of this direction (for averaging).
    *   `neqv` (integer): Number of equivalent directions.

*   **`eqv_bz_directions`**: A set of equivalent Brillouin Zone directions.
    *   `neqv` (integer): Number of equivalent directions in this set.
    *   `eqv_dir` (type(bz_direction), allocatable, dim(:)): Array of equivalent directions.

*   **`irr_bz_directions`**: A set of irreducible Brillouin Zone directions.
    *   `irr_dir` (type(bz_direction), allocatable, dim(:)): Array of irreducible directions.
    *   `neqv` (integer): Number of equivalent directions this irreducible set represents in PC.
    *   `neqv_SCBZ` (integer): Number of equivalent directions this irreducible set represents in SCBZ.
    *   `ncompl_dirs` (integer): Number of complementary directions.
    *   `n_irr_compl_dirs` (integer): Number of irreducible complementary directions.

*   **`trial_folding_pckpt`**: Information about a trial pc-kpt for unfolding.
    *   `coords_actual_unfolding_K` (real(dp), dim(3)): Actual SC K-vector used for unfolding this pc-kpt.
    *   `coords` (real(dp), dim(3)): Coordinates of the pc-kpt.
    *   `coords_SCKPT_used_for_coeffs` (real(dp), dim(3)): SC K-vector whose coefficients are used.
    *   `Scoords` (real(dp), dim(3)): Symmetrized coordinates of the pc-kpt.
    *   `Sfolding_vec` (real(dp), dim(3)): Symmetrized folding vector (S*k_pc - K_sc).
    *   `Sorigin_for_spin_proj` (real(dp), dim(3)): Symmetrized origin for spin projection.
    *   `folds` (logical): True if this pc-kpt folds to the current SC K-point.

*   **`list_of_trial_folding_pckpts`**: A list of `trial_folding_pckpt`.
    *   `pckpt` (type(trial_folding_pckpt), allocatable, dim(:)): Array of trial pc-kpts.

*   **`needed_pcbz_dirs_for_EBS`**: Holds lists of trial pc-kpts for needed PC BZ directions for EBS.
    *   `needed_dir` (type(list_of_trial_folding_pckpts), allocatable, dim(:)): Array of lists, one for each needed direction.

*   **`selected_pcbz_directions`**: Holds `needed_pcbz_dirs_for_EBS` for user-selected PC BZ directions.
    *   `selec_pcbz_dir` (type(needed_pcbz_dirs_for_EBS), allocatable, dim(:)): Array, one for each selected PC BZ direction.

*   **`GUR_indices`**: Indices for navigating Geometric Unfolding Relations (GUR).
    *   `i_SCKPT` (integer): Index for SC K-point.
    *   `i_selec_pcbz_dir` (integer): Index for selected PC BZ direction.
    *   `i_needed_dirs` (integer): Index for needed (symmetry-equivalent to selected) PC BZ direction.
    *   `ipc_kpt` (integer): Index for pc-kpt along a direction.

*   **`geom_unfolding_relations_for_each_SCKPT`**: Stores the geometric unfolding relations.
    *   `SCKPT` (type(selected_pcbz_directions), allocatable, dim(:)): GUR for each SC K-point.
    *   `SCKPT_used_for_unfolding` (logical, allocatable, dim(:)): Flag indicating if an SC K-point is used for unfolding.
    *   `n_pckpts` (integer): Total number of pc-kpts considered.
    *   `n_folding_pckpts` (integer): Number of pc-kpts that actually fold.
    *   `list_of_SCKPTS` (type(vec3d), allocatable, dim(:)): List of SC K-points.
    *   `current_index` (type(GUR_indices)): Current indices being processed.
    *   `B_matrix_SC` (real(dp), dim(3,3)): Reciprocal lattice of the supercell.
    *   `b_matrix_pc` (real(dp), dim(3,3)): Reciprocal lattice of the primitive cell.

*   **`UnfoldDensityOpContainer`**: Container for the unfolding density operator at a specific energy.
    *   `rho` (complex(kind_cplx_coeffs), allocatable, dim(:)): Stores non-zero elements of the unfolding density operator (sparse).
    *   `band_indices` (type(MatrixIndices), allocatable, dim(:)): Corresponding (SC band m1, SC band m2) indices for elements in `rho`.
    *   `iener_in_full_pc_egrid` (integer): Index of the energy point in the full PC energy grid.
    *   `nbands` (integer): Number of SC bands.

*   **`unfolded_quantities_for_given_pckpt`**: Holds all unfolded quantities for a given pc-kpt over an energy grid.
    *   `dN` (real(dp), allocatable, dim(:)): Unfolded band character (delta_N) at each energy point.
    *   `SF` (real(dp), allocatable, dim(:)): Spectral function at each energy point.
    *   `spin_proj_perp` (real(dp), allocatable, dim(:)): Perpendicular spin projection.
    *   `spin_proj_para` (real(dp), allocatable, dim(:)): Parallel spin projection.
    *   `sigma` (real(dp), allocatable, dim(:,:)): Expected values of Pauli matrices (sigma_x, sigma_y, sigma_z) at each energy point.
    *   `rhos` (type(UnfoldDensityOpContainer), allocatable, dim(:)): Array of unfolding density operators, one for each energy where dN > 0.

*   **`list_of_pckpts_for_unfolded_quantities`**: A list of `unfolded_quantities_for_given_pckpt`.
    *   `pckpt` (type(unfolded_quantities_for_given_pckpt), allocatable, dim(:)): Array of unfolded quantities, one for each pc-kpt in a direction.

*   **`UnfoldedQuantities_info_for_needed_pcbz_dirs`**: Holds lists of unfolded quantities for needed PC BZ directions.
    *   `needed_dir` (type(list_of_pckpts_for_unfolded_quantities), allocatable, dim(:)): Array, one for each needed direction.

*   **`UnfoldedQuantities`**: Main container for all unfolded quantities. Structure: `selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%PROPERTY(iener)`.
    *   `n_SC_bands` (integer): Number of SC bands.
    *   `selec_pcbz_dir` (type(UnfoldedQuantities_info_for_needed_pcbz_dirs), allocatable, dim(:)): Array, one for each selected PC BZ direction.

*   **`UnfoldedQuantitiesForOutput`**: Container for unfolded quantities formatted for output.
    *   `n_SC_bands` (integer): Number of SC bands.
    *   `pcbz_dir` (type(list_of_pckpts_for_unfolded_quantities), allocatable, dim(:)): Array, one for each PC BZ direction to be output.


## Important Variables/Constants

The "Key Components" section above lists all publicly accessible parameters and derived types along with their fields. These are the primary variables and constants provided by this module. The parameters are defined with default values that influence the behavior of the unfolding calculations, such as tolerances for comparisons and flags for optional computations. The derived types structure complex data like wavefunctions, crystal geometries, symmetry operations, and the results of the unfolding process.

## Usage Examples

This module is not typically used directly by end-users. Instead, its parameters and types are imported and utilized by other modules within the BandUP application. For example:

```fortran
! In another module/program
USE constants_and_types, ONLY: dp, crystal_3D, default_symprec

REAL(KIND=dp) :: lattice_parameter
TYPE(crystal_3D) :: my_crystal

lattice_parameter = 5.43_dp ! Angstroms
! ... operations using my_crystal and default_symprec ...
```

## Dependencies and Interactions

This module has no explicit `USE` statements pointing to other custom modules within BandUP, as it provides the foundational constants and types. However, it is **depended upon by almost all other modules** in the BandUP `src` directory. These other modules `USE constants_and_types_mod` to access precision kinds (like `dp`, `sp`), physical constants (like `pi`, `Hartree`), logical flags, and the various derived types for representing and manipulating data.
