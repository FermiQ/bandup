# `symmetry_mod.f90`

## Overview

This module, `symmetry`, authored by Paulo V C Medeiros from Linköping University, is responsible for handling symmetry operations and related concepts within the BandUP application. It leverages the `spglib` library for robust symmetry analysis. Key functionalities include finding primitive cells, determining crystal symmetry operations, generating the star of a k-point, finding irreducible k-points, and analyzing symmetry relations between primitive and supercells. These capabilities are crucial for correctly performing band unfolding and for averaging results over symmetric k-points.

The module's comment block incorrectly states its name as `math` and description related to mathematical functions. This is likely a copy-paste error from another module. The actual content is clearly focused on symmetry.

## Key Components

This module provides several public subroutines and functions for symmetry-related tasks.

### Public Subroutines and Functions:

*   **`get_prim_cell(crystal, symmetryze_and_refine_cell, symprec)`**
    *   **Purpose**: Finds the primitive cell for a given `crystal` structure.
    *   **Arguments**:
        *   `crystal` (INOUT, type `crystal_3D`): The input crystal structure. If it's not already a primitive cell, its `corresp_pc` pointer component will be allocated and will point to the found primitive cell. The `is_prim_cell` flag is updated.
        *   `symmetryze_and_refine_cell` (IN, optional, logical): If `.TRUE.`, `spg_refine_cell` is called to refine the input cell before finding the primitive cell. Defaults to `.FALSE.`.
        *   `symprec` (IN, optional, real(dp)): Symmetry precision for `spglib` calls. Defaults to `default_symprec`.
    *   **Description**: Uses `spg_find_primitive` from `spglib` to determine the primitive cell. If the input `crystal` is already primitive, `crystal%corresp_pc` points to a copy of itself (though the logic seems to imply it would just set `is_prim_cell` and not necessarily create a distinct `corresp_pc` unless refinement happens or it wasn't primitive). If a smaller primitive cell is found, `crystal%corresp_pc` is populated with the new primitive cell structure.

*   **`get_symm(crystal, use_pc_to_get_symm, symprec)`**
    *   **Purpose**: Determines the symmetry operations for a given `crystal`.
    *   **Arguments**:
        *   `crystal` (INOUT, type `crystal_3D`): The crystal structure. Its `nsym` and `symops` components are populated.
        *   `use_pc_to_get_symm` (IN, optional, logical): If `.TRUE.`, the symmetry is determined from the crystal's primitive cell (if `crystal` is not already primitive). Defaults to `.FALSE.`.
        *   `symprec` (IN, optional, real(dp)): Symmetry precision. Defaults to `default_symprec`.
    *   **Description**: Uses `spg_get_dataset` from `spglib` to get the symmetry operations. The rotation and translation parts of each symmetry operation are stored in `crystal%symops` in both fractional and Cartesian coordinates. It also populates `crystal%international_symb` and `crystal%schoenflies` symbol.

*   **`pt_eqv_by_point_group_symop(point, symops, isym, fractional_coords, invert_symop)` (Function)**
    *   **Purpose**: Applies a specific point group symmetry operation to a `point`.
    *   **Arguments**:
        *   `point` (IN, real(dp), dim(3)): The input point (k-point or real-space vector).
        *   `symops` (IN, array of type `symmetry_operation`): Array of symmetry operations for the crystal.
        *   `isym` (IN, integer): Index of the symmetry operation in `symops` to apply.
        *   `fractional_coords` (IN, optional, logical): If `.TRUE.`, uses fractional coordinate rotation matrix; otherwise, uses Cartesian. Defaults to `.FALSE.` (Cartesian).
        *   `invert_symop` (IN, optional, logical): If `.TRUE.`, applies the inverse of the rotation matrix. Defaults to `.FALSE.`.
    *   **Returns**: (real(dp), dim(3)) The transformed point.
    *   **Description**: Multiplies the input `point` by the rotational part of the specified symmetry operation (`symops(isym)`). Note: This function only applies the rotational part, not the translational part, making it suitable for point group operations on k-vectors or vectors relative to an origin.

*   **`get_star(star_of_pt, list_of_all_generated_points, points, crystal, tol_for_vec_equality, symprec, reduce_to_bz, try_to_reduce_to_a_pc)`**
    *   **Purpose**: Generates the star of a set of `points` (k-points). The star is the set of all unique points obtained by applying all symmetry operations of the `crystal` to each input point.
    *   **Arguments**:
        *   `star_of_pt` (OUT, optional, allocatable array of type `star`): For each input point, stores its star (list of unique equivalent points and the generating symmetry operation index).
        *   `list_of_all_generated_points` (OUT, optional, allocatable array of type `vec3d`): A flat list of all unique points generated from all input points.
        *   `points` (IN, array of type `vec3d`): List of input points (Cartesian coordinates).
        *   `crystal` (IN, type `crystal_3D`): The crystal structure providing symmetry operations.
        *   `tol_for_vec_equality` (IN, optional, real(dp)): Tolerance for comparing if two generated points are the same. Defaults to `default_tol_for_vec_equality`.
        *   `symprec` (IN, optional, real(dp)): Symmetry precision.
        *   `reduce_to_bz` (IN, optional, logical): If `.TRUE.`, input points are first reduced to the 1st Brillouin Zone. Defaults to `.FALSE.`.
        *   `try_to_reduce_to_a_pc` (IN, optional, logical): If `.TRUE.`, symmetry operations are obtained from the primitive cell if `crystal` is a supercell.
    *   **Description**: For each input `point`, it applies all symmetry operations (from `crystal%symops`, calling `get_symm` if needed). It then identifies the set of unique generated points, removing duplicates within the specified tolerance.

*   **`get_irr_SC_kpts(n_irr_kpts, irr_kpts_list, irr_kpts_list_frac_coords, kpts_list, crystal, args, reduce_to_bz, symprec)`**
    *   **Purpose**: Finds the set of irreducible k-points from a given `kpts_list` for a supercell (`crystal`).
    *   **Arguments**:
        *   `n_irr_kpts` (OUT, integer): Number of irreducible k-points found.
        *   `irr_kpts_list` (OUT, allocatable array of type `vec3d`): List of irreducible k-points (Cartesian coordinates).
        *   `irr_kpts_list_frac_coords` (OUT, allocatable array of type `vec3d`): List of irreducible k-points (fractional coordinates).
        *   `kpts_list` (IN, array of type `vec3d`): Input list of SC k-points.
        *   `crystal` (IN, type `crystal_3D`): The supercell structure.
        *   `args` (IN, type `comm_line_args`): Command-line arguments (used for `args%no_symm_sckpts`).
        *   `reduce_to_bz` (IN, optional, logical): If `.TRUE.`, k-points are first reduced to the 1st BZ.
        *   `symprec` (IN, optional, real(dp)): Symmetry precision.
    *   **Description**: First, optionally reduces all k-points in `kpts_list` to the first BZ. Then, it uses `get_star` to find symmetry equivalences. It iterates through the `kpts_list` and identifies k-points that are not equivalent to any preceding k-point under symmetry, thus forming the irreducible set.

*   **`analise_symm_pc_SC(crystal_pc, crystal_SC, symprec, verbose)`**
    *   **Purpose**: Analyzes and prints symmetry information for the primitive cell (`crystal_pc`) and supercell (`crystal_SC`).
    *   **Arguments**:
        *   `crystal_pc` (INOUT, optional, type `crystal_3D`): Primitive cell structure. Symmetry info will be populated.
        *   `crystal_SC` (INOUT, optional, type `crystal_3D`): Supercell structure. Symmetry info will be populated.
        *   `symprec` (IN, optional, real(dp)): Symmetry precision.
        *   `verbose` (IN, optional, logical): If `.TRUE.` (default), prints detailed symmetry information.
    *   **Description**: Calls `get_symm` for the provided cell(s) and prints the number of symmetry operations found, the Schoenflies symbol, and the International symbol. It includes notes if the symmetry was determined using an associated primitive cell.

*   **`get_pcbz_dirs_2b_used_for_EBS(all_dirs_used_for_EBS_along_pcbz_dir, input_crystal_pc, input_crystal_SC, k_starts, k_ends, args, verbose, symprec)`**
    *   **Purpose**: Determines all unique primitive cell Brillouin Zone (PC BZ) directions that need to be calculated for a symmetry-averaged Electronic Band Structure (EBS), based on user-selected k-point paths and the symmetries of both the PC and SC.
    *   **Arguments**:
        *   `all_dirs_used_for_EBS_along_pcbz_dir` (OUT, allocatable array of `irr_bz_directions`): For each input k-path segment, this stores a list of irreducible directions (and their weights) that cover all symmetry-equivalent parts.
        *   `input_crystal_pc` (IN, type `crystal_3D`): Primitive cell structure.
        *   `input_crystal_SC` (IN, type `crystal_3D`): Supercell structure.
        *   `k_starts`, `k_ends` (IN, real(dp) array, dim(:,3)): Start and end points of user-selected k-paths in PC BZ.
        *   `args` (IN, type `comm_line_args`): Command-line arguments (used for `args%no_symm_avg`).
        *   `verbose`, `symprec` (IN, optional): Control verbosity and symmetry precision.
    *   **Description**:
        1.  If `args%no_symm_avg` is true, each input path is treated as a single irreducible direction with weight 1.
        2.  Otherwise, for each input k-path segment (`k_starts(idir,:)` to `k_ends(idir,:)`):
            a.  Finds all directions equivalent to this segment under SC symmetry (`get_eqv_bz_dirs` with `crystal_SC`).
            b.  Finds all directions equivalent to this segment under PC symmetry (`get_eqv_bz_dirs` with `crystal_pc`).
            c.  Identifies "complementary" PC BZ directions: those equivalent under PC symmetry but *not* under SC symmetry (using `get_compl_pcbz_direcs`).
            d.  Finds the irreducible set among these complementary directions using SC symmetry (`get_irr_bz_directions` with `crystal_SC`).
            e.  The final set of directions for this input segment includes the original segment (weighted by non-complementary PC symmetries) and the irreducible complementary directions (weighted by their respective degeneracies). Ensures total weight sums to 1.

### Private Functions/Subroutines (Notable examples):

*   **`get_compl_pcbz_direcs(...)`**: Identifies directions equivalent under PC symmetry but not SC symmetry.
*   **`get_irr_bz_directions(...)`**: Finds irreducible set of BZ directions from a list, given a crystal symmetry.
*   **`get_eqv_bz_dirs(...)`**: Finds all BZ directions equivalent to a given k-start/k-end pair under a crystal's symmetry.
*   **`points_are_eqv(v1, v2, symops, symprec)` (Function)**: Checks if two points `v1` and `v2` are equivalent under any of the `symops`.
*   **`direcs_are_eqv(start1, end1, start2, end2, symops, symprec)` (Function)**: Checks if two k-path segments (defined by start/end points) are equivalent under `symops`. They must have the same length and their start/end points must be symmetrically equivalent.

## Important Variables/Constants

This module relies heavily on:
*   Derived types from `constants_and_types_mod`, especially `crystal_3D`, `symmetry_operation`, `star`, `bz_direction`, `eqv_bz_directions`, `irr_bz_directions`, `vec3d`.
*   Constants like `default_symprec` from `constants_and_types_mod`.
*   The `spglib_f08` module for direct calls to `spglib` functions.

## Usage Examples

The routines in this module are primarily called by `main_BandUP.f90` to set up and guide the unfolding calculations.

```fortran
! In main_BandUP.f90
CALL analise_symm_pc_SC(crystal_pc, crystal_SC)

CALL get_pcbz_dirs_2b_used_for_EBS( &
    all_dirs_used_for_EBS_along_pcbz_dir, &
    crystal_pc, crystal_SC, &
    k_starts, k_ends, args &
)

! Later, when determining Geometric Unfolding Relations (in band_unfolding_mod, called by main)
! get_star is used to find equivalent SC K-points.
! CALL get_star(star_of_pt=SKPTS_eqv_to_SKPT, points=list_of_SCKPTS, crystal=crystal_SC, ...)
```

## Dependencies and Interactions

*   **`constants_and_types`**: For all crystal and symmetry related derived types, default symmetry precision.
*   **`math`**: For `inverse_of_3x3_matrix`, `same_vector`, `norm`, `coords_cart_vec_in_new_basis`.
*   **`spglib_f08`**: This is a critical dependency, providing the underlying symmetry-finding capabilities (e.g., `spg_get_dataset`, `spg_find_primitive`, `spg_refine_cell`).
*   **`crystals`**: For `create_crystal` (used in `get_prim_cell`) and `reduce_point_to_bz` (used in `get_star` and `get_irr_SC_kpts`).
*   **`cla_wrappers`**: (Indirectly) `args` (type `comm_line_args`) is passed to some routines.

This module is essential for correctly handling symmetries, which is vital for efficient and accurate band unfolding, especially when averaging or comparing results across different parts of the Brillouin zone.
