# `crystals_mod.f90`

## Overview

This module, `crystals`, authored by Paulo V C Medeiros from Linköping University, provides derived types and routines to define and manipulate crystal structures within the BandUP application. It includes functionalities for creating crystal objects, calculating reciprocal lattice vectors, checking for commensurability between supercell and primitive cell, and operations related to k-points within the Brillouin zone, such as checking if a vector is a lattice vector or reducing a k-point to the first Brillouin zone.

## Key Components

This module defines several public subroutines and functions for crystal-related operations. The primary data structure for representing crystals, `crystal_3D`, is imported from `constants_and_types_mod`.

### Public Subroutines and Functions:

*   **`check_if_pc_and_SC_are_commensurate(commensurate, M, crystal_pc, crystal_SC, tol)`**
    *   **Purpose**: Checks if a primitive cell (`crystal_pc`) and a supercell (`crystal_SC`) are commensurate.
    *   **Arguments**:
        *   `commensurate` (OUT, logical): Set to `.TRUE.` if cells are commensurate, `.FALSE.` otherwise.
        *   `M` (OUT, real(dp), dim(3,3)): The transformation matrix such that `A_i = sum(M_ij * a_j)`, where `A` are SC lattice vectors and `a` are PC lattice vectors.
        *   `crystal_pc` (IN, type `crystal_3D`): The primitive cell structure.
        *   `crystal_SC` (IN, type `crystal_3D`): The supercell structure.
        *   `tol` (IN, optional, real(dp)): Tolerance for checking if elements of `M` are integers. Defaults to `default_tol_for_int_commens_test`.
    *   **Description**: Calculates `M` using the reciprocal lattice vectors of both cells. Then, it checks if the elements of `M` are integers within the specified tolerance.

*   **`create_crystal(crystal, description, latt_vecs, coords_basis_atoms, atomic_symbols_basis_atoms, unconstrained_dof_basis_atoms)`**
    *   **Purpose**: Creates and initializes a `crystal_3D` object.
    *   **Arguments**:
        *   `crystal` (OUT, type `crystal_3D`): The crystal structure object to be created.
        *   `description` (IN, optional, character): A text description for the crystal.
        *   `latt_vecs` (IN, real(dp), dim(3,3)): The lattice vectors as columns (e.g., `latt_vecs(:,1)` is the first lattice vector).
        *   `coords_basis_atoms` (IN, optional, real(dp), dim(:,3)): Cartesian coordinates of basis atoms. If not provided, a single atom at the origin is assumed.
        *   `atomic_symbols_basis_atoms` (IN, optional, character(len=3), dim(:)): Symbols for each atom in the basis.
        *   `unconstrained_dof_basis_atoms` (IN, optional, logical, dim(:,3)): Flags indicating unconstrained degrees of freedom for each atom.
    *   **Description**: Initializes the `crystal_3D` type instance. It sets the lattice vectors, calculates the cell volume, computes reciprocal lattice vectors (using `get_rec_latt`), stores atomic coordinates (both Cartesian and fractional), atomic symbols, and assigns integer types to atoms based on unique symbols.

*   **`get_rec_latt(latt, rec_latt, rec_latt_vol)`**
    *   **Purpose**: Calculates the reciprocal lattice vectors.
    *   **Arguments**:
        *   `latt` (IN, real(dp), dim(3,3)): Direct lattice vectors (as columns or rows, usage implies columns based on cross product order).
        *   `rec_latt` (OUT, real(dp), dim(3,3)): Calculated reciprocal lattice vectors (as rows: `b_i = 2pi * (a_j x a_k) / (a_i . (a_j x a_k))`).
        *   `rec_latt_vol` (OUT, optional, real(dp)): Volume of the reciprocal lattice unit cell.
    *   **Description**: Computes the standard reciprocal lattice vectors `b1, b2, b3` from the direct lattice vectors `a1, a2, a3`.

*   **`vec_in_latt(vec, latt, tolerance)` (Function)**
    *   **Purpose**: Checks if a given vector `vec` is a lattice vector of the lattice defined by `latt`.
    *   **Arguments**:
        *   `vec` (IN, real(dp), dim(3)): The vector to check (Cartesian coordinates).
        *   `latt` (IN, real(dp), dim(3,3)): The direct lattice vectors.
        *   `tolerance` (IN, optional, real(dp)): Tolerance for checking if fractional coordinates are integers. Defaults to `default_tol_for_vec_equality`.
    *   **Returns**: (logical) `.TRUE.` if `vec` is a lattice vector, `.FALSE.` otherwise.
    *   **Description**: Converts `vec` to fractional coordinates in the basis of `latt`. If these fractional coordinates are integers (within tolerance), the vector is a lattice vector. It checks combinations of `vec` with `+/- latt_i` to handle points near cell boundaries.

*   **`reduce_point_to_bz(point, rec_latt, point_reduced_to_bz, max_index_lin_comb_RL_vecs)` (Recursive Subroutine)**
    *   **Purpose**: Reduces a given k-point `point` to its equivalent in the first Brillouin Zone (BZ) defined by `rec_latt`.
    *   **Arguments**:
        *   `point` (IN, real(dp), dim(3)): The k-point to reduce (Cartesian coordinates).
        *   `rec_latt` (IN, real(dp), dim(3,3)): The reciprocal lattice vectors.
        *   `point_reduced_to_bz` (OUT, real(dp), dim(3)): The k-point reduced to the first BZ (Cartesian coordinates).
        *   `max_index_lin_comb_RL_vecs` (IN, optional, integer): Controls the search depth for the nearest reciprocal lattice vector. The routine recursively increases this if reduction fails.
    *   **Description**: Finds the reciprocal lattice vector `G` such that `point - G` lies within the first BZ. This is done by finding the `G` that is closest to `point`. The subroutine calls `point_is_in_bz` to verify the reduction and may recursively call itself with an increased search range if the initial attempt fails.

### Private Functions (Notable examples):

*   **`point_is_in_bz(point, rec_latt, origin_point, tolerance, verbose)` (Function)**
    *   **Purpose**: Checks if a `point` lies within the first Brillouin Zone.
    *   **Description**: A point `k` is in the first BZ if it is closer to the origin (or `origin_point`) than to any other reciprocal lattice point `G`. That is, `|k - origin| <= |k - G|` for all `G`. The function checks this condition for `G` vectors formed by `+/- rec_latt(i,:)`.

## Important Variables/Constants

This module does not define its own public parameters. It relies on:
*   The `crystal_3D` derived type (and others like `dp`, `default_tol_for_vec_equality`) imported from `constants_and_types_mod`.
*   Mathematical functions and constants from `math_mod` and `constants_and_types_mod`.

## Usage Examples

Other modules, particularly `main_BandUP.f90` and `symmetry_mod.f90`, use the routines from `crystals_mod`.

```fortran
! Example of creating a crystal and getting its reciprocal lattice
USE constants_and_types, ONLY: dp, crystal_3D
USE crystals, ONLY: create_crystal, get_rec_latt

TYPE(crystal_3D) :: silicon_cell
REAL(KIND=dp), DIMENSION(3,3) :: Si_latt_vecs, Si_rec_latt_vecs
REAL(KIND=dp) :: Si_alat

Si_alat = 5.43_dp
Si_latt_vecs(:,1) = (/ Si_alat/2.0_dp, Si_alat/2.0_dp, 0.0_dp /)
Si_latt_vecs(:,2) = (/ Si_alat/2.0_dp, 0.0_dp, Si_alat/2.0_dp /)
Si_latt_vecs(:,3) = (/ 0.0_dp, Si_alat/2.0_dp, Si_alat/2.0_dp /)

CALL create_crystal(silicon_cell, latt_vecs=Si_latt_vecs, description="Silicon FCC")
! silicon_cell%rec_latt_vecs is now populated by create_crystal,
! or one could call get_rec_latt explicitly:
! CALL get_rec_latt(silicon_cell%latt_vecs, Si_rec_latt_vecs)
```

## Dependencies and Interactions

*   **`constants_and_types`**: Essential for the `crystal_3D` type, precision kinds (`dp`), and default tolerance parameters (e.g., `default_tol_for_vec_equality`, `default_tol_for_int_commens_test`, `twopi`, `pi`).
*   **`math`**: For mathematical operations like `coords_cart_vec_in_new_basis` (coordinate transformation), `norm` (vector norm), `same_vector` (vector comparison), `cross` (cross product), `triple_product`, `inverse_of_3x3_matrix`.
*   **`lists_and_seqs`**: For the `append` utility, used internally in `create_crystal` for handling atomic symbols.

This module provides fundamental crystal manipulation routines that are used by higher-level modules for setting up calculations and performing symmetry analyses.
