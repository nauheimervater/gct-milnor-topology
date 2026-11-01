"""
gct_topology: Symbolic and finite-field verification of topological and singularity
barriers for Permanent vs. Determinant Milnor fibers in Geometric Complexity Theory.

Author: Thomas Nauheimer (2026)
License: MIT
"""

from .milnor_permanent import (
    BETTI_PERM3,
    BETTI_SL3,
    count_points_fp,
    verify_jacobian_15_components,
    verify_padding_fibration_and_hamm_lefschetz,
    verify_padding_singularity_and_bs_counterexample,
    verify_stratified_cohomology_and_betti,
    verify_torus_stabilizers_and_fixed_points,
)

__all__ = [
    "BETTI_PERM3",
    "BETTI_SL3",
    "count_points_fp",
    "verify_jacobian_15_components",
    "verify_torus_stabilizers_and_fixed_points",
    "verify_stratified_cohomology_and_betti",
    "verify_padding_fibration_and_hamm_lefschetz",
    "verify_padding_singularity_and_bs_counterexample",
]
