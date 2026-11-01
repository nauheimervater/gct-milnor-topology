"""
Reproduce all symbolic and finite-field tables in the paper:
"On the Topological Collapse of the Milnor-Fiber Route for Permanent vs. Determinant
and the Exact Cohomology of F_{Perm_3}"

Run with: python examples/reproduce_all.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from gct_topology import (
    count_points_fp,
    verify_jacobian_15_components,
    verify_padding_fibration_and_hamm_lefschetz,
    verify_padding_singularity_and_bs_counterexample,
    verify_stratified_cohomology_and_betti,
    verify_torus_stabilizers_and_fixed_points,
)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print("=" * 78)
    print("1. THEOREM 1: PADDING FIBRATION & HAMM-LEFSCHETZ COLLAPSE")
    print("=" * 78)
    t1 = verify_padding_fibration_and_hamm_lefschetz()
    print(f"b_k(F_Perm3)            = {list(t1['betti_perm3'])} (sum = {sum(t1['betti_perm3'])})")
    print(f"b_k(F_pad(Perm3, m))    = {list(t1['betti_padded_perm3'])} (sum = {t1['sum_betti_padded']}, chi = {t1['euler_padded']})")
    print(f"Top degree of SL_4 (k=15): b_15(SL_4) = {t1['b_top_sl4']} vs. b_15(F_pad) = {t1['b_top_fpad']} -> Contradiction: {t1['hamm_lefschetz_contradiction']}")

    print("\n" + "=" * 78)
    print("2. THEOREM 2: SCALING TORUS STABILIZERS & BOREL EQUIVARIANT TORSION")
    print("=" * 78)
    print(f"{'n':>3} | {'dim(T)':>6} | {'dim Stab(I_n)':>13} | {'dim Stab(dense)':>15} | {'Free?':>6} | {'(F_Perm_n)^T == empty':>21}")
    print("-" * 78)
    for n in (2, 3, 4, 5):
        t2 = verify_torus_stabilizers_and_fixed_points(n)
        print(
            f"{n:>3} | {t2['torus_dim']:>6} | {t2['dim_stab_identity']:>13} | "
            f"{t2['dim_stab_dense']:>15} | {str(t2['action_is_free']):>6} | {str(t2['fixed_point_set_empty']):>21}"
        )

    print("\n" + "=" * 78)
    print("3. THEOREM 3: PRIMARY DECOMPOSITION OF SING(Perm_3 = 0) (15 PRIME COMPONENTS)")
    print("=" * 78)
    t3 = verify_jacobian_15_components()
    print(f"Verified prime components : {t3['num_components']} (all vanish: {t3['all_vanish']}, all Krull dim 3: {t3['all_dim_3']})")
    print(f"Singular locus dimension  : s = {t3['max_singular_dim_s']}")
    print(f"Kato connectivity bound   : N - 2 - s = 9 - 2 - 3 = {t3['kato_connectivity']} (=> b_1 = b_2 = b_3 = b_4 = 0)")
    print(f"Generic rank-1 minor (1,1): {t3['rank1_minor_11']} != 0 (rank-1 matrices are smooth for Perm_3)")

    print("\n" + "=" * 78)
    print("4. THEOREM 4A: EXACT FINITE-FIELD ENUMERATION OVER F_p (p = 3, 5, 7, 11)")
    print("=" * 78)
    print(f"{'p':>3} | {'|V(F_p)|':>10} | {'|F_Perm3(F_p)|':>14} | {'|SL_3(F_p)|':>14} | {'Diff = p^2(p-1)^4':>18} | {'Match':>5}")
    print("-" * 78)
    for p in (3, 5, 7, 11):
        c = count_points_fp(p)
        ok = (
            c["exact_V"] == c["formula_V"]
            and c["exact_Perm3"] == c["formula_Perm3"]
            and c["exact_SL3"] == c["formula_SL3"]
            and c["diff_Perm3_minus_SL3"] == c["formula_diff"]
        )
        print(
            f"{p:>3} | {c['exact_V']:>10} | {c['exact_Perm3']:>14} | "
            f"{c['exact_SL3']:>14} | {c['diff_Perm3_minus_SL3']:>18} | {str(ok):>5}"
        )

    print("\n" + "=" * 78)
    print("5. THEOREM 4B: STRATIFIED BOUNDARY RANKS & UNCONDITIONAL BETTI SPECTRUM")
    print("=" * 78)
    t4 = verify_stratified_cohomology_and_betti()
    print(f"Boundary ranks: rank(delta_4: Q^3->Q^9) = {t4['rank_delta_4']}, rank(delta_3: Q^6->Q^9) = {t4['rank_delta_3']}, rank(delta_2: Q^2->Q^3) = {t4['rank_delta_2']}")
    print(f"dim H^m_c(V) for m=6..2: {[t4['hc_V'][m] for m in (6, 5, 4, 3, 2)]}")
    print(f"Exact Betti numbers b_0..b_8(F_Perm3): {list(t4['betti_perm3'])} (Total = {t4['total_betti_perm3']}, chi = {t4['euler_characteristic']})")
    print(f"Exact Betti numbers b_0..b_8(SL_3)   : {list(t4['betti_sl3'])} (Total = {t4['total_betti_sl3']}, chi = 0)")

    print("\n" + "=" * 78)
    print("6. THEOREM 5: PADDING SINGULARITY & BERNSTEIN-SATO COUNTEREXAMPLE")
    print("=" * 78)
    t5 = verify_padding_singularity_and_bs_counterexample()
    print(f"Gradient of z^2 * Perm_2 vanishes identically on z=0: {t5['grad_at_z0_all_zero']} (codim = {t5['codimension_of_z0']})")
    print(f"BS roots of generic cusp f_t = x^2 - t*y^3 : {t5['roots_ft_cusp']}")
    print(f"BS roots of special fiber f_0 = x^2        : {t5['roots_f0_x2']}")
    print(f"Naive subset containment Roots(f_0) c Roots(f_t) holds? {t5['naive_root_containment_holds']}")


if __name__ == "__main__":
    main()
