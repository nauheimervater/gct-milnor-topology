"""
Unit tests verifying all 5 theorems of the GCT Milnor-fiber paper.
Run with: python -m unittest discover -s tests
"""

import unittest
from gct_topology import (
    BETTI_PERM3,
    BETTI_SL3,
    count_points_fp,
    verify_jacobian_15_components,
    verify_padding_fibration_and_hamm_lefschetz,
    verify_padding_singularity_and_bs_counterexample,
    verify_stratified_cohomology_and_betti,
    verify_torus_stabilizers_and_fixed_points,
)


class TestGCTMilnorTopology(unittest.TestCase):
    def test_theorem1_padding_betti_and_hamm_lefschetz(self):
        res = verify_padding_fibration_and_hamm_lefschetz()
        self.assertEqual(res["betti_perm3"], (1, 0, 0, 0, 0, 5, 6, 3, 1))
        self.assertEqual(res["betti_padded_perm3"], (1, 1, 0, 0, 0, 5, 11, 9, 4, 1))
        self.assertEqual(res["sum_betti_padded"], 32)
        self.assertEqual(res["euler_padded"], 0)
        self.assertTrue(res["hamm_lefschetz_contradiction"])

    def test_theorem2_torus_stabilizers_and_fixed_points(self):
        for n in (2, 3, 4, 5):
            res = verify_torus_stabilizers_and_fixed_points(n)
            self.assertEqual(res["torus_dim"], 2 * n - 1)
            self.assertEqual(res["dim_stab_identity"], n)
            self.assertEqual(res["dim_stab_dense"], 1)
            self.assertFalse(res["action_is_free"])
            self.assertTrue(res["fixed_point_set_empty"])

    def test_theorem3_jacobian_15_prime_components(self):
        res = verify_jacobian_15_components()
        self.assertEqual(res["num_components"], 15)
        self.assertTrue(res["all_vanish"])
        self.assertTrue(res["all_dim_3"])
        self.assertEqual(res["max_singular_dim_s"], 3)
        self.assertEqual(res["kato_connectivity"], 4)

    def test_theorem4_finite_field_point_counts(self):
        for p in (3, 5, 7):
            res = count_points_fp(p)
            self.assertEqual(res["exact_V"], res["formula_V"])
            self.assertEqual(res["exact_Perm3"], res["formula_Perm3"])
            self.assertEqual(res["exact_SL3"], res["formula_SL3"])
            self.assertEqual(res["diff_Perm3_minus_SL3"], res["formula_diff"])

    def test_theorem4_stratified_cohomology_and_betti(self):
        res = verify_stratified_cohomology_and_betti()
        self.assertEqual(res["rank_delta_4"], 3)
        self.assertEqual(res["rank_delta_3"], 6)
        self.assertEqual(res["rank_delta_2"], 2)
        self.assertEqual(res["betti_perm3"], BETTI_PERM3)
        self.assertEqual(res["total_betti_perm3"], 16)
        self.assertEqual(res["total_betti_sl3"], 4)
        self.assertEqual(res["euler_characteristic"], 0)

    def test_theorem5_padding_singularity_and_bs_counterexample(self):
        res = verify_padding_singularity_and_bs_counterexample()
        self.assertTrue(res["grad_at_z0_all_zero"])
        self.assertEqual(res["codimension_of_z0"], 1)
        self.assertFalse(res["naive_root_containment_holds"])


if __name__ == "__main__":
    unittest.main()
