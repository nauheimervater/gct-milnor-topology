# Topological Barriers in Geometric Complexity Theory: Milnor-Fiber Padding Collapse and the Exact Cohomology of the 3×3 Permanent

**Thomas Nauheimer — November 2026 — preprint, version 1.0**

Status: five proved theorems with symbolic and finite-field verification (`tests/`, `examples/reproduce_all.py`). No claim of resolving VP vs. VNP is made; the paper establishes structural barriers against naive Betti-number, torus-quotient, and bare Bernstein–Sato root obstructions for padded permanent orbit closures, and computes the exact unconditional cohomology of the 3×3 permanent Milnor fiber.

## Abstract

In Valiant's algebraic complexity theory, separating VP from VNP reduces to showing that the padded permanent pad(Perm_n, m) = z^{m-n} Perm_n(X) does not lie in the GL_{m^2}(ℂ)-orbit closure of the determinant Det_m(Y). Following the refutation of naive representation-theoretic occurrence obstructions for padded polynomials (Ikenmeyer & Panova, 2017; Bürgisser, Ikenmeyer & Panova, 2019), a natural question is whether topological invariants of the Milnor fibers F_{Det_m} = SL_m(ℂ) and F_{pad(Perm_n, m)}—combined with Hamm–Lefschetz subspace restriction and the (2n−1)-dimensional row/column scaling torus action—can obstruct orbit-closure containment.

We prove five results clarifying the exact scope and limits of this topological route:

1. **Padding fibration and m-independence of Betti numbers (Theorem 1).** For all m > n, projection onto the padding variable exhibits F_{pad(Perm_n, m)} → ℂ* as a locally trivial fiber bundle with fiber F_{Perm_n}. Its Leray–Serre spectral sequence degenerates at E_2, yielding b_k(F_{pad(Perm_n, m)}) = b_k(F_{Perm_n}) + b_{k-1}(F_{Perm_n}), completely independent of m. Consequently, no naive Hamm–Lefschetz cohomology injection H^k(SL_m) ↪ H^k(F_{pad}) can hold in degrees k ≥ n^2 + 1 (for instance, b_{m^2-1}(SL_m) = 1 > 0 = b_{m^2-1}(F_{pad}) for all m > n), because the degenerate padding subspace fails Whitney transversality to the projective closure of SL_m at infinity.
2. **Non-freeness of the scaling torus and Borel equivariant torsion (Theorem 2).** The scaling torus T ≅ (ℂ*)^{2n-1} does not act freely on F_{Perm_n}: the stabilizer has dimension n at the identity and permutation matrices, and dimension 1 on dense matrices X ∈ (ℂ*)^{n×n}. The naive quotient F_{Perm_n}/T is therefore not a principal bundle. However, we prove that for all n ≥ 2 the fixed-point locus F_{Perm_n}^T is empty. By Atiyah–Bott–Segal localization, the Borel equivariant cohomology H*_T(F_{Perm_n}, ℂ) is a pure torsion module over H*(BT, ℂ), forcing the top-degree Betti bound b_{n^2-1}(F_{Perm_n}) ≤ ∑_{k=0}^{n^2-2} b_k(F_{Perm_n}).
3. **Primary decomposition of Sing(Perm_3 = 0) (Theorem 3).** Over ℂ, the Jacobian ideal J(Perm_3) of the 3×3 permanent decomposes into **15 prime ideals of Krull dimension s = 3** (3 two-zero-row planes R_{ij}, 3 two-zero-column planes C_{ij}, and 9 cross-minor quadric cones D_{ij}). Unlike the determinant, generic rank-1 matrices are smooth points of {Perm_3 = 0}. By Kato's connectivity theorem, F_{Perm_3} ⊂ ℂ^9 is **4-connected** (b_1 = b_2 = b_3 = b_4 = 0).
4. **Exact point count, mixed Tate stratification, and Betti numbers of F_{Perm_3} (Theorem 4).** Projecting onto the first two rows exhibits F_{Perm_3} as an affine 𝔸^2-bundle over U = ℂ^6 \ V, where V = {C_{31} = C_{32} = C_{33} = 0} admits a disjoint stratification into linear and toric pieces: ℂ^3 ⊔ (ℂ^3 \ {0}) ⊔ _{l=1}^3 (ℂ*)^2 ⊔ _{k=1}^3 (ℂ*)^3. This proves unconditionally that F_{Perm_3} is mixed Tate, yields the exact finite-field point count |F_{Perm_3}(𝔽_q)| = q^8 − 5q^5 + 6q^4 − 3q^3 + q^2 = |SL_3(𝔽_q)| + q^2(q−1)^4 for all char(𝔽_q) ≠ 2, and—via full-rank stratified boundary maps δ_4, δ_3, δ_2—establishes the **exact Betti spectrum**:
   b(F_{Perm_3}) = (1, 0, 0, 0, 0, 5, 6, 3, 1), total Betti number ∑ b_k = 16 (versus 4 for SL_3(ℂ)), and Euler characteristic χ = 0.
5. **Codimension-1 padding singularity and Bernstein–Sato degeneration counterexample (Theorem 5).** For m − n ≥ 2, the factor z^{m-n} creates a codimension-1 singular locus {z = 0}, driving generic slice Milnor numbers to infinity and depressing the log canonical threshold at the origin to lct_0 = 1/(m−n) < 1 = lct_0(Det_m), which is compatible with lower semicontinuity. While disjoint-variable factorization preserves intrinsic Bernstein–Sato roots via Roots(b_{g·h}) = Roots(b_g) ∪ Roots(b_h), bare Bernstein–Sato roots do not satisfy subset containment under specialization, as shown by the explicit counterexample f_t = x^2 − t y^3 → x^2.

## 1. Background and problem setting

Let Det_m(Y) denote the determinant of an m×m matrix of indeterminates and Perm_n(X) the permanent of an n×n matrix. In multiplicative padding form, Valiant's VP ≠ VNP conjecture asks whether for m = poly(n),

```
pad(Perm_n, m)(X, z) = z^{m-n} Perm_n(X)  ∉  overline{GL_{m^2}(ℂ) · Det_m}.
```

The Geometric Complexity Theory (GCT) program (Mulmuley & Sohoni, 2001, 2007) originally proposed separating these orbit closures via representation-theoretic occurrence obstructions. After Ikenmeyer & Panova (2017) and Bürgisser, Ikenmeyer & Panova (2019) proved that occurrence obstructions cannot even surpass the classical determinantal complexity bounds due to the padding factor, attention turned toward multiplicity obstructions and geometric/topological invariants of the hypersurfaces and their Milnor fibers:

```
F_{Det_m} = { Y ∈ M_m(ℂ) : Det_m(Y) = 1 } = SL_m(ℂ) ≃ SU(m),
F_{Perm_n} = { X ∈ M_n(ℂ) : Perm_n(X) = 1 },
F_{pad(Perm_n, m)} = { (X, z) ∈ M_n(ℂ) × ℂ : z^{m-n} Perm_n(X) = 1 }.
```

The cohomology ring of F_{Det_m} = SL_m(ℂ) is the exterior algebra Λ_ℂ[y_3, y_5, ..., y_{2m-1}], so its Betti numbers b_k(SL_m) = P_{odd}(k, m) count partitions of k into distinct odd integers in {3, 5, ..., 2m−1}, with total Betti number 2^{m-1}. A natural topological separation strategy attempts to bound b_k(F_{Perm_n}) or b_k(F_{pad(Perm_n, m)}) using the (2n−1)-dimensional row/column scaling torus of the permanent and compare it against b_k(SL_m) via Hamm's affine Lefschetz theorem (Hamm, 1983). Sections 2–6 analyze each step of this strategy rigorously.

## 2. Padding fibration and failure of naive Hamm–Lefschetz restriction

### 2.1 Locally trivial fiber bundle over ℂ*

**Lemma 1.** *For all integers m > n ≥ 1, the projection π: F_{pad(Perm_n, m)} → ℂ*, (X, z) ↦ z, is a locally trivial algebraic fiber bundle with fiber F_{Perm_n}.*

*Proof.* On F_{pad(Perm_n, m)}, the equation z^{m-n} Perm_n(X) = 1 forces z ≠ 0. Because Perm_n(X) is homogeneous of degree n, locally on any analytic or étale chart of ℂ* where a branch of z^{(m-n)/n} is chosen, the rescaling (X, z) ↦ (z^{(m-n)/n} X, z) trivializes π with fiber {Perm_n(W) = 1} = F_{Perm_n}. Globally, passing to the n-fold cyclic cover w ↦ w^n = z of ℂ* shows that the monodromy of the bundle is of finite order n, hence acts trivially on rational/complex cohomology because it is isotopic to the identity inside the connected group C* acting by scalar multiplication. ∎

**Theorem 1 (Padding Betti formula and Hamm–Lefschetz collapse).** *For all m > n and all degrees k ≥ 0 (with b_{-1} = 0), the Betti numbers of the padded permanent Milnor fiber are independent of the padding dimension m and satisfy*

```
b_k(F_{pad(Perm_n, m)}) = b_k(F_{Perm_n}) + b_{k-1}(F_{Perm_n}).
```

*Consequently, for any m > n, there is no injective cohomology restriction H^k(SL_m(ℂ), ℂ) ↪ H^k(F_{pad(Perm_n, m)}, ℂ) in top degree k = m^2 − 1.*

*Proof.* In the Leray–Serre spectral sequence E_2^{p,q} = H^p(ℂ*, R^q π_* ℂ) ⇒ H^{p+q}(F_{pad}, ℂ), the local system R^q π_* ℂ is trivial because the scalar monodromy ∈ μ_n ⊂ ℂ* extends to the connected group ℂ*. Since ℂ* ≃ S^1 has cohomology only in p ∈ {0, 1}, all differentials d_r (r ≥ 2) vanish for degree reasons, so E_2 = E_∞ and H^k(F_{pad}) ≅ H^k(F_{Perm_n}) ⊕ H^{k-1}(F_{Perm_n}).

Since F_{Perm_n} ⊂ ℂ^{n^2} is a smooth affine hypersurface of complex dimension n^2 − 1, Andreotti–Frankel implies b_j(F_{Perm_n}) = 0 for all j > n^2 − 1, so b_k(F_{pad}) = 0 for all k > n^2. By contrast, SL_m(ℂ) ≃ SU(m) has top non-vanishing cohomology b_{m^2-1}(SL_m) = 1. For m > n we have m^2 − 1 > n^2, so an injection H^{m^2-1}(SL_m) ↪ H^{m^2-1}(F_{pad}) would require 1 ≤ 0, a contradiction. The affine Hamm–Lefschetz theorem does not apply to the specific linear subspace L cutting out a padded projection because L is parallel to asymptotic degenerations of SL_m at infinity and violates Whitney transversality on the projective closure. ∎

## 3. Non-freeness of the scaling torus and Borel equivariant torsion

Let T = {(U, V) ∈ (ℂ*)^n × (ℂ*)^n : ∏_{i=1}^n u_i v_i = 1} ≅ (ℂ*)^{2n-1} act on M_n(ℂ) by (U, V) · X = U X V. Since Perm_n(UXV) = (∏ u_i v_i) Perm_n(X), T preserves F_{Perm_n}.

**Theorem 2 (Stabilizer stratification and Borel equivariant torsion).** *For all n ≥ 2:*
1. *The action of T on F_{Perm_n} is **not free**: at the identity matrix I_n (and at every permutation matrix P_σ), the stabilizer Stab_T(I_n) = {(U, U^{-1}) : U ∈ (ℂ*)^n} ≅ (ℂ*)^n has dimension n; on the open dense torus (ℂ*)^{n×n} ∩ F_{Perm_n}, the stabilizer {(c I_n, c^{-1} I_n) : c ∈ ℂ*} ≅ ℂ* has dimension 1.*
2. *The fixed-point locus F_{Perm_n}^T is **empty**: F_{Perm_n}^T = ∅.*
3. *The Borel equivariant cohomology H*_T(F_{Perm_n}, ℂ) is a pure torsion module over R = H*(BT, ℂ) ≅ ℂ[t_1, ..., t_{2n-1}], and the top-degree Betti number satisfies*

```
b_{n^2-1}(F_{Perm_n}) ≤ ∑_{k=0}^{n^2-2} b_k(F_{Perm_n}).
```

*Proof.* (1) At X = I_n, U I_n V = I_n ⇔ v_i = u_i^{-1} for all i = 1..n, which automatically satisfies ∏ u_i v_i = 1 and gives a stabilizer isomorphic to (ℂ*)^n. If x_{ij} ≠ 0 for all i, j, then u_i v_j = 1 for all i, j implies u_1 = ... = u_n = c and v_1 = ... = v_n = c^{-1}, a 1-dimensional stabilizer. Thus the action is never free (in fact, the 1-dimensional diagonal subgroup {(c I_n, c^{-1} I_n)} acts trivially on the entire space M_n(ℂ)).

(2) Suppose X ∈ M_n(ℂ) is fixed by all of T. If some entry x_{ij} ≠ 0, then the character χ(U, V) = u_i v_j must be identically 1 on T = ker(Φ), where Φ(U, V) = ∏_{k=1}^n u_k v_k. Hence u_i v_j = Φ(U, V)^a for some a ∈ ℤ. Comparing exponents of u_i gives a = 1, whereas comparing exponents of u_k for k ≠ i (which exists since n ≥ 2) gives a = 0, a contradiction. Thus X = 0; since Perm_n(0) = 0 ≠ 1, F_{Perm_n}^T = ∅.

(3) By the Atiyah–Bott–Segal localization theorem, localizing with respect to S = R \ {0} gives S^{-1} H*_T(F_{Perm_n}, ℂ) ≅ S^{-1} H*_T(F_{Perm_n}^T, ℂ) = 0, so H*_T(F_{Perm_n}, ℂ) is a pure torsion R-module of rank 0. In the Borel Leray–Serre spectral sequence E_2^{p,q} = R^p ⊗ H^q(F_{Perm_n}, ℂ) ⇒ H^{p+q}_T(F_{Perm_n}, ℂ), every row E_2^{*,q} is a free R-module of rank b_q(F_{Perm_n}) that must be killed in rank by differentials d_r (r ≥ 2) before E_∞. Since no differential can enter the top row q = n^2 − 1 from above, the rank b_{n^2-1}(F_{Perm_n}) must be injected by outgoing differentials into rows q ≤ n^2 − 2, yielding the stated inequality. ∎

| Matrix size n | dim(T) = 2n−1 | dim Stab_T(I_n) | dim Stab_T(dense) | Action free? | F_{Perm_n}^T = ∅ |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 2 | 3 | 2 | 1 | False | True |
| 3 | 5 | 3 | 1 | False | True |
| 4 | 7 | 4 | 1 | False | True |
| 5 | 9 | 5 | 1 | False | True |

## 4. Primary decomposition of Sing(Perm_3 = 0)

**Theorem 3 (15 prime components of the 3×3 permanent singular locus).** *Over ℂ, the singular locus Σ = Sing({Perm_3 = 0}) = V(J(Perm_3)) ⊂ ℂ^9 is the union of **15 irreducible prime components**, each of Krull dimension s = 3:*
- *3 two-zero-row linear subspaces R_{12}, R_{13}, R_{23} ≅ ℂ^3;*
- *3 two-zero-column linear subspaces C_{12}, C_{13}, C_{23} ≅ ℂ^3;*
- *9 cross-minor quadric cones D_{ij} (1 ≤ i, j ≤ 3), where row i and column j vanish and the complementary 2×2 permanent minor vanishes: D_{ij} ≅ { (y_{11}, y_{12}, y_{21}, y_{22}) ∈ ℂ^4 : y_{11} y_{22} + y_{12} y_{21} = 0 }.*

*Consequently, dim_ℂ Sing({Perm_3 = 0}) = 3, and by Kato's connectivity theorem (Kato, 1975) the Milnor fiber F_{Perm_3} ⊂ ℂ^9 is at least (9 − 2 − 3) = **4-connected**, so b_1 = b_2 = b_3 = b_4 = 0.*

*Proof.* The 9 partial derivatives ∂Perm_3/∂x_{ij} are the complementary 2×2 permanent minors C_{ij}. Each R_{ab} and C_{ab} is defined by 6 linear variables, so ℂ[X]/P ≅ ℂ[t_1, t_2, t_3] is a domain of Krull dimension 3, and every 2×2 minor contains at least one row in {a, b} or column in {a, b}, hence vanishes. Each D_{ij} is defined by 5 linear variables (the 5 entries in row i ∪ col j) together with the non-degenerate rank-4 quadric C_{ij} = y_{11} y_{22} + y_{12} y_{21} = 0, which is irreducible of Krull dimension 4 − 1 = 3. Conversely, if X ∈ Σ has a row (say row 1) and a column (say col 1) vanishing, C_{11} = 0 places X in D_{11}; if no two rows and no row-column pair vanish, having all 2×2 permanent minors zero in characteristic ≠ 2 forces two rows or two columns to vanish (see Section 5.1). Note that each R_{ab} and C_{ab} lies in the boundary of the cross components D_{ij} where a free coordinate degenerates to 0, which explains why naive affine elimination solvers report 13 branches unless homogeneous primary decomposition is used. Finally, on a generic rank-1 matrix X = u v^T with u_i, v_j ≠ 0, ∂Perm_3/∂x_{11} = 2 u_2 u_3 v_2 v_3 ≠ 0, so unlike Det_3, generic rank-1 matrices are smooth points of {Perm_3 = 0}. ∎

## 5. Exact point count and unconditional cohomology of F_{Perm_3}

### 5.1 Affine 𝔸^2-bundle projection and toric stratification

Expanding Perm_3(X) along the third row gives

```
Perm_3(X) = x_{31} C_{31}(a, b) + x_{32} C_{32}(a, b) + x_{33} C_{33}(a, b) = 1,
```

where a = (x_{11}, x_{12}, x_{13}), b = (x_{21}, x_{22}, x_{23}) ∈ ℂ^3 are the first two rows, and

```
C_{31} = a_2 b_3 + a_3 b_2,    C_{32} = a_1 b_3 + a_3 b_1,    C_{33} = a_1 b_2 + a_2 b_1.
```

Let V = {(a, b) ∈ ℂ^6 : C_{31} = C_{32} = C_{33} = 0} and U = ℂ^6 \ V. For every (a, b) ∈ U, the vector (C_{31}, C_{32}, C_{33}) ∈ ℂ^3 \ {0} is non-zero, so the fiber of π: F_{Perm_3} → U over (a, b) is an affine 2-plane 𝔸^2. Hence π: F_{Perm_3} → U is an affine 𝔸^2-bundle, and **F_{Perm_3} is homotopy equivalent to U = ℂ^6 \ V**.

**Lemma 2 (Disjoint toric stratification of V).** *Over any field k of characteristic ≠ 2, V ⊂ k^6 is the union of 5 irreducible 3-dimensional components—two row 3-planes R_1 = {a = 0}, R_2 = {b = 0} and three quadric cones Q_k = {a_k = b_k = 0, C_{3k} = 0} (k = 1, 2, 3)—and decomposes into four pairwise disjoint locally closed strata:*
1. *R_1 ≅ 𝔸^3;*
2. *R_2 \ {0} ≅ 𝔸^3 \ {0};*
3. *Three 2-tori C_l^∘ = {(a, b) : a_k = b_k = a_m = b_m = 0, a_l ≠ 0, b_l ≠ 0} ≅ (𝔾_m)^2 for {k, l, m} = {1, 2, 3} (where two columns vanish and the remaining column has both entries non-zero);*
4. *Three 3-tori Q_k^∘ = Q_k \ (R_1 ∪ R_2 ∪ C_i ∪ C_j) = {(a, b) : a_k = b_k = 0, a_i, a_j, b_i ∈ 𝔾_m, b_j = −a_j b_i / a_i} ≅ (𝔾_m)^3 for k = 1, 2, 3.*

*Proof.* Suppose (a, b) ∈ V. If all three columns are non-zero, having any entry zero forces the entire corresponding row to be zero (e.g., a_1 = 0 ⇒ b_1 ≠ 0 ⇒ a_2 = a_3 = 0 via C_{33} = C_{32} = 0), while if all six entries are non-zero, C_{33} = C_{32} = 0 implies b_2/a_2 = −b_1/a_1 = b_3/a_3, and C_{31} = 0 then gives 2 a_2 a_3 (b_2/a_2) = 0, impossible since char(k) ≠ 2. Thus every point of V outside R_1 ∪ R_2 has at least one zero column (a_k, b_k) = (0, 0). If two columns vanish and the third has both entries non-zero, it lies in C_l^∘ ≅ (𝔾_m)^2; if exactly one column (a_k, b_k) vanishes and neither row vanishes, all four remaining entries must be non-zero (else a row or a second column would vanish) and b_j = −a_j b_i / a_i, giving Q_k^∘ ≅ (𝔾_m)^3. ∎

**Theorem 4 (Exact point count, mixed Tate structure, and Betti numbers of F_{Perm_3}).**
1. *For every finite field 𝔽_q with char(𝔽_q) ≠ 2,*

```
|V(𝔽_q)| = q^3 + (q^3 − 1) + 3(q − 1)^2 + 3(q − 1)^3 = 5q^3 − 6q^2 + 3q − 1,
|F_{Perm_3}(𝔽_q)| = q^2 (q^6 − |V(𝔽_q)|) = q^8 − 5q^5 + 6q^4 − 3q^3 + q^2 = |SL_3(𝔽_q)| + q^2(q − 1)^4.
```

2. *The compactly supported cohomology groups H^m_c(V, ℚ) are pure of a single Tate weight in every degree: H^6_c(V) = ℚ(−3)^5, H^5_c(V) = ℚ(−2)^6, H^4_c(V) = ℚ(−1)^3, H^3_c(V) = ℚ(0)^1, and H^m_c(V) = 0 for m ≤ 2.*
3. *The exact Betti numbers of F_{Perm_3} are*

```
b_0 = 1,  b_1 = b_2 = b_3 = b_4 = 0,  b_5 = 5,  b_6 = 6,  b_7 = 3,  b_8 = 1,
∑_{k=0}^8 b_k(F_{Perm_3}) = 16,      χ(F_{Perm_3}) = 0.
```

*Proof.* Part (1) follows immediately by summing the point counts of the four disjoint strata in Lemma 2. For (2) and (3), let Z = R_1 ∪ R_2 (two ℂ^3 planes meeting at {0}) and W = Z ∪ C_1 ∪ C_2 ∪ C_3 (adding the three column 2-planes C_l ≅ ℂ^2). Since Z, C = ⋃ C_l, and Z ∩ C (six coordinate lines ℂ^1 meeting at {0}) are unions of linear subspaces meeting at {0}, the Mayer–Vietoris sequence for W = Z ∪ C yields H^6_c(W) = ℚ(−3)^2, H^5_c(W) = 0, H^4_c(W) = ℚ(−2)^3, H^3_c(W) = ℚ(−1)^6, H^2_c(W) = ℚ(0)^2, and H^1_c(W) = 0.

Next, V \ W = _{k=1}^3 Q_k^∘ ≅ _{k=1}^3 (ℂ*)^3 has H^6_c = ℚ(−3)^3, H^5_c = ℚ(−2)^9, H^4_c = ℚ(−1)^9, H^3_c = ℚ(0)^3. In the long exact sequence of compactly supported cohomology for W ⊂ V, the connecting homomorphisms preserve Hodge weights and reduce to three explicit linear maps (`gct_topology/milnor_permanent.py`):
- δ_4: H^4_c(W) = ℚ(−2)^3 → H^5_c(V \ W) = ℚ(−2)^9 has full column rank 3;
- δ_3: H^3_c(W) = ℚ(−1)^6 → H^4_c(V \ W) = ℚ(−1)^9 has full column rank 6;
- δ_2: H^2_c(W) = ℚ(0)^2 → H^3_c(V \ W) = ℚ(0)^3 has full column rank 2.

Because δ_4, δ_3, δ_2 are all injective, ker(δ_4) = ker(δ_3) = ker(δ_2) = 0, so no mixed-weight extensions or cancellations occur in H^m_c(V): dim H^6_c(V) = 3 + 2 = 5, dim H^5_c(V) = 9 − 3 = 6, dim H^4_c(V) = 9 − 6 = 3, dim H^3_c(V) = 3 − 2 = 1, and H^2_c(V) = 0. Finally, since F_{Perm_3} ≃ U = ℂ^6 \ V is a smooth complex 6-fold homotopy equivalent to U, Poincaré duality on U and the excision sequence for V ⊂ ℂ^6 give b_j(F_{Perm_3}) = dim H^j(U) = dim H^{12-j}_c(U) = dim H^{11-j}_c(V) for all j ≥ 1. ∎

| Prime p | Exact \|V(𝔽_p)\| | Exact \|F_{Perm_3}(𝔽_p)\| | Exact \|SL_3(𝔽_p)\| | Difference p^2(p−1)^4 | Formula match |
|:---:|---:|---:|---:|---:|:---:|
| 3 | 89 | 5,760 | 5,616 | 144 | Exact |
| 5 | 489 | 378,400 | 372,000 | 6,400 | Exact |
| 7 | 1,441 | 5,694,192 | 5,630,688 | 63,504 | Exact |
| 11 | 5,961 | 213,637,600 | 212,427,600 | 1,210,000 | Exact |

| Degree j | b_j(F_{Perm_3}) | b_j(SL_3(ℂ)) | b_j(F_{pad(Perm_3, m)}) | Cohomological origin in V ⊂ ℂ^6 |
|:---:|:---:|:---:|:---:|---|
| 0 | 1 | 1 | 1 | Connectedness of U = ℂ^6 \ V |
| 1 | 0 | 0 | 1 | 4-connected (dim_ℂ V = 3) |
| 2 | 0 | 0 | 0 | 4-connected (dim_ℂ V = 3) |
| 3 | 0 | 1 | 0 | 4-connected (dim_ℂ V = 3) |
| 4 | 0 | 0 | 0 | 4-connected (dim_ℂ V = 3) |
| 5 | 5 | 1 | 5 | H^6_c(V) = ℚ(−3)^5 (5 irreducible 3-folds R_1, R_2, Q_1, Q_2, Q_3) |
| 6 | 6 | 0 | 11 | H^5_c(V) = coker(δ_4: ℚ^3 ↪ ℚ^9) = ℚ(−2)^6 |
| 7 | 3 | 0 | 9 | H^4_c(V) = coker(δ_3: ℚ^6 ↪ ℚ^9) = ℚ(−1)^3 |
| 8 | 1 | 1 | 4 | H^3_c(V) = coker(δ_2: ℚ^2 ↪ ℚ^3) = ℚ(0)^1 |
| 9 | 0 | 0 | 1 | Top fiber-bundle shift b_8(F_{Perm_3}) in F_{pad} |
| **Sum** | **16** | **4** | **32** | **Euler characteristic χ = 0 in all three cases** |

## 6. Singularity barriers under padding and Bernstein–Sato specialization

Beyond Betti numbers, one may ask whether classical singularity invariants at the origin—such as the log canonical threshold lct_0(f), restricted Milnor numbers μ^{(d)}(f), or the roots of the Bernstein–Sato polynomial b_f(s)—can separate pad(Perm_n, m) from the determinant orbit closure.

**Theorem 5 (Padding singularity collapse and Bernstein–Sato root non-monotonicity).**
1. *For m − n ≥ 2, the padded polynomial f_{pad} = z^{m-n} Perm_n(X) has singular locus containing the entire hyperplane {z = 0} of codimension 1. Consequently, the restriction of f_{pad} to any generic linear subspace L of dimension d ≥ 2 through the origin has a non-isolated singularity along L ∩ {z = 0}, so the Teissier restricted Milnor number is μ^{(d)}(f_{pad}) = ∞.*
2. *At the origin, lct_0(Det_m) = 1 (since b_{Det_m}(s) = ∏_{k=1}^m (s + k)), whereas by the disjoint-variable product rule,*

```
lct_0(z^{m-n} Perm_n(X)) = min(1/(m − n), lct_0(Perm_n)) = 1/(m − n) < 1   (for m − n ≥ 2).
```

*Because the log canonical threshold is lower semicontinuous under specialization (lct_0(f_0) ≤ lct_0(f_t)), the inequality lct_0(pad(Perm_n, m)) < lct_0(Det_m) is consistent with orbit-closure containment and provides no obstruction.*
3. *For any product f(x, y) = g(x) h(y) in disjoint variable sets, the D-module tensor factorization D_{x,y}(g·h)^s ≅ D_x g^s ⊗_{ℂ[s]} D_y h^s implies*

```
Roots(b_{z^{m-n} Perm_n}(s)) = { −1, −1/(m−n), ..., −(m−n−1)/(m−n) } ∪ Roots(b_{Perm_n}(s)),
```

*so the intrinsic Bernstein–Sato roots of Perm_n survive padding unshifted. However, bare Bernstein–Sato roots do **not** satisfy subset containment Roots(b_{f_0}) ⊂ Roots(b_{f_t}) under 1-parameter degenerations f_t → f_0: for f_t = x^2 − t y^3 → f_0 = x^2, Roots(b_{f_t}) = {−1, −5/6, −7/6} for t ≠ 0, whereas Roots(b_{f_0}) = {−1, −1/2}, and −1/2 ∉ Roots(b_{f_t}).*

*Proof.* (1) Every partial derivative of z^k Perm_n(X) (k = m−n ≥ 2) contains a factor of z^{k-1}, hence vanishes identically on {z = 0}. Any linear subspace L ⊂ ℂ^{n^2+1} of dimension d ≥ 2 intersects {z = 0} in dimension ≥ d − 1 ≥ 1, so 0 is a non-isolated singularity of f_{pad}|_L. (2) and (3) follow from standard D-module properties of b-functions (Kashiwara, 1976; Saito, 1994) and the explicit b-functions b_{x^2}(s) = (s+1)(s+1/2) and b_{x^2-y^3}(s) = (s+1)(s+5/6)(s+7/6). ∎

## 7. Conclusion and open problem

The results above close three naive shortcuts in the topological approach to GCT:
- Betti numbers of the padded Milnor fiber F_{pad(Perm_n, m)} decouple completely from the padding size m (Theorem 1), and naive Hamm–Lefschetz injections fail at infinity.
- The row/column scaling torus T^{2n-1} does not act freely on F_{Perm_n} (Theorem 2), though its fixed-point freeness makes H*_T(F_{Perm_n}) a pure torsion H*(BT)-module.
- For n = 3, F_{Perm_3} is 4-connected with exact Betti vector (1, 0, 0, 0, 0, 5, 6, 3, 1) and total Betti number 16 versus 4 for SL_3(ℂ) (Theorems 3 and 4): the permanent Milnor fiber is topologically **richer**, not simpler, than the determinant Milnor fiber, because the 5-component singular locus V ⊂ ℂ^6 introduces 12 extra toric cohomology classes across degrees 5, 6, 7.

Any viable singularity or topological obstruction for VP vs. VNP must therefore work either in the **de-padded slice geometry** Z(Perm_n) = Z(Det_m) ∩ L (requiring universal upper bounds over arbitrary non-generic subspaces L ⊂ M_m(ℂ)) or with semicontinuous invariants of the **Kashiwara–Malgrange V-filtration / Hodge spectrum** on vanishing cycles rather than bare Betti numbers or bare b-function roots.

## References

Bürgisser, P., Ikenmeyer, C., & Panova, G. (2019). No occurrence obstructions in geometric complexity theory. *Journal of the American Mathematical Society*, 32(1), 163–193.

Hamm, H. A. (1983). Lefschetz theorems for singular varieties. *Proceedings of Symposia in Pure Mathematics*, 40(1), 547–557.

Ikenmeyer, C., & Panova, G. (2017). Rectangular Kronecker coefficients and plethysms in geometric complexity theory. *Advances in Mathematics*, 319, 40–66.

Kashiwara, M. (1976). B-functions and holonomic systems: Rationality of roots of B-functions. *Inventiones Mathematicae*, 38(1), 33–53.

Kato, M. (1975). Topology of closed manifolds and singularities of hypersurfaces. *Topology*, 14(2), 107–113.

Mulmuley, K., & Sohoni, M. (2001). Geometric complexity theory I: An approach to the P vs. NP and related problems. *SIAM Journal on Computing*, 31(2), 496–561.

Mulmuley, K., & Sohoni, M. (2007). Geometric complexity theory II: Towards explicit obstructions for embeddings among class varieties. *SIAM Journal on Computing*, 37(2), 326–390.

Saito, M. (1994). On microlocal b-function. *Bulletin de la Société Mathématique de France*, 122(2), 163–184.
