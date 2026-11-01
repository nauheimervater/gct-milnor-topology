# gct-milnor-topology: Topological Barriers in Geometric Complexity Theory & Exact Cohomology of $F_{\text{Perm}_3}$

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Paper: PDF](https://img.shields.io/badge/Paper-PDF-red.svg)](paper/Topological_Barriers_in_GCT_Permanent_vs_Determinant.pdf)

> **Topological Barriers in Geometric Complexity Theory: Milnor-Fiber Padding Collapse and the Exact Cohomology of the 3×3 Permanent**  
> Thomas Nauheimer (`nauheimer.t@gmail.com`), November 2026.  
> Paper: [`paper/Topological_Barriers_in_GCT_Permanent_vs_Determinant.pdf`](paper/Topological_Barriers_in_GCT_Permanent_vs_Determinant.pdf) (Markdown source: [`paper/Topological_Barriers_in_GCT_Permanent_vs_Determinant.md`](paper/Topological_Barriers_in_GCT_Permanent_vs_Determinant.md), LaTeX source: [`paper/Topological_Barriers_in_GCT_Permanent_vs_Determinant.tex`](paper/Topological_Barriers_in_GCT_Permanent_vs_Determinant.tex)).

---

## Summary of Proved Results

In Valiant's algebraic complexity theory ($\mathsf{VP}$ vs. $\mathsf{VNP}$), the multiplicative padding formulation asks whether $\text{pad}(\text{Perm}_n, m) = z^{m-n}\text{Perm}_n(X)$ lies in the $GL_{m^2}(\mathbb{C})$-orbit closure of the determinant $\text{Det}_m(Y)$. After the refutation of naive representation-theoretic occurrence obstructions for padded polynomials (Ikenmeyer & Panova, 2017; Bürgisser, Ikenmeyer & Panova, 2019), a natural question is whether topological invariants of the Milnor fibers $F_{\text{Det}_m} = SL_m(\mathbb{C})$ and $F_{\text{pad}(\text{Perm}_n, m)}$ can obstruct orbit-closure containment.

This repository provides the paper and self-contained symbolic (`SymPy`) and exact finite-field (`NumPy`) verification scripts for five theorems:

1. **Theorem 1 (Padding Fibration & $m$-Independence of Betti Numbers):** For all $m > n$, $F_{\text{pad}(\text{Perm}_n, m)} \to \mathbb{C}^*$ is a locally trivial fiber bundle with fiber $F_{\text{Perm}_n}$ whose Leray–Serre spectral sequence degenerates at $E_2$, yielding
   $$b_k(F_{\text{pad}(\text{Perm}_n, m)}) = b_k(F_{\text{Perm}_n}) + b_{k-1}(F_{\text{Perm}_n}),$$
   completely independent of the padding dimension $m$. Consequently, no naive Hamm–Lefschetz cohomology injection $H^k(SL_m) \hookrightarrow H^k(F_{\text{pad}})$ can hold in high degrees ($k > n^2$, e.g. $b_{m^2-1}(SL_m) = 1 > 0 = b_{m^2-1}(F_{\text{pad}})$).
2. **Theorem 2 (Non-Freeness of the Scaling Torus & Borel Equivariant Torsion):** The row/column scaling torus $T \cong (\mathbb{C}^*)^{2n-1}$ does **not** act freely on $F_{\text{Perm}_n}$: $\dim \text{Stab}_T(I_n) = n$ at permutation matrices and $\dim \text{Stab}_T(X) = 1$ on dense matrices $X \in (\mathbb{C}^*)^{n\times n}$. However, $F_{\text{Perm}_n}^T = \emptyset$ for all $n \ge 2$, so by Atiyah–Bott–Segal localization $H^*_T(F_{\text{Perm}_n}, \mathbb{C})$ is a pure torsion $H^*(BT)$-module, forcing $b_{n^2-1}(F_{\text{Perm}_n}) \le \sum_{k=0}^{n^2-2} b_k(F_{\text{Perm}_n})$.
3. **Theorem 3 (Primary Decomposition of $\text{Sing}(\text{Perm}_3 = 0)$):** Over $\mathbb{C}$, the Jacobian ideal $J(\text{Perm}_3)$ decomposes into **15 prime ideals of Krull dimension $s = 3$** ($3$ two-zero-row planes $R_{ij}$, $3$ two-zero-column planes $C_{ij}$, and $9$ cross-minor quadric cones $D_{ij}$). By Kato's connectivity theorem, $F_{\text{Perm}_3} \subset \mathbb{C}^9$ is **4-connected** ($b_1 = b_2 = b_3 = b_4 = 0$).
4. **Theorem 4 (Exact Point Count, Mixed Tate Structure, and Betti Numbers of $F_{\text{Perm}_3}$):** Projecting onto the first two rows exhibits $F_{\text{Perm}_3}$ as an affine $\mathbb{A}^2$-bundle over $U = \mathbb{C}^6 \setminus V$, where $V$ admits a disjoint stratification into $\mathbb{A}^3 \sqcup (\mathbb{A}^3 \setminus \{0\}) \sqcup \bigsqcup_{l=1}^3 (\mathbb{G}_m)^2 \sqcup \bigsqcup_{k=1}^3 (\mathbb{G}_m)^3$. This gives:
   - Exact finite-field point count for all $\text{char}(\mathbb{F}_q) \neq 2$:
     $$|F_{\text{Perm}_3}(\mathbb{F}_q)| = q^8 - 5q^5 + 6q^4 - 3q^3 + q^2 = |SL_3(\mathbb{F}_q)| + q^2(q-1)^4.$$
   - Full-rank stratified boundary maps ($\text{rank}(\delta_4) = 3, \text{rank}(\delta_3) = 6, \text{rank}(\delta_2) = 2$), proving unconditionally that $H^m_c(V, \mathbb{Q})$ is pure of a single Tate weight in every degree and yielding the **exact Betti numbers**:
     $$b(F_{\text{Perm}_3}) = (1, 0, 0, 0, 0, 5, 6, 3, 1), \quad \sum_{k=0}^8 b_k(F_{\text{Perm}_3}) = 16 \;\; (\text{vs. } 4 \text{ for } SL_3(\mathbb{C})), \quad \chi(F_{\text{Perm}_3}) = 0.$$
5. **Theorem 5 (Padding Singularity Collapse & Bernstein–Sato Specialization Counterexample):** For $m - n \ge 2$, $z^{m-n}\text{Perm}_n(X)$ has a codimension-1 singular locus $\{z=0\}$ ($\mu^{(d)} = \infty$, $\text{lct}_0 = \frac{1}{m-n} < 1$), and bare Bernstein–Sato roots fail subset containment under specialization ($f_t = x^2 - ty^3 \to x^2$).

---

## Quickstart & Verification

```bash
git clone https://github.com/nauheimervater/gct-milnor-topology.git
cd gct-milnor-topology
pip install -r requirements.txt

# Run unit tests for all 5 theorems (< 0.2s)
python -m unittest discover -s tests -v

# Reproduce all verification tables in the paper
python examples/reproduce_all.py
```

### Verified Finite-Field Enumeration ($\mathbb{F}_p$, $p \in \{3, 5, 7, 11\}$)

| Prime $p$ | Exact $\|V(\mathbb{F}_p)\|$ | Exact $\|F_{\text{Perm}_3}(\mathbb{F}_p)\|$ | Exact $\|SL_3(\mathbb{F}_p)\|$ | Difference $p^2(p-1)^4$ | Formula Match |
|:---:|---:|---:|---:|---:|:---:|
| 3 | 89 | 5,760 | 5,616 | 144 | Exact |
| 5 | 489 | 378,400 | 372,000 | 6,400 | Exact |
| 7 | 1,441 | 5,694,192 | 5,630,688 | 63,504 | Exact |
| 11 | 5,961 | 213,637,600 | 212,427,600 | 1,210,000 | Exact |

### Exact Betti Numbers ($F_{\text{Perm}_3}$ vs. $SL_3(\mathbb{C})$ vs. $F_{\text{pad}(\text{Perm}_3, m)}$)

| Degree $j$ | $b_j(F_{\text{Perm}_3})$ | $b_j(SL_3(\mathbb{C}))$ | $b_j(F_{\text{pad}(\text{Perm}_3, m)})$ | Cohomological Origin in $V \subset \mathbb{C}^6$ |
|:---:|:---:|:---:|:---:|---|
| 0 | 1 | 1 | 1 | Connectedness of $U = \mathbb{C}^6 \setminus V$ |
| 1 | 0 | 0 | 1 | 4-connected ($\dim_{\mathbb{C}} V = 3$) |
| 2 | 0 | 0 | 0 | 4-connected ($\dim_{\mathbb{C}} V = 3$) |
| 3 | 0 | 1 | 0 | 4-connected ($\dim_{\mathbb{C}} V = 3$) |
| 4 | 0 | 0 | 0 | 4-connected ($\dim_{\mathbb{C}} V = 3$) |
| 5 | 5 | 1 | 5 | $H^6_c(V) = \mathbb{Q}(-3)^5$ ($5$ irreducible 3-folds $R_1, R_2, Q_1, Q_2, Q_3$) |
| 6 | 6 | 0 | 11 | $H^5_c(V) = \text{coker}(\delta_4: \mathbb{Q}^3 \hookrightarrow \mathbb{Q}^9) = \mathbb{Q}(-2)^6$ |
| 7 | 3 | 0 | 9 | $H^4_c(V) = \text{coker}(\delta_3: \mathbb{Q}^6 \hookrightarrow \mathbb{Q}^9) = \mathbb{Q}(-1)^3$ |
| 8 | 1 | 1 | 4 | $H^3_c(V) = \text{coker}(\delta_2: \mathbb{Q}^2 \hookrightarrow \mathbb{Q}^3) = \mathbb{Q}(0)^1$ |
| 9 | 0 | 0 | 1 | Top fiber-bundle shift $b_8(F_{\text{Perm}_3})$ in $F_{\text{pad}}$ |
| **Sum** | **16** | **4** | **32** | **Euler characteristic $\chi = 0$ in all three cases** |

---

## Citation

```bibtex
@misc{nauheimer2026gct,
  author       = {Thomas Nauheimer},
  title        = {Topological Barriers in Geometric Complexity Theory: Milnor-Fiber Padding Collapse and the Exact Cohomology of the 3x3 Permanent},
  year         = {2026},
  month        = {November},
  url          = {https://github.com/nauheimervater/gct-milnor-topology},
  note         = {Preprint}
}
```

## License

MIT. Copyright © 2026 Thomas Nauheimer.
