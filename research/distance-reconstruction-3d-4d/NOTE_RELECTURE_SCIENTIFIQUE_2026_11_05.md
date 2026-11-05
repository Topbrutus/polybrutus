# A seven-point distance multiset with rigid reconstruction in R^3 and a noncongruent realization in R^4
## Note for independent mathematical review — 2026-11-05

**Author's project:** Topbrutus / Polybrutus v0.1.0.  
**Status:** computationally certified example; not peer-reviewed; novelty and priority NOT established.  
**Scope:** purely Euclidean distance geometry; no physics or boson claim.

### Abstract (English)
Let M be the multiset of squared distances
`{1^6, (7/4)^2, (9/4)^3, 3^6, (13/4)^2, 4^2}`,
with superscripts denoting multiplicities, not powers.
Two separately implemented exact enumerations (Gram-based and Cayley–Menger-based) find that every seven-point configuration in Euclidean 3-space realizing M is congruent to the Polybrutus v0.1.0 configuration, up to relabeling. In Euclidean 4-space an explicit, noncongruent realization of the same M exists, with affine rank four. We provide exact coordinates for this witness and runnable verification code. This is an instance-specific computational certificate, not a claim of a new general theorem.

### Résumé (français)
Pour un multiensemble précis de 21 distances au carré associées à sept points, deux énumérateurs exacts distincts établissent une seule classe de congruence dans R³, sans hypothèse de cosphéricité. Une construction explicite en R⁴ réalise les mêmes distances non étiquetées avec un rang affine de quatre, donc n'est pas congruente à l'original. La nouveauté dans la littérature n'a pas été déterminée.

## 1. Definitions and target
Consider n=7 distinct points. A *distance multiset* is the collection of all 21 pairwise Euclidean distances, independent of which pair carries each value. It is equivalent to use squared distances because all distances are nonnegative.

The target multiset of squared distances, including multiplicities, is:
```text
value       1   7/4   9/4    3   13/4   4
count       6     2     3    6      2   2
sum of counts = 21
```

For exact computations use the integer colors `k_ij=4 d_ij^2` with multiplicities `{4:6,7:2,9:3,12:6,13:2,16:2}`.

**Definition of congruence:** Two configurations are congruent iff some Euclidean isometry followed by a permutation of the seven points maps one to the other. Affine rank is a congruence invariant.

## 2. Three-dimensional statement (computer-assisted)
**Claim A (certified by two enumerations).** Among configurations of seven distinct points in R³ whose 21 pairwise squared distances have the target multiset, there is exactly one congruence class. The representative is the Polybrutus v0.1.0 configuration. No spherical constraint is imposed on competitors.

**Gram enumeration.** For each affine dimension r in {3,2,1}, choose r+1 affinely independent anchors, and use permutation symmetry to put them first. Enumerate all assignments of available edge colors to their complete graph, modulo reordering of anchors. Build an exact, positive-definite Gram matrix from the assigned squared distances. Enumerate all assignments of remaining distances to the remaining vertices. A vertex's distances to r+1 independent anchors uniquely determine its position in that affine hull; compatibility and pairwise distances between new vertices are tested by integer Gram identities. All 21 distances must exhaust the target multiset. Compare each successful labeled distance matrix with the original under all permutations.

**Independent Cayley–Menger enumeration.** Start from the same distance multiset, but separately enumerate anchor edge assignments and reject inadmissible simplices using signs of exact Cayley–Menger determinants. Candidate vertex-anchor edges satisfy a Cayley–Menger vanishing relation on r+2 vertices; each pair of additional vertices is tested by related determinants. Every actual R^r realization satisfies these necessary tests, so they do not exclude genuine counterexamples. Completed matrices are checked against the full target multiset and the original's relabelings.

**Observed exact certificate:**
```text
                                   Gram enumerator       Cayley-Menger enumerator
Completed affine-rank-3 cases               27                       27
Of these congruent to original              27                       27
Completed affine-rank-2 cases                0                        0
Completed affine-rank-1 cases                0                        0
Different 3D configurations                 0                        0
```
The 27 enumerated labeled completions are not 27 distinct geometric classes. Both methods use exact integer arithmetic, canonical root labelings, and the same original multiset, but their geometric admissibility tests are separately implemented. The auxiliary SymPy checker validates Gram rank and PSD for all 27 completed outputs of the first enumerator.

**Scope of proof:** The certificate is a finite search whose completeness is justified by considering every affine rank 1, 2 and 3 and every assignment of the distance colors up to the symmetries of anchor relabeling. A peer reviewer should inspect the canonicalization and extension logic before treating it as a formal theorem.

## 3. Explicit four-dimensional witness
Let `s3=sqrt(3), s6=sqrt(6), s30=sqrt(30)`. The points below lie in an orthonormal coordinate system in R⁴:
```text
q1 = (   0,       0,        0,       0)
q2 = (   1,       0,        0,       0)
q3 = ( 1/2,   s3/2,        0,       0)
q4 = ( 1/2,  -s3/6,   2*s6/3,       0)
q5 = ( 3/8,11*s3/24, 13*s6/24,  s30/8)
q6 = (-1/8, -s3/24,  7*s6/24, -s30/8)
q7 = (   1,    s3/3,   2*s6/3,       0)
```
The exact matrix `K_ij=4||q_i-q_j||²` is:
```text
       0  4  4 12 12  4 16
       4  0  4 12 13  9 12
       4  4  0 16  9  9 12
      12 12 16  0  7  7  4
      12 13  9  7  0 13  4
       4  9  9  7 13  0 12
      16 12 12  4  4 12  0
```
Its 21 upper-triangular entries have exactly the target multiplicities. The 6x6 Gram matrix of `q2-q1,...,q7-q1` has **rank 4**, and its leading 4x4 determinant is exactly **15/16**. Thus `q1,...,q5` form a nondegenerate 4-simplex of volume `sqrt(15)/96`.

The source Polybrutus points have affine rank **3**. Since affine rank is invariant under all Euclidean isometries and point permutations, the witness is not congruent to Polybrutus. This part is a direct constructive proof, independent of any enumeration.

## 4. Falsifiability, known limitations
- A single R³ noncongruent configuration with the same multiset falsifies Claim A.
- Any nonzero error in one of the 21 exact witness distances falsifies the 4D construction.
- No uniqueness statement is claimed in dimensions 4, 5, 6 or higher.
- This is not a stability result: perturbed/noisy measurements require a different analysis.
- Neither mathematical priority nor any physical interpretation follows from these computations.
- The example is not a claim that unlabeled-distance reconstruction was first invented here.

## 5. Relation to previous work, sources consulted
1. M. Boutin and G. Kemper, 'On reconstructing n-point configurations from the distribution of distances or areas', *Advances in Applied Mathematics* **32**(4) (2004), 709–735. DOI: 10.1016/S0196-8858(03)00101-5.
2. M. Boutin and G. Kemper, 'Which point configurations are determined by the distribution of their pairwise distances?', *International Journal of Computational Geometry & Applications* **17**(1) (2007), 31–43. DOI: 10.1142/S0218195907002239.
3. 'Graph Reconstruction from Unlabeled Edge Lengths', *Discrete & Computational Geometry* (2021), DOI: 10.1007/s00454-021-00275-7.

The generic reconstruction results of Boutin–Kemper do not by themselves settle this degenerate example with repeated distances. Their framework supplies the correct terminology and highlights known non-reconstructible cases. The specific 3D/4D pair above was **not matched** in this targeted search; this does not establish originality.

## 6. Reproduction & requested peer review
Files in the reviewer package:
- `laboratoire_spectral.py`: exact input coordinates and distances.
- `polybrutus_unicite_sans_sphere_2026_11_05.py`: Gram-based exhaustive enumerator.
- `polybrutus_independant_cayley_menger_2026_11_05.py`: separate Cayley–Menger enumerator.
- `polybrutus_controle_sympy_global_2026_11_05.py`: exact verification of returned 3D cases.
- `polybrutus_contre_exemple_4d_exact_2026_11_05.py`: direct exact calculation from the seven R⁴ coordinates.
- Original execution logs, SHA-256 manifest and this note.

Suggested reviewer questions:
A. Does canonical anchor deduplication preserve every unlabeled realization and affine rank?
B. Does the Cayley–Menger extension test eliminate only impossible candidates?
C. Are the 21 upper-triangle entries and 4D Gram rank correct, independently computed?
D. Is an equivalent seven-point example with this multiset already in the literature?
E. Can a short non-computational proof replace or shorten the enumerative 3D part?

**Research status:** Ready for external *review request*, not ready for a claim of a novel physics result. No submission or contact has been made.

## Provenance
The research note and five verification programs originate from the local Polybrutus review package dated 2026-11-05. Public repository integration preserves Polybrutus v0.1.0 and adds an independent research supplement. No claim of peer review, discovery priority, Zenodo publication or new DOI is made here.
