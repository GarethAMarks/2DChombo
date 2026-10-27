2D version of GRChombo. Includes my modifications of the code for constructing and evolving boson stars in axisymmetry, based on ExoZvezda.

## Dirac stars (branch `dirac-stars`)

`Examples/DiracStar_2D` evolves Dirac stars in axisymmetry: two spinors coupled to CCZ4 with the modified cartoon
method, with initial data from a port of the Spherical_BS_Evolver Dirac star solver. The vacuum CCZ4 cartoon RHS and
constraints are now in `CCZ4CartoonBase` / `ConstraintsCartoonBase`, shared by the scalar-field and Dirac versions; the
scalar-field results are unchanged bit for bit. See `Examples/DiracStar_2D/README_Dirac.md`.

## Proca stars (branch `proca-stars`, built on `dirac-stars`)

`Examples/ProcaStar_2D` evolves Proca stars in axisymmetry: a complex Proca field (Zilhão–Witek–Cardoso variables
with Gauss-constraint damping) coupled to CCZ4 with the modified cartoon method. Initial data come from a port of the
Spherical_BS_Evolver Proca star solver, and stars can be boosted along the axis or set up as head-on binaries. The
branch reuses the cartoon base classes of `dirac-stars`. See `Examples/ProcaStar_2D/README_Proca.md`.
