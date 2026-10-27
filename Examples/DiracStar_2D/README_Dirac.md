# Dirac stars in axisymmetry (DiracStar_2D)

The axisymmetric (modified cartoon) version of the Dirac stars of Spherical_BS_Evolver (`dirac-stars` branch): two
Dirac spinors of mass mu in the j = 1/2, m_j = +-1/2 states, coupled to CCZ4. Conventions are those of the 1D code,
signature (-,+,+,+).

## Files

| File | Contents |
|---|---|
| `Source/Matter/DiracFieldCartoon.hpp` | Spinor variables, 3+1 Dirac right-hand side and stress tensor on the cartoon plane |
| `Source/Cartoon/CCZ4CartoonDirac.hpp` | CCZ4 cartoon RHS with Dirac matter |
| `Source/Cartoon/ConstraintsCartoonDirac.hpp` | Ham/Mom with Dirac sources; also the diagnostics `n`, `n_psi`, `N` and `mod_psi` (see below) |
| `Source/Cartoon/CCZ4CartoonBase.hpp`, `ConstraintsCartoonBase.hpp` | Vacuum CCZ4 cartoon RHS/constraints, now shared by the scalar and Dirac versions (moved verbatim out of `CCZ4Cartoon` / `ConstraintsCartoon`; the scalar-field RHS and constraints are bit-identical to before) |
| `Source/InitialConditions/DiracStars/DiracStarSolution.hpp` | Port of the 1D Dirac star solver (shooting, isotropic transformation, Gaussian perturbation) |
| `Source/InitialConditions/DiracStars/DiracStar.hpp` | Sets the 2D initial data from the 1D solution: a single star at rest or boosted along the axis, or a head-on binary along the axis |
| `Source/InitialConditions/DiracStars/DiracStarBoost.hpp` | Boosted star: lab metric, lapse, shift and K_ij as for the boson stars, plus the spinor transformation (shared with ExoZvezda) |
| `Source/TaggingCriteria/DiracAndChiExtractionTaggingCriterion.hpp` | Tagging on the spinor and chi, as for the scalar field |
| `Examples/DiracStar_2D/` | Level, parameters (`params.txt`), variables |
| `derivations/` | Python derivation checks of the equations (`check_dirac_cartoon.py`) and the reference implementation (`dirac_ref.py`) |
| `tests/` | C++ vs Python test of the RHS (`test_dirac_cpp.*`, no Chombo needed) and a single-box test of the assembled code (`box_test/`) |
| `family/` | Ground-state family (`dirac_family.dat`) and the mass / binding-energy vs central-amplitude plot (`plot_family.py`) |

The scalar field variables are removed from this example by `#define CARTOON_NO_SCALAR_FIELD` in `UserVariables.hpp`
(guards in `CCZ4CartoonUserVariables.hpp` and `ADMConformalVars.hpp`; the other examples are unaffected).

## Equations

**Matter content.** Only the J_x = +1/2 spinor psi is evolved (x is the symmetry axis). The J_x = -1/2 spinor is its
mirror image in the plane z = 0, which is a symmetry of any twist-free axisymmetric spacetime. On the cartoon plane its
stress tensor therefore equals that of psi with the z-odd components reversed. The sources are 2 x the z-even part of
T[psi].

**Variables.** There are 4 complex components (8 real variables), which are frame components on the plane. The spin is
quantised along x.

| Component | Meaning | Parity in y | Parity in x |
|---|---|---|---|
| `psi0` | upper, spin up | even | even |
| `psi1` | upper, spin down | odd | odd |
| `psi2` | lower, spin up | even | odd |
| `psi3` | lower, spin down | odd | even |

The x parities apply only for equatorial symmetry, i.e. a reflective boundary at x = 0. For the spherical star:

    psi = e^{-i omega t} (g, 0, -i f x/r, -i f y/r) / sqrt(4 pi)

**Tetrad.** Time gauge e_0 = n and the symmetric triad e = gamma^{-1/2} (Cartesian components). The triad is an
algebraic function of gamma_ij, so it needs no extra variables, and it is smooth on the axis. On the plane the 2x2
square root has a closed form. The spin connection includes the rotation of the triad as gamma_ij evolves.

**Evolution equation.**

    d_t psi = beta^k d_k psi - alpha Gamma_n psi + alpha gamma^0 gamma^b e_b^k (d_k + Gamma_k) psi + alpha mu gamma^0 psi

**Cartoon terms.** The z-derivatives on the plane follow from axisymmetry:

    d_z psi = (0, i psi1, 0, i psi3) / y

and d_z of the tensors is [G, T] / y.

**Sources.** The time derivative of psi in T_ab is eliminated with the equation itself, so rho, S_i and S_ij need only
slice data.

**Verification.** `derivations/check_dirac_cartoon.py` checks:

- the 3+1 form against the covariant 4D Dirac operator for arbitrary metric functions (2e-13);
- the cartoon reduction against the full 3D equations on axisymmetric data, the equivariance that fixes the spinor
  rotation generator, and the mirror-spinor doubling of the sources (all <= 3e-13);
- the reduction to the 1D (Spherical_BS_Evolver) equations and sources for spherical data (1e-14).

## Tests run so far (no 2D evolutions)

- **C++ vs Python reference** at 200 random points: 3e-15.
- **Ported 1D solver vs Spherical_BS_Evolver:**
  - M, N, omega and r_99 agree at A = 0.05 to 1.0;
  - the isotropic profiles, including the perturbed ones, agree to the precision of the 1D output (5e-10).
- **Single-box test (`tests/box_test`)**, 2D initial data plus one evaluation of the full RHS and constraints at
  dx = 0.1 and 0.05:

  | Star | Quantity | dx = 0.1 | dx = 0.05 |
  |---|---|---|---|
  | A = 0.15 | \|d_t psi + i omega psi\| / \|omega psi\| | 4.5e-7 | 2.8e-8 |
  | A = 0.15 | metric RHS | 1.5e-7 | 1.0e-8 |
  | A = 0.15 | Ham / 16 pi rho | 1.5e-6 | 9.9e-8 |
  | A = 0.30 (unstable) | \|d_t psi + i omega psi\| / \|omega psi\| | 2.1e-5 | 1.3e-6 |
  | A = 0.30, perturbed | Ham / 16 pi rho | 1.4e-5 | 9.1e-7 |

  All of these converge at 4th order. For the perturbed star, d_t psi stays at a resolution-independent 1.6e-3, as it
  should.
- **Boosted star (rapidity 0.3, `box_test` mode 1, 1D solution with 80001 points):** d_t psi equals the exact
  derivative of the moving solution, and Ham and Mom converge:

  | Quantity | dx = 0.1 | dx = 0.05 | factor |
  |---|---|---|---|
  | d_t psi error | 6.3e-7 | 3.9e-8 | 16.4 |
  | Ham / 16 pi rho | 1.3e-6 | 8.1e-8 | 16.0 |
  | Mom | 6.1e-8 | 3.8e-9 | 15.9 |

  The default 20001 points leave a floor of about 1e-8 from the interpolated psi' and alpha'.
- **Head-on binary (mode 2):** for all three `id_choice` options the data are exactly mirror-symmetric,
  psi(-x) = diag(1, -1, -1, 1) psi(x), with the metric even and the shift odd.
- **Aspherical perturbation (mode 3, m = 3, delta = 0.05):**
  - the spinor equals exp(-2 delta cos(m phi)) times the unperturbed one to round-off, and the metric is unchanged;
  - Mom stays zero (1e-16);
  - Ham / 16 pi rho = 0.18 at both resolutions, the expected O(delta) violation since the constraints are not
    re-solved.
- **Compilation:** both examples (DiracStar_2D and ScalarField_2D) compile against Chombo (2D, MPI, gcc).
- **Scalar-field regression:** the scalar RHS and constraints on a test box are bit-identical before and after the
  refactor.

## Boosted stars and head-on binaries

These follow the boson stars of `BosonStar.impl.hpp`.

- **Single star:** `DS_rapidity` boosts the star along the symmetry axis; rapidity > 0 moves it towards -x.
- **Binary:** `DS_binary = 1` gives a head-on binary along the axis, with `DS_separation`, `DS_rapidity2`,
  `central_amplitude_dirac2`, `phase2` and `mass_ratio`.
- **Metric:** `id_choice` 0, 1 or 2 selects plain superposition, the equal-mass fix or the unequal-mass fix.
- **Spinor:** psi_lab = U(eta) e^{-i omega t~} psi_0(x~), with

      U = cosh(eta/2) - sinh(eta/2) alpha_x,   tanh(eta) = tanh(xi) alpha / psi^2,   alpha_x block off-diagonal sigma_3

  U is the local Lorentz boost between the star's rest tetrad and the lab tetrad. A boost along the axis preserves
  axisymmetry and the mirror relation between the two spinors.
- **Derivation:** checked numerically in ExoZvezda's `Examples/SingleDiracStar/derivations/check_dirac_boost.py`
  (dirac-stars branch).
- **Interference:** both stars occupy the same spinor fields, so where they overlap the fields interfere (as for
  boson stars) and the relative phase matters.

## The family (mu = G = 1)

From `family/plot_family.py`:

- **Maximum mass and particle number:** at A_c = 0.2288, with M = 0.70891 and N = 0.73686.
  - Below this the ground state is perturbatively stable (the fundamental radial mode).
  - Beyond it one radial mode is unstable. Further turning points at A_c = 1.240 and 2.363 add one more unstable mode
    each.
- **Binding energy:** E_b = M - mu N is minimal (-0.0279) at the maximum mass and changes sign at A_c = 0.460.
- **1D runs:**
  - S1 (A = 0.05) and S2 (0.15) are stationary.
  - U1 (0.30) and U2 (0.40) collapse for a perturbation a = +0.01 and migrate to the stable branch for a = -0.01.
  - U3 (0.60, unbound) collapses.

## Diagnostics

All of these are plottable (`plot_vars`). The `n` variables are the particle number density measured by the normal
observers, n = -n_mu J^mu = psi^dag psi summed over the spinors.

| Name | Meaning |
|---|---|
| `n` | number density of both spinor fields |
| `n_psi` | number density of each spinor field: the evolved spinor and its mirror image have equal densities on the cartoon plane, so n = 2 n_psi |
| `N` | n chi^{-3/2} 2 pi y, the integrand summed into `ParticleNumber` |
| `mod_psi` | sqrt(2 pi n), which equals A_central at the centre of a star at rest |
| `rho` | energy density |

## Perturbations

`dirac_perturb = 1` perturbs a single star. In this generalised formalism the azimuthal index
`dirac_perturb_m` selects the perturbation, and the spherical perturbation of Spherical_BS_Evolver is its m = 0 member.

- **m = 0 (spherical):** the 1D solver applies (f, g) -> (1 + a exp(-(r - r0)^2 / w^2)) (f, g) and re-solves the
  constraints (metric), as before. Here a = `dirac_perturb_amp`, w = `dirac_perturb_spread` and
  r0 = `dirac_perturb_center`.
- **m >= 1 (aspherical):** the perturbation of arXiv:2609.30021, Eq. (11). The matter fields (the spinor components) are multiplied
  by exp(-2 delta cos(m phi)), with delta = `dirac_perturb_amp`.
  - phi is the azimuthal angle about the star centre with respect to the symmetry axis and the cylindrical radius,
    i.e. the polar angle from the axis.
  - In the cartoon code the perturbation is therefore axisymmetric. The two-axis perturbations of the paper are only
    in ExoZvezda.
  - As in the paper, the constraints are not re-solved, so Ham is violated at O(delta).
  - F = exp(-2 delta cos(m theta)), theta the polar angle, is smooth away from the centre. ExoZvezda's smoothed option
    (`dirac_perturb_smooth`), which removes the discontinuity of F along the line normal to the perturbation plane,
    is therefore not needed here; on the cartoon plane it would give the same F anyway.
- **Switching off:** `dirac_perturb_aspherical` (default 0) must be 1 for m >= 1. With 0 the aspherical perturbations
  are disabled and the perturbation is spherical whatever m is.
- **Restrictions:**
  - aspherical perturbations are for single stars only (not with the `binary` option);
  - odd m is rejected with a reflective x boundary (equatorial symmetry), since cos(m phi) is then odd.

The code is in `Source/utils/AsphericalPerturbation.hpp`, shared with ExoZvezda.

## Running

    cd Examples/DiracStar_2D && make all && mpirun ./Main_DiracStar2D*.ex params.txt

**Initial-data parameters** (`central_amplitude_dirac`, `dirac_*`) mirror Spherical_BS_Evolver's `BSParams.par`:

- `dirac_perturb*` give the spherical and aspherical perturbations described below.
- `dirac_solver_radius` should be at least about 3 r_99.
- The frequency search brackets automatically.

**Outputs:**

- `ParticleNumber` (conserved apart from outflow);
- `mod_psi_max`, which equals the central amplitude A_c for a star at rest;
- `constraint_norms` and `min_chi`, as for the boson stars. The ADM integrands `M_adm` and `P_adm` are plottable, for
  integration over an extraction sphere.

**Apparent horizons.** The example compiles with or without `USE_AHFINDER` (PETSc), like the scalar-field examples.
With it, `AH_activate = 1` puts a horizon finder on the star (on each star plus their merger for a binary), e.g. to
detect the collapse of an unstable star. The initial guesses are `AH_1_initial_guess` / `AH_2_initial_guess`
(default 0.3), contracted along the axis by sqrt(1 - v^2) for boosted stars.

The grid, AMR and tagging values in `params.txt` are untested starting points.
