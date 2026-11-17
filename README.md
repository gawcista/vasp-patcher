# Wannier models and electrostatic interactions in VASP

[English](#english) · [中文说明](#中文说明) · [References / 参考文献](#references--参考文献)

| Feature / 功能 | English | 中文 |
|---|---|---|
| Compact-band SCDM / 紧凑选带 | [Motivation and use](#en-scdm) | [动机与使用](#zh-scdm) |
| `HFRCUT=-2`, BvK WS | [Native treatment and WS](#en-ws) | [原生处理与 WS](#zh-ws) |
| `SLAB_WS` | [Definition and derivation](#en-slab) | [定义与推导](#zh-slab) |
| Bare, Yukawa, RK, gates | [Physical kernels](#en-kernels) | [各物理核](#zh-kernels) |
| PAW and frozen-Wannier projection | [Workflow](#en-use) | [投影与工作流](#zh-use) |
| Physical and numerical convergence | [Interpretation](#en-convergence) | [收敛与解释](#zh-convergence) |
| Apply / reverse / update | [Installation](#en-install) | [安装升级](#zh-install) |

## English

### 1. The physical problem: from Bloch states to an interacting Wannier model

A low-energy model needs both an orbital basis and an interaction appropriate to
its environment. In a moiré material, a small set of bands may describe the
physics near the Fermi level even though a plane-wave calculation contains
thousands of bands. Localized Wannier orbitals turn these bands into hopping
amplitudes and matrix elements between spatial charge distributions.

For a scalar, spin-independent interaction, a useful building block is the
bilinear form below. Its density-direct element is the electrostatic interaction
between two normalized Wannier charge clouds:

$$
\mathcal B_K[\rho,\sigma]
=\int d^3r\,d^3r'\,\rho^*(\mathbf r)
K(\mathbf r,\mathbf r')\sigma(\mathbf r'),\qquad
D_{mn}(\mathbf R)=\mathcal B_K[|w_{m0}|^2,|w_{n\mathbf R}|^2].
$$

Onsite repulsion depends on the extent and internal structure of the orbitals;
off-center terms describe their spatial interaction range. Comparing these
energy scales with the band width helps assess an interacting model. General
transition-density matrix elements also contain exchange and orbital-changing
processes, and need not be positive. The native file convention is specified
in the [VASP VRijkl documentation](https://vasp.at/wiki/VRijkl).

This patch addresses three parts of that construction:

1. **Band-compact SCDM:** make the initial Wannier projection affordable in a
   large-band calculation.
2. **`HFRCUT=-2`:** represent the bare interaction on the lattice-dependent
   Wigner–Seitz domain of the finite k mesh.
3. **Layered interactions:** retain the normal extent of the Wannier densities
   while introducing an ideal polarizable sheet or metallic gates.

The physical kernel, the finite integration domain and the Wannier basis are
separate choices. Changing the kernel at a fixed WAVECAR and final WANPROJ
compares interactions in the same basis. A PBE calculation can supply these
orbitals; choosing a screened interaction afterward does not turn that
calculation into a hybrid-functional calculation.

<a id="en-scdm"></a>

### 2. Band-compact SCDM: selecting the relevant Hilbert space without storing all bands

#### Motivation and principle

SCDM constructs localized orbitals by selecting informative real-space columns
of a density matrix or a spectrally filtered projector. QR with column pivoting
identifies a compact set of anchor points; the selected columns provide an
initial localized basis. The original SCDM, periodic SCDM-k and entangled-band
formulations are described in Refs. [1–3](#references--参考文献).

A spectral filter gives greater weight to the states one wants to represent:

$$
P_f(\mathbf k)=\sum_{n\in\mathcal A}
f(E_{n\mathbf k})\,|u_{n\mathbf k}\rangle\langle u_{n\mathbf k}|,
\qquad
M_{\mathrm{SCDM}}\propto
N_{\mathrm{grid}}N_{\mathrm{spinor}}N_{\mathrm{active}}.
$$

Here the weights specify a Wannier subspace, rather than necessarily being
thermodynamic occupations. For example, an ERFC filter smoothly suppresses
states above an energy scale:

$$
f(E)=\frac12\mathrm{erfc}\!\left(\frac{E-\mu}{\sigma}\right).
$$

In this example, `CUTOFF_MU` sets the filter center and `CUTOFF_SIGMA` its energy
width. These are subspace-selection parameters; the center need not coincide
with the chemical potential.

The memory problem arises when SCDM work arrays retain the full VASP band
count even though Wannier90 has excluded most bands. Large moiré cells make
both the real-space grid and the total band count large. This patch applies
the existing `exclude_bands` mask to SCDM's working wavefunction, eigenvalue,
filter and transformation arrays. The retained original band indices remain
in ascending order, so each AMN row keeps its intended band meaning.

The dominant selected-band buffers scale with `N_active` instead of `NBANDS`.
This is a storage and indexing improvement to VASP's SCDM implementation.
Global wavefunction/WAVECAR dimensions and the native MMN temporary allocation
retain their existing sizes; other parts of VASP therefore still contribute
to total memory.

#### How to use it

Use the normal VASP/Wannier90 interface and choose the bands in the Wannier90
input. For an isolated set with exactly as many retained bands as Wannier
functions, unit weights provide a natural projector:

```text
LSCDM = .TRUE.
LWANNIER90 = .TRUE.
NUM_WANN = 8
CUTOFF_TYPE = none
```

The example assumes that the accompanying `exclude_bands` selection retains
eight suitable bands. For entangled bands, retain a larger window and supply
an explicit filter.

| Input | Meaning in the compact subspace |
|---|---|
| Wannier90 `exclude_bands` | Original VASP bands omitted from the SCDM working space |
| `NUM_WANN` | Number of Wannier functions to construct |
| `CUTOFF_TYPE=none` | Unit weights on the retained states |
| `CUTOFF_MU`, `CUTOFF_SIGMA` | Explicit spectral-filter parameters; ERFC uses the center/width convention above |

With exclusions and more active bands than Wannier functions, ERFC/Gaussian
require both parameters, Fermi requires `CUTOFF_MU`, and `num_wann` requires
`CUTOFF_SIGMA`. Other filter forms retain the conventions of this VASP version;
the ERFC width formula above should not be transferred to them automatically.
With exactly `NUM_WANN` active bands, an omitted filter also gives unit weights.

SCDM produces the initial projection. Symmetry-adapted Wannier localization
and export of the final WANPROJ are external steps; this patch consumes that
final transformation for the interaction calculation.

<a id="en-ws"></a>

### 3. `HFRCUT=-2`: why a lattice-shaped Coulomb cutoff is useful

#### 3.1 What the original VASP treatment represents

The reciprocal bare kernel is `4πC/Q²`, with an integrable singularity at
`Q=0`. A finite k mesh requires a prescription for this contribution. Native
VASP uses a probe-charge Ewald correction for a regular unscreened mesh;
related supported cases use its Gaussian auxiliary-function treatment.
Alternatively, `HFRCUT` selects a spherical cutoff. See the official
[HFRCUT](https://vasp.at/wiki/HFRCUT) and
[Coulomb singularity](https://vasp.at/wiki/Coulomb_singularity) descriptions.

For a full, unreduced uniform mesh, the relevant real-space repetition is the
Born–von Karman (BvK) superlattice:

$$
\mathbf A_i=N_i\mathbf a_i,\qquad
N_k=N_1N_2N_3,\qquad V_B=N_k\Omega.
$$

BvK describes the periodicity associated with sampling Bloch states. Ewald is
an electrostatic summation/convergence prescription on that periodic system.
The boundary one wants for a projected interaction must be specified in
addition to choosing a summation prescription.

For the regular unscreened case, the periodic Green function and the local
Madelung correction can be expressed as

$$
-\nabla^2G_B(\mathbf r)
=4\pi\left[\sum_{\mathbf S\in\Lambda_B}
\delta(\mathbf r-\mathbf S)-\frac1{V_B}\right],
\qquad \langle G_B\rangle=0.
$$

$$
K_{\mathrm{per}}(\mathbf r)=C\,[G_B(\mathbf r)-v_M],
\qquad
v_M=\lim_{r\to0}\left[G_B(\mathbf r)-\frac1r\right].
$$

The uniform term in the Poisson equation is the compensating background.
The correction cancels the regular image/background contribution at the
origin. Periodic images still change the potential at finite separation. Consequently, a
long-range density-direct entry referenced to this periodic potential can
be negative even though the isolated kernel `C/r` is positive. Such a result
calls for identifying the intended electrostatic operator; it does not by
itself demonstrate an error in the underlying DFT orbitals.

In native VASP 6.6.1, `HFRCUT=0` selects that periodic/convergence treatment;
positive values specify a spherical radius. The documented `HFRCUT=-1`
chooses a volume-matched sphere:

$$
K_{\mathrm{sph}}(r)=\frac{C}{r}\Theta(R_c-r),\qquad
R_c=\left(\frac{3N_k\Omega}{4\pi}\right)^{1/3}.
$$

The radius formula assumes the full uniform mesh without additional mesh
reductions. The native code sends all negative values, including an unpatched
`-2`, through its automatic-sphere branch. The WS meaning of `-2` is introduced
by this patch.

#### 3.2 Why a two-dimensional moiré geometry is demanding

Vacuum sufficient to suppress orbital overlap between layers can still be
short compared with the in-plane moiré length. The relevant aspect ratio is
that of **`N_i a_i`**, not only that of the primitive simulation cell. With one
normal k point, the normal BvK period remains `Lz`; refining the in-plane k mesh
only enlarges the in-plane extent.

A single spherical radius cannot follow these very different lengths. For
example, an approximately orthogonal `41 × 39 × 29 Å` cell with a `6 × 6 × 1` mesh has BvK lengths
near `246 × 233 × 29 Å`. Its automatic spherical radius is about `73.5 Å`,
whereas the normal WS half-width is about `14.5 Å` and the in-plane half-widths
are about `117–123 Å`. The sphere overlaps repeated normal cells while ending
well inside the useful in-plane range. Its sampled Fourier components describe
a periodized spherical kernel; overlap can therefore retain slab-image
contributions.

This shape problem and the periodic/FSG reference above are distinct. The
negative tails of a calculation using `HFRCUT=0` should not be attributed to a
spherical cutoff that was not selected.

#### 3.3 The WS construction

The WS domain follows the actual BvK lattice metric:

$$
W_B=\{\mathbf s:\ |\mathbf s|\leq|\mathbf s-\mathbf S|
\text{ for every }\mathbf S\in\Lambda_B\},
\qquad
K_{\mathrm{WS}}(\mathbf s)=\frac{C}{|\mathbf s|}\chi_{W_B}(\mathbf s).
$$

$$
\widetilde K_{\mathrm{WS}}(\mathbf Q)
=\int_{W_B}d^3s\,\frac{C\,e^{-i\mathbf Q\cdot\mathbf s}}{|\mathbf s|}.
$$

For orthogonal vectors it is the box with half-widths `|A_i|/2`; oblique cells
use the true nearest-image geometry. All Fourier components, including the
finite head, come from this same truncated kernel. The old Ewald/FSG head is
absent in this mode.

This is the WS-supercell construction of
[Sundararaman and Arias (2013)](https://doi.org/10.1103/PhysRevB.87.165122).
Their analysis includes partially periodic systems. The point of the WS
shape is to match the repetition lattice: it is useful for both bulk and
layered systems. The paper's exchange-energy convergence result has its own
localization assumptions; accuracy of an individual long-range Wannier
integral still depends on orbital tails, spatial range and numerical resolution.

| Input | Definition |
|---|---|
| `HFRCUT=-2` | Enable this patch's WS representation for the Wannier-integral workflow |
| `WINT_BOUNDARY=BVK_WS` | Full BvK WS domain; default boundary for explicit BARE/YUKAWA |
| `WS_NGRID` | Three-component kernel FFT grid; automatic by default |
| `WS_TOL` | Internal kernel construction/quadrature control; default `1e-8` |

With all `WINT_*` tags omitted, `HFRCUT=-2` preserves the earlier bare-WS
evaluator. An explicit `WINT_KERNEL=BARE` selects the complete reconstructed-PAW
density route described below. They use the same bare law through different
finite-basis evaluations.

<a id="en-slab"></a>

### 4. `SLAB_WS`: a finite-height representation of layered electrostatics

#### 4.1 Why keep two normal coordinates?

A metal plate or a polarizable sheet occupies a particular height. Two pairs
of charges with the same displacement can experience different screening if
one pair lies closer to that sheet or plate. Their interaction is therefore
`V(rho−rho'; z,z')`, with separate source and observation heights. A scalar
multiplier indexed only by a three-dimensional momentum cannot express this
general dependence.

`SLAB_WS` combines a finite in-plane WS range with such a height-resolved Green
function. It is the name of the construction defined here. Published Coulomb
cutoff and layered Green-function methods provide its ingredients; there is
no claim that a paper presents this exact combination under that name.

#### 4.2 Specify the physical Green function, then its finite representation

First choose the infinite-plane interaction `V_infinity(s;z,z')` appropriate to
the environment: free Coulomb, a polarizable sheet, or grounded plates. Let
`W_parallel` be the WS cell generated by the two in-plane BvK vectors. Then
define the actual finite-domain operator:

$$
I=[z_0-L_z/2,z_0+L_z/2],\qquad
K_W(\boldsymbol\rho,z;\boldsymbol\rho',z')
=\chi_{W_\parallel}(\boldsymbol\rho-\boldsymbol\rho')
V_\infty(\boldsymbol\rho-\boldsymbol\rho';z,z'),
\qquad z,z'\in I.
$$

The planar WS window limits relative in-plane separations. Each density is
represented over **its own** height coordinate in the same interval `I`.
The interval endpoints are the boundaries of the density representation, not
conducting walls. Free and sheet kernels use an unbounded normal Green
function. In the gate model the actual grounded walls are at `z0±d`.

This differs from the published relative-height slab cutoff
`Theta(zc−|z−z'|)/|r−r'|` of
[Ismail-Beigi (2006)](https://doi.org/10.1103/PhysRevB.73.233103).
That construction cuts a function of the displacement. The present one selects
a continuous normal representative of each density and permits a Green
function depending on both heights.

A Fourier transform in the plane makes another distinction explicit. If
`G_infinity(p;z,z')` is the analytic infinite-plane Green function, the finite
WS operator is

$$
\widetilde V_W(\mathbf q;z,z')
=\int\frac{d^2p}{(2\pi)^2}\,
\widetilde\chi_{W_\parallel}(\mathbf q-\mathbf p)
G_\infty(\mathbf p;z,z').
$$

Multiplication by a real-space window becomes convolution in momentum.
Evaluating the infinite-plane analytic expression at a single `q` would omit
this window. The implementation integrates the windowed operator itself.

#### 4.3 Normal Fourier–Galerkin projection

The WAVECAR-derived pair densities are available as Fourier coefficients.
A normal Fourier basis can represent them inside `I` while the electrostatic
Green function remains open. In centered coordinates, define

$$
x=z-z_0,\quad g_m=\frac{2\pi m}{L_z},\quad
\phi_m(x)=\frac{e^{ig_mx}}{\sqrt{L_z}},
$$

$$
T_{mn}(\mathbf q)=\frac1{L_z}
\int_{-L_z/2}^{L_z/2}\!dx\int_{-L_z/2}^{L_z/2}\!dx'\,
e^{-ig_mx}\,
\widetilde V_W(\mathbf q;z_0+x,z_0+x')\,e^{ig_nx'}.
$$

$$
K^{\mathrm{global}}_{mn}(\mathbf q)
=e^{i(g_n-g_m)z_0}T_{mn}(\mathbf q).
$$

The last line expresses the matrix in the global-coordinate Fourier convention
of the supplied density. It is why `WINT_Z0` has an actual role in the integral,
including for bare Coulomb.

Even the free kernel, which depends only on `z−z'`, generally produces a dense
matrix on a finite open interval. Convolution over that interval is not the
circular convolution diagonalized by a periodic Green function. Fourier basis
functions provide a representation; they do not impose conducting or periodic
electrostatic boundary conditions on the chosen Green function.

With `C` included, the planar Green function has units `eV Å²`, and this normal
matrix has units `eV Å³`. The native Fourier/volume normalization converts the
projected contraction to energies. The purely geometric matrix obtained after
removing `C` has units `Å²`.

The computational separation follows the identity

$$
\frac1{\sqrt{\rho^2+(x-x')^2}}
=\frac2{\sqrt\pi}\int_0^\infty dt\,
e^{-t^2\rho^2}e^{-t^2(x-x')^2}.
$$

It separates a plane integral over the WS window from normal matrix elements.
The latter are evaluated as one-dimensional moments. Sheet-induced and gate
image terms admit corresponding representations, so the calculation retains
orbital thickness without sampling a singular interaction on an arbitrary
extra normal grid.

#### 4.4 Physical assumptions and the choice of reference plane

`WINT_Z0` is an absolute Cartesian coordinate along the layer normal, in Å.
For bare/Yukawa it selects the numerical interval center. For RK_SHEET it also
locates the ideal polarizable plane; for DUAL_GATE it locates the plate midpoint.
A geometric layer center is a useful starting point, but choosing an effective
polarization plane or an actual device geometry is a physical model choice.

Choose `I` to contain the layer and its significant density tails, placing its
ends in vacuum. Shifting those ends within adequate vacuum should change a
bare integral only through tail and numerical errors. A material crossing the
periodic coordinate origin should first be represented as one continuous
layer. The arithmetic center of wrapped fractional coordinates can be wrong.

The finite plane window is also a range choice. It does not by itself give the
infinite-area interaction at every stored lattice vector. WS faces, normal
vacuum, Fourier resolution and physical screening parameters require separate
assessment. `BVK_WS` remains meaningful for layers; the two finite-domain bare
operators can agree in the interior when all important normal separations lie
within the BvK normal cutoff and density tails are small. Their complete
finite-basis matrices are not generally identical.

| Input / condition | Meaning |
|---|---|
| `WINT_BOUNDARY=SLAB_WS` | In-plane BvK WS window plus open-normal density interval |
| `WINT_Z0` | Required interval center for finite-height models; also the sheet/gate reference plane when applicable |
| `WINT_NORMAL` | Optional Cartesian normal, default along the third lattice vector |
| Lattice and k mesh | Third vector perpendicular to the plane; one normal k point; oblique in-plane cells are supported |
| `WS_NGRID` | First two components refine the plane kernel FFT; third component is unused by this boundary |

RK_SHEET, RK2D and DUAL_GATE currently use this boundary explicitly. BARE and
YUKAWA support either boundary. RK2D integrates out height and does not use Z0.

<a id="en-kernels"></a>

### 5. Choosing a physical interaction

All models below are static operators projected onto the same chosen Wannier
basis. The dielectric constant `epsilon` is relative permittivity; `epsilon_0`
is the vacuum constant. A screening length is a physical input, distinct from
an FFT convergence parameter.

#### 5.1 `BARE`: the repulsion of the Wannier charge clouds

$$
K_0(\mathbf r,\mathbf r')=\frac{C}{|\mathbf r-\mathbf r'|},
\qquad C=\frac{e^2}{4\pi\epsilon_0}.
$$

This is the vacuum reference for an unscreened interacting model. Localization
sets its onsite scale: concentrating a normalized orbital generally increases
its electrostatic self-interaction. Well-separated localized densities approach
a point-charge `C/distance` interaction away from finite-domain boundaries.
A PBE Wannier orbital changes the shape of that charge cloud; it does not
insert dielectric screening into the bare kernel.

| Keyword | Meaning |
|---|---|
| `WINT_KERNEL=BARE` | Explicit full-PAW bare projection |
| `WINT_BOUNDARY` | BVK_WS by default, or SLAB_WS |
| `WINT_Z0` | Required only with SLAB_WS; sets its numerical interval |

For explicit BARE, omit `WINT_EPS` and the screening lengths. For a comparison
with layered screened models, use the same SLAB_WS domain and Z0 as their
finite-height reference.

#### 5.2 `YUKAWA`: homogeneous three-dimensional screening

A mobile screening charge in a homogeneous medium can be represented by the
modified Poisson equation:

$$
K_Y(r)=\frac{C}{\epsilon}\frac{e^{-\kappa r}}r,
\qquad
(-\nabla^2+\kappa^2)K_Y
=\frac{4\pi C}{\epsilon}\delta(\mathbf r),
\qquad
\widetilde K_Y(\mathbf Q)=\frac{4\pi C}{\epsilon(Q^2+\kappa^2)}.
$$

The additional inverse length `kappa` suppresses the distant field with
`lambda=1/kappa`. The functional form originates with
[Yukawa (1935)](https://doi.org/10.11429/ppmsj1919.17.0_48); its use here is the
static electronic-screening model defined by this equation. It is useful for
exploring how interaction range changes a low-energy model. Metallic plates
instead impose spatial boundary conditions and generate a different kernel.

| Keyword | Meaning |
|---|---|
| `WINT_KERNEL=YUKAWA` | Exponentially screened radial interaction |
| `WINT_KAPPA` | Required inverse screening length, Å⁻¹; finite and nonnegative |
| `WINT_EPS` | Homogeneous relative permittivity, default 1; finite and positive |
| `WINT_BOUNDARY` | BVK_WS or SLAB_WS; the latter also requires Z0 |

At `kappa=0` the result is bare/epsilon on the same domain. With `epsilon=1`
it equals explicit BARE. Increasing kappa emphasizes shorter length scales,
so the density resolution must also be checked.

#### 5.3 `RK_SHEET`: a polarizable sheet with finite-height Wannier charges

A two-dimensional insulating layer polarizes mainly within the material while
its electric field extends into the surrounding three-dimensional dielectric.
This makes screening depend on in-plane wavelength. An ideal-sheet model
captures that response without collapsing the two charge clouds onto one plane.

Place the sheet at `z0` and define `r0` through the following Poisson problem:

$$
\left[-\partial_z^2+q^2+2r_0q^2\delta(z-z_0)\right]
G_s(q;z,z')=\frac{4\pi C}{\epsilon}\delta(z-z').
$$

Solving it with a decaying free-space solution on either side gives

$$
G_s(q;z,z')=\frac{2\pi C}{\epsilon q}
\left[e^{-q|z-z'|}
-\frac{r_0q}{1+r_0q}
 e^{-q(|z-z_0|+|z'-z_0|)}\right],\qquad q>0.
$$

The first term is the background Coulomb field; the second is generated by
induced sheet charge. Continuity of the potential and the derivative jump

$$
\partial_zG_s(z_0^+,z')-\partial_zG_s(z_0^-,z')
=2r_0q^2G_s(z_0,z'),\qquad z'\ne z_0.
$$

specify the sheet response. This finite-height expression is derived here from
the ideal-sheet problem. Its on-sheet limit is the Rytova–Keldysh interaction
reviewed and derived by [Cudazzo, Tokatly and Rubio (2011)](https://doi.org/10.1103/PhysRevB.84.085406),
following the earlier work of Rytova and Keldysh [7–9](#references--参考文献).

For charges on the sheet, `r0 q/(1+r0 q)` increasingly suppresses
short-wavelength fields. Away from the sheet, the exponential also suppresses
the response according to the charge heights. For two charges on the same side
away from the sheet, the relative sheet correction vanishes at sufficiently
large q. At large in-plane separation the
`C/(epsilon rho)` tail survives. For charges coinciding **away from the sheet**,
the local three-dimensional Coulomb singularity remains; the logarithmic RK
singularity belongs to the on-sheet limit. This model keeps charge thickness
but represents the polarizable medium as one ideal plane, rather than a film
with a spatially varying dielectric tensor.

| Keyword | Meaning |
|---|---|
| `WINT_KERNEL=RK_SHEET` | Finite-height ideal polarizable-sheet Green function |
| `WINT_BOUNDARY=SLAB_WS` | Required boundary representation |
| `WINT_R0` | Required effective screening length, Å; finite and nonnegative |
| `WINT_Z0` | Required sheet position and integration-interval center, Å |
| `WINT_EPS` | Uniform background relative permittivity, default 1; finite and positive |

Here `r0` is defined by the displayed equation. Converting a 2D polarizability
into r0 requires consistent units and the same background-dielectric convention.
At `r0=0` the sheet disappears and the finite-height interaction is bare/epsilon.

#### 5.4 `RK2D`: the strict zero-thickness approximation

If the charge thickness is negligible on the length scales being studied, each
pair density can first be integrated over height and represented on the sheet:

$$
G_{2D}(q)=\frac{2\pi C}{\epsilon q(1+r_0q)},\qquad
K_{2D}(\rho)=\frac{\pi C}{2\epsilon r_0}
\left[H_0(\rho/r_0)-Y_0(\rho/r_0)\right],
$$

$$
\rho_{2D}(\boldsymbol\rho)=\int_I dz\,\rho_{3D}(\boldsymbol\rho,z).
$$

The displayed real-space special-function expression applies for r0>0; zero
is taken through the limit below. `H0` and `Y0` are the Struve function and
Bessel function of the second kind.
The asymptotic forms show what a polarizable sheet does:

$$
K_{2D}(\rho)\simeq
\begin{cases}
\displaystyle\frac{C}{\epsilon r_0}
\left[\ln\!\left(\frac{2r_0}{\rho}\right)-\gamma\right],&\rho\ll r_0,\\[6pt]
\displaystyle\frac{C}{\epsilon\rho},&\rho\gg r_0.
\end{cases}
$$

The short-distance logarithm reflects the sheet response; the distant Coulomb
tail reflects electric field propagating outside the layer. This is the same
ideal-sheet screening law as the on-sheet limit above, applied to projected
planar densities. It provides a comparison for neglecting Wannier thickness.
Even if thickness is unimportant for distant interactions, onsite matrix
elements probe large momenta where that approximation can matter. Extending
a long-wavelength sheet model into the microscopic onsite range is a model
assumption to assess for the chosen orbitals.

| Keyword | Meaning |
|---|---|
| `WINT_KERNEL=RK2D` | RK operator acting on height-integrated pair densities |
| `WINT_BOUNDARY=SLAB_WS` | Required in-plane finite-domain representation |
| `WINT_R0` | Required screening length, Å; finite and nonnegative; same convention as RK_SHEET |
| `WINT_EPS` | Background relative permittivity, default 1; finite and positive |
| `WINT_Z0` | Omit: the normal coordinate has been integrated out |

At `r0=0`, this gives planar `C/(epsilon rho)`, rather than the finite-height
BARE operator. The electric field remains three-dimensional even though the
charge/polarization model is two-dimensional.

#### 5.5 `DUAL_GATE`: screening by two grounded metal plates

Gates screen through induced surface charge. With two perfect plates at
`z=z0−d` and `z=z0+d`, define `x=z−z0`. The appropriate electrostatic problem is

$$
(-\partial_x^2+q^2)G_g(q;x,x')
=\frac{4\pi C}{\epsilon}\delta(x-x'),
\qquad G_g(q;\pm d,x')=0.
$$

Its finite-height Dirichlet Green function is

$$
G_g(q;x,x')=\frac{4\pi C}{\epsilon q}
\frac{\sinh[q(x_<+d)]\sinh[q(d-x_>)]}{\sinh(2qd)},
\quad x_<=\min(x,x'),\quad x_>=\max(x,x').
$$

$$
G_g(q;0,0)=\frac{2\pi C}{\epsilon q}\tanh(qd),
\qquad G_g(0;0,0)=\frac{2\pi Cd}{\epsilon}.
$$

This follows from the plate Green function and image construction in
[Takae and Onuki (2013)](https://doi.org/10.1063/1.4821085).
The parallel-plate model also appears in the interacting-electron setting of
[Throckmorton and Vafek (2012)](https://doi.org/10.1103/PhysRevB.86.115447).
Here `d` always means the distance from the midpoint to **each** plate.

For point charges in the midplane, the equivalent images and transverse modes
make the spatial screening especially transparent:

$$
K_g(\rho;0,0)=\frac{C}{\epsilon}
\sum_{n=-\infty}^{\infty}
\frac{(-1)^n}{\sqrt{\rho^2+(2nd)^2}}
=\frac{2C}{\epsilon d}\sum_{j=0}^\infty
K_0\!\left(\frac{(2j+1)\pi\rho}{2d}\right).
$$

The lowest transverse mode produces a distant tail proportional to
`exp[−pi rho/(2d)]/sqrt(rho)`, with scale `2d/pi`. A single Yukawa length therefore
does not reproduce the complete gate interaction. The finite-height Green
function also distinguishes charges closer to either plate.

| Keyword | Meaning |
|---|---|
| `WINT_KERNEL=DUAL_GATE` | Two grounded plates, retaining both charge heights |
| `WINT_BOUNDARY=SLAB_WS` | Required boundary representation |
| `WINT_GATE_D` | Required midpoint-to-each-plate distance d in Å; plate spacing is 2d |
| `WINT_Z0` | Required plate midpoint and integration-interval center, Å |
| `WINT_EPS` | Homogeneous spacer relative permittivity, default 1; finite and positive |

The current interval must lie between the plates, requiring `d>Lz/2`.
This is a restriction of the chosen integration domain, not a rule that a
physical gate must sit outside a DFT vacuum cell. At large d the operator
approaches finite-height bare/epsilon. The individual two-electron integral
contains no extra factor 1/2; the many-body Hamiltonian supplies its pair
counting convention. One-electron image potentials are a separate contribution
to its hopping/onsite Hamiltonian.

<a id="en-use"></a>

### 6. PAW projection and practical use

#### The same orbital density must meet the chosen operator

PAW wavefunctions contain smooth plane-wave and atom-centered reconstruction
parts. The interaction must act consistently on both. In the explicit WINT
route the complete fixed pair density is reconstructed as

$$
\rho_{ij}=\rho^{\mathrm{PW}}_{ij}
+\sum_a\left(\rho^{\mathrm{AE}}_{ij,a}
-\rho^{\mathrm{PS}}_{ij,a}\right).
$$

The operator is then contracted with this density, including smooth–atomic
and interatomic cross terms. This uses the supplied PAW partial waves and
projectors; a separate all-electron solid calculation is unnecessary.
The PAW principle and its VASP implementation are described in
[Blöchl (1994)](https://doi.org/10.1103/PhysRevB.50.17953) and
[Kresse and Joubert (1999)](https://doi.org/10.1103/PhysRevB.59.1758).

This common density construction lets the physical Green function change
without inheriting an operator-specific bare local correction. The older
omitted-tag bare route retains its native compensation/local-correction
construction. Their finite-cutoff agreement is a numerical question; the
explicit route is the common reference for the new screened models.

#### From a finalized Wannier basis to interaction files

1. Obtain a converged electronic structure and consistent WAVECAR/k mesh.
2. Use compact SCDM for the initial Wannier projection.
3. Finalize the Wannier gauge, including external SAWF if required, and export
   its VASP-compatible WANPROJ.
4. Run the fixed-orbital `ALGO=2E4WA` calculation with the chosen operator.

The supported interaction path uses CPU NCL, `ISPIN=1` and `NCORE=1`, including
noncollinear SOC spinors. A typical finite-height gate example, for a material
reference plane at Cartesian z=0, adds:

```text
ALGO = 2E4WA
HFRCUT = -2
NCORE = 1
WINT_KERNEL = DUAL_GATE
WINT_BOUNDARY = SLAB_WS
WINT_Z0 = 0.0
WINT_GATE_D = 100.0
```

Use the actual reference plane of the structure. Here 100 Å is 10 nm to each
gate, an illustrative geometric parameter. Preserve the fixed orbital basis
and PAW datasets when comparing kernels. VASP writes VIJKL and VRijkl in its
native index/file convention; model parameters are recorded in their headers
and OUTCAR. The unchanged one-electron SAWF Hamiltonian and these interactions
must describe the same Wannier basis.

For an fRG model, epsilon, r0 and kappa specify a chosen environmental or
partially screened interaction. Low-energy target-band screening that fRG
will generate should be separated from any screening already included in the
input. The current scalar epsilon describes a homogeneous isotropic medium;
an anisotropic substrate or a multilayer device needs an appropriate effective
model or a more general Poisson Green function. Sheet polarization and gates
can coexist physically; the current selectors describe separate models rather
than their combined Poisson problem.

<a id="en-convergence"></a>

### 7. What must be converged for a physical interpretation?

- **Orbital model:** band selection, Wannier gauge and the one-electron band
  representation determine which charge clouds interact.
- **Normal vacuum and placement:** contain the significant density tails in I.
  For bare, moving the boundaries within vacuum checks the finite-domain
  representation; moving a sheet or gate relative to a fixed layer changes
  the physical environment.
- **Spatial range:** the BvK k mesh fixes the WS range. Increasing WS_NGRID at
  fixed k mesh refines the integral but does not enlarge that range.
- **Kernel resolution:** refine WS_NGRID. WS_TOL controls construction and
  internal quadrature; it is not a final error bar on every matrix element.
- **Density resolution:** converge VCUTOFF together with NGX/NGY/NGZ. Short-range
  and onsite interactions probe high-momentum density components. Raising
  VCUTOFF alone can encounter a density-grid ceiling.

When an enlarged density grid is incompatible with the old CHGCAR,
`ICHARG=12` has been checked for this fixed-WAVECAR interaction projection.
Its companion one-electron/FCIDUMP quantities do not replace the established
SAWF Hamiltonian. Physical screening parameters remain distinct from these
numerical refinements.

Independent physical Green-function checks, normal-matrix tests and native
projection/writer comparisons validate the implemented finite-basis operator.
Material-specific interaction accuracy still follows from the convergence
questions above. For a positive pointwise kernel and nonnegative exact densities,
density-direct entries should be positive; general transition-density entries
have no such sign requirement.

<a id="en-install"></a>

### Installation and updating the patch

The current source profiles support VASP 6.6.1 and 6.6.0; native local execution
uses 6.6.1 CPU NCL. A valid VASP license and your own VASP source are required.
For clean vendor or recognized SCDM-only sources, use the compiler/MPI/
Wannier90 setup of your working installation:

```bash
sha256sum -c manifest.sha256
VASP_WINT_SRC=/absolute/path/to/vasp-source
bash apply_patch.sh --status "$VASP_WINT_SRC"
bash apply_patch.sh --check "$VASP_WINT_SRC"
bash apply_patch.sh --apply "$VASP_WINT_SRC"
make -C "$VASP_WINT_SRC" ncl DEPS=1 -j4
```

For an upgrade, **first use the currently installed version's matching script
and patch** to restore the vendor source:

```bash
bash /path/to/current-package/apply_patch.sh --reverse "$VASP_WINT_SRC"
```

Then replace the package by the new version, run its check/apply operations and
rebuild. Reverse before overwriting the old pair or running `git pull` in its directory:
`--reverse` verifies the postimages of its own payload. The repository needs no separate collection of
historical rollback packages. Recognized clean vendor or SCDM-only sources can
use the current installer directly. Source backups are retained by the script.

Deployment requires the adjacent `apply_patch.sh` and unified `.patch` file.
`SOURCE_IDENTITIES.json` records exact source profiles; `manifest.sha256`
checks package integrity. The historical patch filename covers all features
above. The repository distributes a third-party patch, not VASP binaries,
PAW datasets or a complete VASP source distribution.

## 中文说明

### 1. 物理目标：从 Bloch 能带到相互作用 Wannier 模型

低能模型既需要合适的轨道基底，也需要与环境相符的电子相互作用。moiré 体系可能
用少数几条费米能附近的能带描述低能物理，但平面波计算同时需要数千条能带。
Wannier 函数把这些 Bloch 态转换为局域轨道，从而得到跃迁和轨道电荷分布之间的
相互作用矩阵元。

对与自旋无关的标量相互作用，可用下面的双线性形式表示投影。密度直接项是两个
归一化 Wannier 电荷云之间的静电相互作用：

$$
\mathcal B_K[\rho,\sigma]
=\int d^3r\,d^3r'\,\rho^*(\mathbf r)
K(\mathbf r,\mathbf r')\sigma(\mathbf r'),\qquad
D_{mn}(\mathbf R)=\mathcal B_K[|w_{m0}|^2,|w_{n\mathbf R}|^2].
$$

轨道的大小及内部结构决定 onsite 排斥的尺度；off-center 项给出空间作用范围。
将这些尺度与带宽比较，有助于判断低能模型的相互作用强度。一般的跃迁密度矩阵元
还包含交换及轨道变化过程，不要求全部为正。文件指标约定见
[VASP VRijkl 说明](https://vasp.at/wiki/VRijkl)。

本补丁分别处理三个问题：

1. **Band-compact SCDM：** 降低大能带计算中初始 Wannier 投影的内存需求。
2. **`HFRCUT=-2`：** 用与有限 k 网格对应超晶格匹配的 WS 区域表示裸库仑相互作用。
3. **层状相互作用：** 保留 Wannier 密度的法向厚度，并引入理想极化片或金属栅极。

轨道基底、物理核和有限积分区域是不同的选择。固定 WAVECAR 与最终 WANPROJ 后
改变 kernel，比较的是同一轨道基底中的不同相互作用。PBE 可以提供这组轨道；
随后选用屏蔽核，并不会把原来的 PBE 计算变成杂化泛函计算。

<a id="zh-scdm"></a>

### 2. Band-compact SCDM：让工作空间对应真正保留的能带

#### 动机和原理

SCDM 从密度矩阵或经过谱过滤的投影算符中选择具有代表性的实空间列，借助带列主元
的 QR 分解确定一组锚点，并由这些列构造局域轨道的初始基底。原始 SCDM、周期
SCDM-k 及纠缠能带的处理见参考文献 [1–3](#references--参考文献)。

谱过滤让希望保留的态获得更大权重：

$$
P_f(\mathbf k)=\sum_{n\in\mathcal A}
f(E_{n\mathbf k})\,|u_{n\mathbf k}\rangle\langle u_{n\mathbf k}|,
\qquad
M_{\mathrm{SCDM}}\propto
N_{\mathrm{grid}}N_{\mathrm{spinor}}N_{\mathrm{active}}.
$$

这里的权重用于定义 Wannier 子空间，不必等于热力学占据数。例如 ERFC 过滤可以
在某一能量附近平滑降低高能态的权重：

$$
f(E)=\frac12\mathrm{erfc}\!\left(\frac{E-\mu}{\sigma}\right).
$$

在这个例子中，`CUTOFF_MU` 是过滤中心，`CUTOFF_SIGMA` 是能量宽度。它们用于
选择子空间；过滤中心不必等于体系的化学势。

内存问题来自：Wannier90 已排除了大部分能带，SCDM 工作矩阵仍可能按 VASP 的
完整能带数分配。moiré 大晶胞同时增大实空间网格和总能带数，使这种额外分配尤其
昂贵。本补丁把已有 `exclude_bands` 掩码用于 SCDM 的波函数、能量、过滤权重和
变换矩阵，按原始能带编号的升序组织保留态，使 AMN 行与原始能带保持对应。

主导工作矩阵的维度由 `NBANDS` 缩小为 `N_active`。这是对现有 SCDM 存储及指标
组织的优化；全局波函数、WAVECAR 和原有 MMN 临时矩阵仍保留原来的维度，因此
总 VASP 内存还包括其他部分。

#### 怎样使用

使用常规 VASP/Wannier90 接口，并在 Wannier90 输入中指定选带。若孤立子空间
恰好包含与 Wannier 函数相同数量的能带，单位权重对应自然的子空间投影：

```text
LSCDM = .TRUE.
LWANNIER90 = .TRUE.
NUM_WANN = 8
CUTOFF_TYPE = none
```

此例需要配套的 `exclude_bands` 保留八条合适的能带。纠缠子空间则保留更大的窗口，
并显式指定过滤方式。

| 输入 | 在紧凑子空间中的作用 |
|---|---|
| Wannier90 `exclude_bands` | 从 SCDM 工作空间排除的原始 VASP 能带 |
| `NUM_WANN` | 需要构造的 Wannier 函数数目 |
| `CUTOFF_TYPE=none` | 保留态均使用单位权重 |
| `CUTOFF_MU`、`CUTOFF_SIGMA` | 显式谱过滤参数；ERFC 的中心与宽度定义见上式 |

有选带且保留带数多于 Wannier 数时，ERFC/Gaussian 需同时给出两个参数，Fermi
需给出 `CUTOFF_MU`，`num_wann` 需给出 `CUTOFF_SIGMA`。其他过滤形式保留本版本
VASP 的约定，不能直接套用上面的 ERFC 宽度公式。若保留带数恰好为 `NUM_WANN`，
不设置过滤也会使用单位权重。

SCDM 提供初始投影；最终 SAWF 局域化及 WANPROJ 转换由外部工具完成。本补丁在
相互作用阶段读取最终的 WANPROJ。

<a id="zh-ws"></a>

### 3. `HFRCUT=-2`：为什么需要随晶格形状变化的库仑截断

#### 3.1 原生 VASP 的处理对应什么物理对象

裸库仑核的倒空间形式为 `4πC/Q²`，在 `Q=0` 有可积奇点。有限 k 网格需要为该项
指定处理方法。原生 VASP 对规则、未屏蔽网格使用 probe-charge Ewald 修正，其他
支持的情况使用相应高斯辅助函数处理；也可通过 `HFRCUT` 选择球形截断。定义见
[HFRCUT](https://vasp.at/wiki/HFRCUT) 和
[Coulomb singularity](https://vasp.at/wiki/Coulomb_singularity)。

对于没有额外约化的完整均匀网格，实空间周期性对应 Born–von Karman 超晶格：

$$
\mathbf A_i=N_i\mathbf a_i,\qquad
N_k=N_1N_2N_3,\qquad V_B=N_k\Omega.
$$

BvK 是 Bloch 态离散采样对应的周期性；Ewald 是这个周期体系上的静电求和及收敛
处理。要把投影相互作用解释为某个物理模型，还需要明确其电势和边界的定义。

对规则、未屏蔽情形，周期 Green 函数及局域 Madelung 修正可写成

$$
-\nabla^2G_B(\mathbf r)
=4\pi\left[\sum_{\mathbf S\in\Lambda_B}
\delta(\mathbf r-\mathbf S)-\frac1{V_B}\right],
\qquad \langle G_B\rangle=0.
$$

$$
K_{\mathrm{per}}(\mathbf r)=C\,[G_B(\mathbf r)-v_M],
\qquad
v_M=\lim_{r\to0}\left[G_B(\mathbf r)-\frac1r\right].
$$

Poisson 方程中的均匀项对应补偿背景。该修正抵消原点处的规则镜像/背景贡献，
但周期镜像仍会改变有限距离的电势。因此，以这一周期电势为参考的长程密度直接项
可以出现负值，即使孤立的 `C/r` 处处为正。此时应先辨认实际计算的静电算符，
不能仅凭负值认定 DFT 轨道本身错误。

原生 VASP 6.6.1 中，`HFRCUT=0` 选择上述周期/收敛处理；正值指定球形半径。
文档中的 `HFRCUT=-1` 选择与 BvK 体积相同的球：

$$
K_{\mathrm{sph}}(r)=\frac{C}{r}\Theta(R_c-r),\qquad
R_c=\left(\frac{3N_k\Omega}{4\pi}\right)^{1/3}.
$$

该半径公式对应无额外网格约化的完整均匀网格。原生源码会把所有负数，包括未打
补丁时的 `-2`，送入自动球形截断分支。这里的 WS 含义由本补丁赋予 `-2`。

#### 3.2 二维 moiré 几何为什么特别敏感

真空足以使层间轨道重叠很小，并不意味着法向晶胞尺寸大于面内 moiré 尺度。
更关键的是 **`N_i a_i` 的比例**，而不只是 DFT 晶胞自身的比例。法向仅一个 k 点
时，BvK 的法向周期仍是 `Lz`；增加面内 k 点只会扩展面内范围。

一个球形半径难以同时匹配这些长度。例如约 `41 × 39 × 29 Å` 的近似正交晶胞配合
`6 × 6 × 1` 网格，其 BvK 尺寸约为 `246 × 233 × 29 Å`。自动球形半径约
`73.5 Å`，法向 WS 半宽却只有约 `14.5 Å`，面内半宽则约为 `117–123 Å`。
球在法向跨越多个重复晶胞，同时过早截掉面内的长距离。离散 Fourier 分量表示
周期重复的球形截断核，因此这种重叠可以保留层间镜像贡献。

这种形状不匹配，与前面的周期/FSG 电势参考问题应分别讨论：使用 `HFRCUT=0`
得到的负尾部，不能归因于该计算没有选择的球形截断。

#### 3.3 WS 构造怎样改变它

WS 区域由实际 BvK 晶格度量确定：

$$
W_B=\{\mathbf s:\ |\mathbf s|\leq|\mathbf s-\mathbf S|
\text{ for every }\mathbf S\in\Lambda_B\},
\qquad
K_{\mathrm{WS}}(\mathbf s)=\frac{C}{|\mathbf s|}\chi_{W_B}(\mathbf s).
$$

$$
\widetilde K_{\mathrm{WS}}(\mathbf Q)
=\int_{W_B}d^3s\,\frac{C\,e^{-i\mathbf Q\cdot\mathbf s}}{|\mathbf s|}.
$$

正交晶格中，它是各方向半宽为 `|A_i|/2` 的长方体；斜晶格使用真实最近镜像几何。
所有 Fourier 分量，包括有限的零模，都来自同一个截断核。此模式不再加入旧的
Ewald/FSG head。

这采用 [Sundararaman–Arias (2013)](https://doi.org/10.1103/PhysRevB.87.165122)
的 WS 超胞构造，其讨论包含部分周期体系。WS 的作用是匹配重复晶格，因此同样
适用于体材料和层状材料。论文中交换能的收敛结论具有相应的局域化假设；单个
长程 Wannier 矩阵元的精度，仍取决于密度尾部、作用范围和数值分辨率。

| 输入 | 定义 |
|---|---|
| `HFRCUT=-2` | 在 Wannier 积分工作流中启用本补丁的 WS 表示 |
| `WINT_BOUNDARY=BVK_WS` | 完整 BvK WS 区域；显式 BARE/YUKAWA 的默认边界 |
| `WS_NGRID` | 三分量核 FFT 网格；默认自动选择 |
| `WS_TOL` | 核构造及内部积分控制，默认 `1e-8` |

省略全部 `WINT_*` 标签时，`HFRCUT=-2` 保留原来的 bare-WS 求值路径。显式设置
`WINT_KERNEL=BARE` 则使用下文的完整重建 PAW 密度。二者采用相同裸库仑律，但
有限基组求值方式不同。

<a id="zh-slab"></a>

### 4. `SLAB_WS`：保留层状静电问题中的两个高度坐标

#### 4.1 为什么需要分别保留 z 和 z'？

金属板和极化片处在确定的高度。两对电荷即使具有相同相对位移，只要到板或极化片
的距离不同，屏蔽就可能不同。因此相互作用一般为 `V(rho−rho'; z,z')`，分别依赖
源与观察点高度。只按三维动量作标量乘法，不能表达这种一般的高度依赖。

`SLAB_WS` 把有限的面内 WS 范围与高度分辨的 Green 函数结合起来。这个名字指本
项目在此定义的构造。已发表的库仑截断和层状 Green 函数方法提供其组成部分，
具体组合则在这里明确给出，不能说某篇论文原样提出了这一同名方法。

#### 4.2 先定义物理环境，再定义有限表示

先根据环境选择无限平面上的相互作用 `V_infinity(s;z,z')`，例如自由库仑、极化片
或接地金属板。令 `W_parallel` 为面内两根 BvK 矢量生成的 WS 区域，实际有限算符
定义为

$$
I=[z_0-L_z/2,z_0+L_z/2],\qquad
K_W(\boldsymbol\rho,z;\boldsymbol\rho',z')
=\chi_{W_\parallel}(\boldsymbol\rho-\boldsymbol\rho')
V_\infty(\boldsymbol\rho-\boldsymbol\rho';z,z'),
\qquad z,z'\in I.
$$

面内 WS 窗口限制的是相对面内位移。两个密度则各自使用同一 `I` 中独立的高度。
**区间端点是密度表示的边界，不是金属壁。** 自由和极化片模型采用法向无界的
Green 函数；双栅极模型的真实接地边界位于 `z0±d`。

这与 [Ismail-Beigi (2006)](https://doi.org/10.1103/PhysRevB.73.233103)
的 `Theta(zc−|z−z'|)/|r−r'|` 相对高度截断不同：后者截断的是位移的函数，
这里则为每个密度选择连续的法向代表区间，并允许核分别依赖两个高度。

在面内作 Fourier 变换，可进一步看出有限 WS 窗口的作用。若解析无限平面 Green
函数为 `G_infinity(p;z,z')`，则有限算符对应

$$
\widetilde V_W(\mathbf q;z,z')
=\int\frac{d^2p}{(2\pi)^2}\,
\widetilde\chi_{W_\parallel}(\mathbf q-\mathbf p)
G_\infty(\mathbf p;z,z').
$$

实空间乘窗口，在动量空间表现为卷积。直接在单个 q 上代入无限平面的解析公式，
会漏掉这个窗口。本实现积分的是带窗口的算符本身。

#### 4.3 法向 Fourier–Galerkin 投影

WAVECAR 提供的跃迁密度可以表示为 Fourier 系数。法向 Fourier 基底用来表示
区间 `I` 内的密度，而静电 Green 函数仍可采用开放形式。取中心化坐标：

$$
x=z-z_0,\quad g_m=\frac{2\pi m}{L_z},\quad
\phi_m(x)=\frac{e^{ig_mx}}{\sqrt{L_z}},
$$

$$
T_{mn}(\mathbf q)=\frac1{L_z}
\int_{-L_z/2}^{L_z/2}\!dx\int_{-L_z/2}^{L_z/2}\!dx'\,
e^{-ig_mx}\,
\widetilde V_W(\mathbf q;z_0+x,z_0+x')\,e^{ig_nx'}.
$$

$$
K^{\mathrm{global}}_{mn}(\mathbf q)
=e^{i(g_n-g_m)z_0}T_{mn}(\mathbf q).
$$

最后一式把中心化矩阵变为原密度使用的全局坐标 Fourier 约定。这也说明即使在
bare 情形，`WINT_Z0` 也实际参与积分。

自由核虽然只依赖 `z−z'`，在有限开放区间上通常仍得到非对角矩阵，因为有限区间
卷积不是周期 Green 函数所对应的循环卷积。使用 Fourier 基底，并不强迫选定的
静电问题具有周期边界或导体边界。

包含 C 时，面内 Green 函数的量纲为 `eV Å²`，法向矩阵为 `eV Å³`，再由原有
Fourier/体积归一化得到能量；去掉 C 的纯几何矩阵量纲为 `Å²`。

计算中的变量分离来自恒等式

$$
\frac1{\sqrt{\rho^2+(x-x')^2}}
=\frac2{\sqrt\pi}\int_0^\infty dt\,
e^{-t^2\rho^2}e^{-t^2(x-x')^2}.
$$

它把面内 WS 积分与法向矩阵元分离，后者化为一维矩积分。极化片诱导项和栅极
镜像项也有相应表示，因此可以保留轨道厚度，而不必在额外的任意法向网格上直接
采样奇异势。

#### 4.4 参考面及物理假设

`WINT_Z0` 是沿法向的绝对笛卡尔坐标，单位 Å。对 bare/Yukawa，它只确定数值
积分区间的中心；对 RK_SHEET，它还确定理想极化片位置；对 DUAL_GATE，它还
确定两板中面。几何层中心可作为初始选择，但有效极化面及实际器件几何属于物理
模型的选择。

应使材料和有意义的密度尾部完整位于 `I` 中，把端点放在真空内。bare 计算中，
边界在足够真空内移动，应只产生尾部和数值误差范围内的变化。跨越周期坐标原点
的材料层需要先展开成连续的一层；直接平均被周期包裹的分数坐标可能出错。

面内窗口还定义了有限作用范围，不能据此认定所有存储 R 上的相互作用都已达到
无限面积极限。WS 面、法向真空、Fourier 分辨率和物理屏蔽参数需分别判断。
BVK_WS 同样适用于层状体系；当重要的法向电荷间距都在 BvK 法向截断内、尾部又
足够小时，两种 bare 表示可在内部给出一致积分，但完整有限基组矩阵一般不相同。

| 输入/条件 | 含义 |
|---|---|
| `WINT_BOUNDARY=SLAB_WS` | 面内 BvK WS 窗口，加上各密度的开放法向区间 |
| `WINT_Z0` | 有限厚度模型必需的区间中心；相应模型中同时是极化片或两板参考面 |
| `WINT_NORMAL` | 可选法向，默认沿第三晶格矢量 |
| 晶格和 k 网格 | 第三矢量垂直面内；法向一个 k 点；面内可以为斜晶格 |
| `WS_NGRID` | 前两项细化面内核 FFT；本边界不使用第三项 |

RK_SHEET、RK2D 和 DUAL_GATE 当前显式使用此边界；BARE、YUKAWA 可选两种
边界。RK2D 已积分掉高度，不使用 Z0。

<a id="zh-kernels"></a>

### 5. 选择物理相互作用

以下都是投影到同一组 Wannier 轨道上的静态算符。`epsilon` 是相对介电常数，
`epsilon_0` 是真空常数。屏蔽长度属于物理输入，与 FFT 收敛参数不同。

#### 5.1 `BARE`：Wannier 电荷云之间的裸排斥

$$
K_0(\mathbf r,\mathbf r')=\frac{C}{|\mathbf r-\mathbf r'|},
\qquad C=\frac{e^2}{4\pi\epsilon_0}.
$$

它是无屏蔽模型的真空参考。轨道局域程度决定 onsite 的尺度：归一化轨道越集中，
其静电自相互作用通常越强。相距足够远且远离有限域边界时，局域电荷云接近点电荷
的 `C/距离` 相互作用。PBE 轨道改变电荷云的形状，不会自动把介电屏蔽带进裸核。

| 关键词 | 定义 |
|---|---|
| `WINT_KERNEL=BARE` | 显式完整 PAW bare 投影 |
| `WINT_BOUNDARY` | 默认 BVK_WS，也可选 SLAB_WS |
| `WINT_Z0` | 仅 SLAB_WS 必需，用于确定数值积分区间 |

显式 BARE 不设置 `WINT_EPS` 和屏蔽长度。与层状屏蔽模型比较时，使用相同的
SLAB_WS 区域和 Z0，作为有限厚度的共同参考。

#### 5.2 `YUKAWA`：均匀三维介质中的指数屏蔽

均匀介质中的可移动屏蔽电荷，可用修正 Poisson 方程表示：

$$
K_Y(r)=\frac{C}{\epsilon}\frac{e^{-\kappa r}}r,
\qquad
(-\nabla^2+\kappa^2)K_Y
=\frac{4\pi C}{\epsilon}\delta(\mathbf r),
\qquad
\widetilde K_Y(\mathbf Q)=\frac{4\pi C}{\epsilon(Q^2+\kappa^2)}.
$$

额外的逆长度 kappa 使远处电场在 `lambda=1/kappa` 尺度上衰减。函数形式源于
[Yukawa (1935)](https://doi.org/10.11429/ppmsj1919.17.0_48)，这里采用的是上式
明确规定的静态电子屏蔽模型，可用来研究相互作用范围对低能物理的影响。
金属板通过空间边界条件产生屏蔽，其 kernel 有不同的形式。

| 关键词 | 定义 |
|---|---|
| `WINT_KERNEL=YUKAWA` | 指数屏蔽的径向相互作用 |
| `WINT_KAPPA` | 必需，逆屏蔽长度 Å⁻¹，有限非负 |
| `WINT_EPS` | 均匀背景相对介电常数，默认 1，有限正值 |
| `WINT_BOUNDARY` | BVK_WS 或 SLAB_WS；后者还需 Z0 |

kappa=0 时恢复同一积分区域上的 bare/epsilon；epsilon=1 时与显式 BARE 相同。
增大 kappa 会强调短尺度，因此还需检查密度的分辨率。

#### 5.3 `RK_SHEET`：极化局限于平面，轨道电荷保留厚度

二维绝缘体的诱导极化主要集中在材料附近，而电场延伸到周围三维介质，因此屏蔽
随面内波长变化。理想极化片模型描述这种响应，同时保留两团 Wannier 电荷的高度。

把极化片放在 z0，并用下面的 Poisson 问题定义 r0：

$$
\left[-\partial_z^2+q^2+2r_0q^2\delta(z-z_0)\right]
G_s(q;z,z')=\frac{4\pi C}{\epsilon}\delta(z-z').
$$

在两侧取衰减的自由空间解，得到

$$
G_s(q;z,z')=\frac{2\pi C}{\epsilon q}
\left[e^{-q|z-z'|}
-\frac{r_0q}{1+r_0q}
 e^{-q(|z-z_0|+|z'-z_0|)}\right],\qquad q>0.
$$

第一项是背景库仑场，第二项是极化片诱导电荷产生的场。电势连续，法向导数满足

$$
\partial_zG_s(z_0^+,z')-\partial_zG_s(z_0^-,z')
=2r_0q^2G_s(z_0,z'),\qquad z'\ne z_0.
$$

这些条件确定片的响应。这里从理想极化片方程推导有限高度形式；源和观察点均在
片上时，它退化为 [Cudazzo、Tokatly 与 Rubio (2011)](https://doi.org/10.1103/PhysRevB.84.085406)
推导的 Rytova–Keldysh 相互作用，其历史来源为 Rytova 和 Keldysh 的工作
[7–9](#references--参考文献)。

对片上的电荷，`r0 q/(1+r0 q)` 使短波长场受到更强屏蔽；片外电荷的响应还受到
高度相关的指数因子抑制。对位于同一侧、远离极化片的两团电荷，相对片修正在
足够大的 q 下消失。长距离面内相互作用仍保留
`C/(epsilon rho)` 尾部。若电荷在**片外**彼此重合，局域三维库仑奇点仍在；
对数 RK 奇点对应源和观察点均在理想片上的极限。本模型保留电荷厚度，同时把
极化介质近似为单个理想平面，尚未引入空间变化的介电张量。

| 关键词 | 定义 |
|---|---|
| `WINT_KERNEL=RK_SHEET` | 保留两个高度的理想极化片 Green 函数 |
| `WINT_BOUNDARY=SLAB_WS` | 所需的边界表示 |
| `WINT_R0` | 必需，有效屏蔽长度 Å，有限非负 |
| `WINT_Z0` | 必需，极化片位置和积分区间中心，Å |
| `WINT_EPS` | 均匀背景相对介电常数，默认 1，有限正值 |

r0 由上述公式定义。从二维极化率换算时，需统一单位和背景介电约定。r0=0 时
极化片消失，恢复有限高度的 bare/epsilon。

#### 5.4 `RK2D`：严格零厚度近似

当电荷厚度相对于关注的长度尺度足够小时，可先把每个跃迁密度沿高度积分，再
把它放到同一个理想平面上：

$$
G_{2D}(q)=\frac{2\pi C}{\epsilon q(1+r_0q)},\qquad
K_{2D}(\rho)=\frac{\pi C}{2\epsilon r_0}
\left[H_0(\rho/r_0)-Y_0(\rho/r_0)\right],
$$

$$
\rho_{2D}(\boldsymbol\rho)=\int_I dz\,\rho_{3D}(\boldsymbol\rho,z).
$$

实空间特殊函数表达式对应 r0>0；r0=0 使用下文的极限。H0 为 Struve 函数，Y0 为
第二类 Bessel 函数。渐近形式直接体现极化片的作用：

$$
K_{2D}(\rho)\simeq
\begin{cases}
\displaystyle\frac{C}{\epsilon r_0}
\left[\ln\!\left(\frac{2r_0}{\rho}\right)-\gamma\right],&\rho\ll r_0,\\[6pt]
\displaystyle\frac{C}{\epsilon\rho},&\rho\gg r_0.
\end{cases}
$$

短距离的对数行为来自片的极化，长距离库仑尾部则来自延伸到层外的电场。
它与上一模型在片上的屏蔽律相同，作用对象变为二维投影密度，因而可用于判断
忽略 Wannier 厚度的影响。即使厚度对远距离作用不重要，onsite 矩阵元仍探测
较大动量，这时零厚度近似可能有明显影响。把长波长有效片模型延伸到微观 onsite
范围，本身也是需要结合轨道判断的物理近似。

| 关键词 | 定义 |
|---|---|
| `WINT_KERNEL=RK2D` | 作用于法向积分后密度的 RK 核 |
| `WINT_BOUNDARY=SLAB_WS` | 所需的面内有限表示 |
| `WINT_R0` | 必需，屏蔽长度 Å，有限非负，约定同 RK_SHEET |
| `WINT_EPS` | 背景相对介电常数，默认 1，有限正值 |
| `WINT_Z0` | 不设置；法向坐标已被积分掉 |

r0=0 时得到平面的 `C/(epsilon rho)`，而非有限厚度 BARE。这里“二维”指电荷
及极化模型，电场仍在三维空间传播。

#### 5.5 `DUAL_GATE`：两块接地金属板的镜像屏蔽

金属栅极通过诱导表面电荷屏蔽相互作用。两块理想接地板位于 `z=z0−d` 与
`z=z0+d`，令 `x=z−z0`，对应的静电问题为

$$
(-\partial_x^2+q^2)G_g(q;x,x')
=\frac{4\pi C}{\epsilon}\delta(x-x'),
\qquad G_g(q;\pm d,x')=0.
$$

其有限高度 Dirichlet Green 函数为

$$
G_g(q;x,x')=\frac{4\pi C}{\epsilon q}
\frac{\sinh[q(x_<+d)]\sinh[q(d-x_>)]}{\sinh(2qd)},
\quad x_<=\min(x,x'),\quad x_>=\max(x,x').
$$

$$
G_g(q;0,0)=\frac{2\pi C}{\epsilon q}\tanh(qd),
\qquad G_g(0;0,0)=\frac{2\pi Cd}{\epsilon}.
$$

它来自 [Takae–Onuki (2013)](https://doi.org/10.1063/1.4821085) 中的平板 Green
函数及镜像构造；相关平行板电子模型亦见
[Throckmorton–Vafek (2012)](https://doi.org/10.1103/PhysRevB.86.115447)。这里
始终用 d 表示中面到**每一块**板的距离。

在中面放置点电荷时，等价的镜像和横向模展开让空间衰减更直观：

$$
K_g(\rho;0,0)=\frac{C}{\epsilon}
\sum_{n=-\infty}^{\infty}
\frac{(-1)^n}{\sqrt{\rho^2+(2nd)^2}}
=\frac{2C}{\epsilon d}\sum_{j=0}^\infty
K_0\!\left(\frac{(2j+1)\pi\rho}{2d}\right).
$$

最低横向模使远处相互作用正比于 `exp[−pi rho/(2d)]/sqrt(rho)`，衰减尺度为
`2d/pi`。因此即使匹配一个名义屏蔽长度，Yukawa 也不能重现完整的 gate 相互作用。
有限高度形式还能区分更靠近上板或下板的电荷。

| 关键词 | 定义 |
|---|---|
| `WINT_KERNEL=DUAL_GATE` | 保留电荷高度的双接地板模型 |
| `WINT_BOUNDARY=SLAB_WS` | 所需的边界表示 |
| `WINT_GATE_D` | 必需，中面到每侧板的距离 d，Å；两板间距为 2d |
| `WINT_Z0` | 必需，两板中面及积分区间中心，Å |
| `WINT_EPS` | 均匀间隔介质的相对介电常数，默认 1，有限正值 |

当前积分区间需完整位于两板之间，因而要求 `d>Lz/2`。这是积分区域的限制，
不是物理栅极必须放在 DFT 真空晶胞之外的规则。d 足够大时恢复有限高度的
bare/epsilon。单个两电子积分没有额外的 1/2；多体 Hamiltonian 另有成对计数
约定。单电子镜像势则是单粒子跃迁/onsite Hamiltonian 的另一项贡献。

<a id="zh-use"></a>

### 6. PAW 投影与实际使用

#### 选定的算符需一致作用于同一轨道密度

PAW 波函数包含平滑平面波部分和原子中心重建部分。相互作用需要一致作用于二者。
显式 WINT 路径重建固定的完整跃迁密度：

$$
\rho_{ij}=\rho^{\mathrm{PW}}_{ij}
+\sum_a\left(\rho^{\mathrm{AE}}_{ij,a}
-\rho^{\mathrm{PS}}_{ij,a}\right).
$$

再用选定算符收缩该密度，包括平滑–原子及不同原子之间的交叉项。它使用现有
PAW 分波与投影子，不需要另做一遍全电子固体计算。PAW 原理及其 VASP 实现见
[Blöchl (1994)](https://doi.org/10.1103/PhysRevB.50.17953) 和
[Kresse–Joubert (1999)](https://doi.org/10.1103/PhysRevB.59.1758)。

统一的密度构造使不同 Green 函数能作用于相同对象，无需把裸核特有的局域修正
加到另一个屏蔽核上。省略 WINT 标签的旧 bare 路径仍使用其原有补偿和局域修正。
有限截断下二者的差异属于数值问题；新屏蔽模型使用显式路径作为共同参考。

#### 从最终 Wannier 基底到相互作用文件

1. 完成电子结构计算，得到一致的 WAVECAR 与 k 网格。
2. 使用紧凑 SCDM 得到初始 Wannier 投影。
3. 固定最终 Wannier gauge；如需 SAWF，先由外部工具完成并导出兼容 VASP 的 WANPROJ。
4. 运行固定轨道的 `ALGO=2E4WA`，投影所选相互作用。

支持的相互作用路径为 CPU NCL、`ISPIN=1`、`NCORE=1`，包含非共线 SOC spinor。
若实际材料参考面位于笛卡尔 z=0，可在已有输入上增加如下双栅极设置：

```text
ALGO = 2E4WA
HFRCUT = -2
NCORE = 1
WINT_KERNEL = DUAL_GATE
WINT_BOUNDARY = SLAB_WS
WINT_Z0 = 0.0
WINT_GATE_D = 100.0
```

参考面需按实际结构设定；例中的 100 Å 表示每侧板距离 10 nm，只是几何示例。
比较核时保持轨道基底与 PAW 数据一致。VASP 依照原有指标约定写出 VIJKL、VRijkl，
模型参数记录在文件头及 OUTCAR。后续使用的 SAWF 单粒子 Hamiltonian 与这些
相互作用必须对应同一 Wannier 基底。

用于 fRG 时，epsilon、r0 和 kappa 定义环境或部分屏蔽相互作用，应区分输入中
已有的屏蔽与后续 fRG 生成的目标低能带屏蔽。当前标量 epsilon 描述均匀各向同性
介质；各向异性基底或多层器件需要合适的有效模型或更一般的 Poisson Green 函数。
片极化与 gate 在物理上可以共存，目前的 selector 分别描述这些模型，尚未求解
二者共存的 Poisson 问题。

<a id="zh-convergence"></a>

### 7. 哪些收敛问题决定物理解释？

- **轨道模型：** 选带、Wannier gauge 及单粒子能带表示决定具体是哪组电荷云在相互作用。
- **法向真空与参考面：** I 应包含有意义的密度尾部。对 bare，在真空中移动区间边界可
  检查有限表示；对固定材料移动极化片或 gate，则是在改变物理环境。
- **空间范围：** BvK k 网格确定 WS 范围。固定 k 网格增加 WS_NGRID，只会细化积分，
  不会扩大相互作用范围。
- **核分辨率：** 检查 WS_NGRID；WS_TOL 控制核构造及内部积分，不直接给出每个矩阵元的误差条。
- **密度分辨率：** 联合检查 VCUTOFF 与 NGX/NGY/NGZ。短程和 onsite 作用涉及较大的
  密度动量；只增大 VCUTOFF 可能触及原 FFT 网格上限。

当扩大密度 FFT 网格后旧 CHGCAR 不再兼容时，`ICHARG=12` 已在固定 WAVECAR 的
相互作用投影中检查过；对应的单电子/FCIDUMP 输出不能替代已有 SAWF Hamiltonian。
物理屏蔽参数与上述数值收敛方向分别选择。

独立 Green 函数检验、法向矩阵和原生投影/输出比较验证了实现的有限基组算符。
具体材料的相互作用精度仍由以上收敛问题确定。若核逐点为正、密度也为非负的精确
密度，密度直接项应为正；一般跃迁密度矩阵元没有这一符号要求。

<a id="zh-install"></a>

### 安装与升级

当前源码配置支持 VASP 6.6.1 和 6.6.0；原生本地运行使用 6.6.1 CPU NCL。使用者需
具有合法 VASP 许可及自己的源码包，并使用已有成功编译环境中的编译器、MPI 和
Wannier90 配置。对干净原版或可识别的仅 SCDM 源码执行：

```bash
sha256sum -c manifest.sha256
VASP_WINT_SRC=/absolute/path/to/vasp-source
bash apply_patch.sh --status "$VASP_WINT_SRC"
bash apply_patch.sh --check "$VASP_WINT_SRC"
bash apply_patch.sh --apply "$VASP_WINT_SRC"
make -C "$VASP_WINT_SRC" ncl DEPS=1 -j4
```

升级时，**先用当前已安装版本配套的脚本和补丁**恢复原版源码：

```bash
bash /path/to/current-package/apply_patch.sh --reverse "$VASP_WINT_SRC"
```

然后替换为新版安装包，再检查、应用并编译。回退必须发生在覆盖旧文件或在原目录执行 `git pull` 之前，
因为 `--reverse` 会核对本版本载荷的源码哈希。仓库无需额外保留多组历史回退包。
可识别的干净原版或仅 SCDM 源码可直接使用新版安装脚本。脚本会保留源码备份。

部署只需相邻的 `apply_patch.sh` 和统一 `.patch`。`SOURCE_IDENTITIES.json`
记录源码身份，`manifest.sha256` 用于核对包完整性。补丁历史文件名包含上文全部
功能；本仓库提供第三方补丁，不提供 VASP 可执行文件、PAW 数据集或完整源码包。

## References / 参考文献

1. A. Damle, L. Lin and L. Ying, “Compressed Representation of Kohn–Sham Orbitals via Selected Columns of the Density Matrix,” *J. Chem. Theory Comput.* **11**, 1463–1469 (2015). [DOI](https://doi.org/10.1021/ct500985f).
2. A. Damle, L. Lin and L. Ying, “SCDM-k: Localized orbitals for solids via selected columns of the density matrix,” *J. Comput. Phys.* **334**, 1–15 (2017). [DOI](https://doi.org/10.1016/j.jcp.2016.12.053), [arXiv](https://arxiv.org/abs/1507.03354).
3. A. Damle and L. Lin, “Disentanglement via Entanglement: A Unified Method for Wannier Localization,” *Multiscale Model. Simul.* **16**, 1392–1410 (2018). [DOI](https://doi.org/10.1137/17M1129696), [arXiv](https://arxiv.org/abs/1703.06958). See also the [VASP SCDM interface](https://vasp.at/wiki/LSCDM).
4. R. Sundararaman and T. A. Arias, “Ideal regularization of the Coulomb singularity in exact exchange by Wigner–Seitz truncated interactions,” *Phys. Rev. B* **87**, 165122 (2013). [DOI](https://doi.org/10.1103/PhysRevB.87.165122), [arXiv](https://arxiv.org/abs/1302.6204).
5. S. Ismail-Beigi, “Truncation of periodic image interactions for confined systems,” *Phys. Rev. B* **73**, 233103 (2006). [DOI](https://doi.org/10.1103/PhysRevB.73.233103), [arXiv](https://arxiv.org/abs/cond-mat/0603448).
6. H. Yukawa, “On the Interaction of Elementary Particles. I,” *Proc. Phys.-Math. Soc. Japan* **17**, 48–57 (1935). [DOI](https://doi.org/10.11429/ppmsj1919.17.0_48).
7. N. S. Rytova, “The screened potential of a point charge in a thin film” (1967). [Publisher record](https://vmu.phys.msu.ru/en/abstract/1967/3/1967-3-030/), [archival English version](https://arxiv.org/abs/1806.00976).
8. L. V. Keldysh, “Coulomb interaction in thin semiconductor and semimetal films,” *JETP Lett.* **29**, 658–661 (1979). The original model is discussed in Ref. 9.
9. P. Cudazzo, I. V. Tokatly and A. Rubio, “Dielectric screening in two-dimensional insulators: Implications for excitonic and impurity states in graphane,” *Phys. Rev. B* **84**, 085406 (2011). [DOI](https://doi.org/10.1103/PhysRevB.84.085406), [arXiv](https://arxiv.org/abs/1104.3346).
10. K. Takae and A. Onuki, “Applying electric field to charged and polar particles between metallic plates: Extension of the Ewald method,” *J. Chem. Phys.* **139**, 124108 (2013), Eqs. 2.17–2.21. [DOI](https://doi.org/10.1063/1.4821085), [arXiv](https://arxiv.org/abs/1305.3503).
11. R. E. Throckmorton and O. Vafek, “Fermions on bilayer graphene: symmetry breaking for B=0 and nu=0,” *Phys. Rev. B* **86**, 115447 (2012). [DOI](https://doi.org/10.1103/PhysRevB.86.115447), [arXiv](https://arxiv.org/abs/1111.2076).
12. P. E. Blöchl, “Projector augmented-wave method,” *Phys. Rev. B* **50**, 17953 (1994). [DOI](https://doi.org/10.1103/PhysRevB.50.17953).
13. G. Kresse and D. Joubert, “From ultrasoft pseudopotentials to the projector augmented-wave method,” *Phys. Rev. B* **59**, 1758 (1999). [DOI](https://doi.org/10.1103/PhysRevB.59.1758).

Native treatment and output definitions: [VASP HFRCUT](https://vasp.at/wiki/HFRCUT), [Coulomb singularity](https://vasp.at/wiki/Coulomb_singularity), [VRijkl](https://vasp.at/wiki/VRijkl).
