# Wannier models and electrostatic interactions in VASP

[English](#english) · [中文说明](#中文说明) · [References / 参考文献](#references--参考文献)

| Feature / 功能 | English | 中文 |
|---|---|---|
| Band-compact SCDM / 紧凑选带 | [Motivation and use](#en-scdm) | [动机与使用](#zh-scdm) |
| Wigner–Seitz Coulomb truncation / WS 库仑截断 | [Principle](#en-ws) · [Representations and use](#en-ws-choice) | [原理](#zh-ws) · [表示与使用](#zh-ws-choice) |
| Physical interaction models / 物理相互作用模型 | [Electrostatic environments](#en-kernels) | [静电环境](#zh-kernels) |
| PAW projection and Wannier workflow / PAW 投影与 Wannier 工作流 | [Projection and use](#en-use) | [投影与使用](#zh-use) |
| Physical and numerical convergence / 物理与数值收敛 | [Convergence](#en-convergence) | [收敛](#zh-convergence) |
| Apply / reverse / update | [Installation](#en-install) | [安装升级](#zh-install) |

## English

### 1. The physical problem: from Bloch states to an interacting Wannier model

A low-energy Wannier model describes electronic motion through hopping
amplitudes and electron interactions through matrix elements between orbital
charge distributions. In a moiré material, only a few bands may govern the
physics near the Fermi level, while the underlying plane-wave calculation
contains thousands. Constructing the model requires both a localized basis
for those bands and an interaction suited to the material's environment.

For a scalar, spin-independent kernel, the interaction is projected through
the bilinear form below. The density-direct element gives the electrostatic
repulsion between two normalized Wannier charge clouds:

$$
\mathcal B_K[\rho,\sigma]
=\int d^3r\,d^3r'\,\rho^*(\mathbf r)
K(\mathbf r,\mathbf r')\sigma(\mathbf r'),\qquad
D_{mn}(\mathbf R)=\mathcal B_K[|w_{m0}|^2,|w_{n\mathbf R}|^2].
$$

The extent and internal structure of an orbital set its onsite repulsion;
off-center matrix elements describe how the interaction decays between
orbitals. Their magnitudes relative to the bandwidth help identify the
interaction regime. General transition-density matrix elements also describe
exchange and orbital-changing processes. VASP's index and file conventions
are given in the [VRijkl documentation](https://vasp.at/wiki/VRijkl).

The patch provides three tools for constructing these models:

1. **Band-compact SCDM:** make the initial Wannier projection affordable in a
   large-band calculation.
2. **Wigner–Seitz Coulomb truncation:** represent the interaction on a domain
   matched to the finite k mesh, using a full-cell or layer-resolved operator.
3. **Physical interaction models:** project vacuum Coulomb, homogeneous
   screening, polarizable-sheet and metallic-gate interactions in that basis.

The Wannier basis fixes the orbital charge distributions. The kernel specifies
their physical interaction, and the finite integration domain determines its
representation on the k mesh. Keeping WAVECAR and the final WANPROJ fixed
allows different interactions to be compared in the same basis, including
bases obtained from PBE or PBE/SOC calculations.

<a id="en-scdm"></a>

### 2. Band-compact SCDM: constructing the initial Wannier basis with less memory

#### Motivation and principle

SCDM constructs an initial localized basis from selected real-space columns
of a density matrix or a spectrally filtered projector. QR decomposition with
column pivoting identifies the anchor points that determine these columns.
The original SCDM method, its periodic extension SCDM-k and the treatment of
entangled bands are described in Refs. [1–3](#references--参考文献).

A spectral filter weights the states used to construct the Wannier subspace:

$$
P_f(\mathbf k)=\sum_{n\in\mathcal A}
f(E_{n\mathbf k})\,|u_{n\mathbf k}\rangle\langle u_{n\mathbf k}|,
\qquad
M_{\mathrm{SCDM}}\propto
N_{\mathrm{grid}}N_{\mathrm{spinor}}N_{\mathrm{active}}.
$$

For example, the ERFC filter smoothly reduces the weight of states above a
chosen energy:

$$
f(E)=\frac12\mathrm{erfc}\!\left(\frac{E-\mu}{\sigma}\right).
$$

Here `CUTOFF_MU` sets the filter center and `CUTOFF_SIGMA` sets its energy
width. Together they control which part of the spectrum contributes to the
initial subspace.

In native VASP, SCDM work arrays can retain the full band count even when
Wannier90 excludes most of those bands. The cost becomes substantial in moiré
cells, where both the real-space grid and the number of bands are large.
This patch applies the existing `exclude_bands` selection to SCDM's wavefunction,
eigenvalue, filter and transformation arrays. Retained states remain ordered
by their original band indices, preserving the correspondence between AMN
rows and VASP bands.

The dominant SCDM buffers then scale with `N_active` in place of `NBANDS`.
The global wavefunction/WAVECAR arrays and the native MMN temporary array keep
their original dimensions. The memory saving therefore comes from the SCDM
working space within the larger VASP calculation.

#### How to use it

Select the bands in the usual Wannier90 input and use the VASP/Wannier90
interface. For an isolated group containing as many retained bands as Wannier
functions, unit weights define the subspace projector:

```text
LSCDM = .TRUE.
LWANNIER90 = .TRUE.
NUM_WANN = 8
CUTOFF_TYPE = none
```

In this example, the accompanying `exclude_bands` selection retains eight
suitable bands. For entangled bands, retain a larger energy window and
specify the filter explicitly.

| Input | Meaning in the compact subspace |
|---|---|
| Wannier90 `exclude_bands` | Original VASP bands omitted from the SCDM working space |
| `NUM_WANN` | Number of Wannier functions to construct |
| `CUTOFF_TYPE=none` | Unit weights on the retained states |
| `CUTOFF_MU`, `CUTOFF_SIGMA` | Explicit spectral-filter parameters; ERFC uses the center/width convention above |

When bands are excluded and the active space contains more bands than
Wannier functions, ERFC/Gaussian require both filter parameters, Fermi requires
`CUTOFF_MU`, and `num_wann` requires `CUTOFF_SIGMA`. Each filter retains the
parameter convention of this VASP version; the center/width formula above
defines ERFC. With exactly `NUM_WANN` active bands, omitting the filter also
gives unit weights.

The resulting SCDM projection is the starting point for Wannier localization.
External tools can impose symmetry through SAWF and export the final WANPROJ,
which the interaction calculation reads as the fixed orbital transformation.

<a id="en-ws"></a>

### 3. Coulomb truncation according to the Wigner–Seitz cell

Projecting the Coulomb interaction onto Wannier orbitals samples a singular
reciprocal-space kernel on a finite k mesh. WS truncation connects the
treatment of that singularity to the shape of the real-space supercell
represented by the mesh.

#### 3.1 Why a finite k mesh needs Coulomb regularization

In a plane-wave calculation, interaction matrix elements involve the Coulomb
factor `4πC/Q²`, where `Q=k′−k+G` is the momentum transferred between Bloch
states and `C=e²/(4πepsilon_0)` is the vacuum Coulomb prefactor. At `Q=0` this
factor diverges. Its continuum integral is finite: the three-dimensional
volume element contributes a factor proportional to Q², which cancels the
inverse square near the origin. A discrete k-point sum instead assigns a
finite weight to the singular point. It therefore needs a prescription for
the finite contribution of the surrounding reciprocal-space region.

The same issue has a useful real-space interpretation. Reconstructing Wannier
orbitals and pair densities from a uniform Gamma-centered mesh repeats them on the
Born–von Karman (BvK) superlattice. For a full mesh without further reductions,

$$
\mathbf A_i=N_i\mathbf a_i,\qquad
N_k=N_1N_2N_3,\qquad V_B=N_k\Omega.
$$

Here `a_i` are the simulation-cell vectors, `N_i` are the mesh dimensions and
Omega is the simulation-cell volume. The relevant electrostatic length scales
are consequently those of this k-mesh supercell. Choosing the finite-mesh
Coulomb operator determines how the projected charges interact with their
artificial supercell copies.

#### 3.2 Native treatments in VASP

VASP offers two approaches described in its
[Coulomb singularity](https://vasp.at/wiki/Coulomb_singularity) and
[HFRCUT](https://vasp.at/wiki/HFRCUT) documentation.

**Periodic kernel with a singularity correction.** With the native default
`HFRCUT=0`, a regular mesh and the unscreened Coulomb law use the probe-charge
Ewald correction. Other supported meshes or kernels use the auxiliary-function
treatment. These prescriptions supply the missing singular contribution
through a model-charge correction. In the regular unscreened case, the
underlying periodic Green function and its correction can be written as

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

The uniform term supplies the neutralizing background. The Madelung
correction removes the regular image and background contribution at the
origin. At finite separation, the potential still contains their
distance-dependent field. This distinction matters when extracting extended
Wannier interactions over a substantial fraction of the supercell.

**Spherical truncation.** Native VASP also replaces the Coulomb law by a
finite-radius kernel. A positive `HFRCUT` gives the radius explicitly;
`HFRCUT=-1` chooses a sphere with the BvK volume:

$$
K_{\mathrm{sph}}(r)=\frac{C}{r}\Theta(R_c-r),\qquad
R_c=\left(\frac{3N_k\Omega}{4\pi}\right)^{1/3}.
$$

This radius assumes the complete uniform mesh. The truncated kernel has a
finite Fourier head and changes the nonzero Fourier components as well.

#### 3.3 Limitations in two-dimensional and moiré cells

An isolated material layer is represented by a periodic cell with vacuum
along its normal. Orbital overlap across the vacuum can become negligible
while Coulomb image interactions remain appreciable, because the electric
field decays much more slowly than the orbital density. Moiré cells add long
in-plane lattice scales and extended Wannier charge clouds to this problem.

With a mesh containing one normal k point, the normal BvK period remains
`Lz`, while the in-plane periods grow as **`N_i a_i`**. The resulting supercell
can be much wider than it is tall. An origin-based correction to the periodic
kernel then leaves a separation-dependent error over the in-plane distances
needed by an extended Hubbard model. In the tSnS calculations motivating this
patch, `HFRCUT=0` produced a negative tail in the density-direct terms. Analysis
of the operator traced this tail to the periodic Green function and its
zero-mode correction. The target here is the interaction within one isolated
layer, so the finite-mesh operator must recover that spatial dependence.

A sphere introduces a further geometric limitation. For an approximately
orthogonal `41 × 39 × 29 Å` cell and a `6 × 6 × 1` mesh, the BvK lengths are
about `246 × 233 × 29 Å`. The automatic sphere has radius about `73.5 Å`,
whereas the normal WS half-width is about `14.5 Å` and the in-plane half-widths
are about `117–123 Å`. The sphere spans repeated normal cells while cutting
off part of the useful in-plane range. Sampling its Fourier transform
periodizes that sphere, retaining contributions where the copies overlap.
Increasing the in-plane mesh at fixed vacuum makes this shape mismatch more
pronounced. A cutoff adapted to the actual repetition lattice resolves the
geometric problem.

#### 3.4 The Wigner–Seitz construction

[Sundararaman and Arias (2013)](https://doi.org/10.1103/PhysRevB.87.165122)
analyze finite-k exchange in terms of localized Wannier pair densities. Their
central observation is that the slowly converging image contribution can be
removed by retaining the Coulomb interaction within the WS cell of the k-mesh
superlattice. The cutoff then follows the nearest-image geometry in every
lattice direction:

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

Inside this domain, the interaction is the vacuum Coulomb law. Its Fourier
transform supplies the entire finite-mesh operator, including a finite head.
For orthogonal vectors the domain is a box with half-widths `|A_i|/2`;
oblique cells use the actual lattice metric. The WS cells tile space, so the
periodized truncated interaction follows `C/r` throughout the interior of
the central cell.

The paper's construction covers bulk and partially periodic systems; for a
layer, the unsampled normal direction has one k point. Its exchange-energy
convergence follows the localization of the occupation-weighted Wannier-like
functions: exponential localization gives exponential convergence with
supercell linear extent, while zero-temperature metals retain algebraic
behavior. In the present application, individual Wannier matrix elements are
converged through their density tails, available interaction range and
numerical resolution.

#### 3.5 Full-cell representation: BVK_WS

The full-cell representation applies this construction to the complete BvK
WS domain. A translation-invariant interaction is a function of the
three-dimensional relative displacement, so its plane-wave representation
is a scalar multiplier at each reciprocal momentum. The current `BVK_WS`
path evaluates bare Coulomb and Yukawa interactions in this form.

The bare numerical kernel follows the range-separated construction in
Appendix A.1 of Sundararaman and Arias. The short-range erfc part is integrated
analytically; the smooth erf/r part is sampled using nearest-image distances
and transformed on the WS Fourier grid. The numerical separation and grid
are controlled by `WS_TOL` and `WS_NGRID`. The physical prefactor and the
native density-contraction normalization are applied by the VASP interface.

For a layer with a perpendicular third vector and one normal k point, the
full WS domain is an in-plane WS cell times a normal interval of half-width
`Lz/2` in the relative displacement. This is a useful bare-interaction
representation when the important normal charge separations lie within that
range. It also supports general three-dimensional lattice metrics. Its use
for layered materials follows directly from the published WS construction.

<a id="en-slab"></a>

#### 3.6 Layer-resolved representation: SLAB_WS

##### Purpose: retain the charge heights and the surrounding environment

A polarizable sheet or a metal plate singles out a height in the device.
Moving two charges together toward a plate changes the induced surface
charge and their interaction, even though their mutual displacement stays
fixed. The Green function must therefore retain the source and observation
heights separately, in the form `V(rho−rho'; z,z')`.

The scalar full-cell path represents an operator depending on a relative
three-dimensional coordinate. To handle the height-dependent environment,
the layered representation keeps the in-plane WS window and integrates one
continuous copy of each orbital density along the normal. This also supplies
an open-normal bare or Yukawa reference with the same finite orbital thickness
as the sheet and gate calculations. The resulting operator, defined below,
combines planar WS truncation with a height-resolved Green function.

The finite-height sheet and gate models use this matrix representation.
The strict RK2D approximation shares its planar WS integration, after first
integrating out the orbital height. Section 4 gives the physical Green
function for each environment.

##### Definition: a planar WS window and a normal density interval

Choose the infinite-plane interaction `V_infinity(s;z,z')` appropriate to free
space, a polarizable sheet or grounded plates. Let `W_parallel` be the WS cell
of the two in-plane BvK vectors. The finite-domain operator used here is

$$
I=[z_0-L_z/2,z_0+L_z/2],\qquad
K_W(\boldsymbol\rho,z;\boldsymbol\rho',z')
=\chi_{W_\parallel}(\boldsymbol\rho-\boldsymbol\rho')
V_\infty(\boldsymbol\rho-\boldsymbol\rho';z,z'),
\qquad z,z'\in I.
$$

The window limits the relative in-plane displacement. Each charge retains
its own height within `I`, an interval containing a single continuous layer
and its density tails. The free and sheet Green functions extend to infinity
along the normal; the gate Green function is grounded at `z0±d`. Those Green
functions provide the electrostatic boundaries, while the endpoints of `I`
specify where the saved orbital densities are integrated.

This layered operator is the extension implemented in this patch. Its
ingredients are WS windowing and layered electrostatic Green functions.
For comparison, the slab cutoff
`Theta(zc−|z−z'|)/|r−r'|` of
[Ismail-Beigi (2006)](https://doi.org/10.1103/PhysRevB.73.233103) limits the
relative normal displacement. The present construction integrates the two
heights independently, allowing the potential to respond to their positions
relative to a sheet or gate.

Writing the infinite-plane Green function as `G_infinity(p;z,z')`, the planar
Fourier transform is

$$
\widetilde V_W(\mathbf q;z,z')
=\int\frac{d^2p}{(2\pi)^2}\,
\widetilde\chi_{W_\parallel}(\mathbf q-\mathbf p)
G_\infty(\mathbf p;z,z').
$$

A real-space window becomes a momentum-space convolution. Each finite-WS
coefficient therefore combines the infinite-plane Green function over the
momenta selected by the window. The same windowed operator supplies the
head and the nonzero planar components.

##### Projection onto the normal Fourier basis

The WAVECAR pair densities are available as Fourier coefficients. Their
height dependence is expanded on `I` and contracted with the open-normal
Green function. In centered coordinates,

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

The final phase converts the centered matrix into the global Fourier
convention of the supplied density. It makes `WINT_Z0` part of the finite
interval projection even for bare Coulomb. Integration over two independently
bounded heights generally gives a dense matrix in the normal indices,
including for the free kernel that depends only on `z−z'`. The planar momentum
remains diagonal because translation symmetry is preserved in the plane.

Including `C`, the planar Green function has units `eV Å²` and the normal
matrix has units `eV Å³`. VASP's Fourier and volume normalization converts
the density contraction into an energy. Factoring out `C` leaves a geometric
matrix with units `Å²`.

The implementation separates the spatial integrations through

$$
\frac1{\sqrt{\rho^2+(x-x')^2}}
=\frac2{\sqrt\pi}\int_0^\infty dt\,
e^{-t^2\rho^2}e^{-t^2(x-x')^2}.
$$

The planar Gaussian factor is integrated over the WS window. The normal
matrix elements reduce to one-dimensional moments, with corresponding
representations for the induced sheet and gate-image contributions. This
retains the orbital thickness and treats the singularity analytically in
the normal integration.

##### Positioning the layer, sheet and gates

`WINT_Z0` is an absolute Cartesian coordinate along the declared normal, in Å.
For bare and Yukawa it centers the numerical density interval. For RK_SHEET
it also places the polarizable sheet, and for DUAL_GATE it places the midpoint
between the plates. The layer's geometric center gives an initial estimate;
the effective sheet and gate positions are physical model choices.

Choose `I` with both endpoints in vacuum and the significant charge tails
inside it. A layer crossing the periodic coordinate origin should first be
unwrapped into one continuous copy. For bare interactions, shifting the
interval boundaries within sufficient vacuum checks tail and numerical
errors. Moving a physical sheet or gate relative to the fixed layer changes
the electrostatic environment.

For thin, localized layer densities whose significant relative heights are
less than `Lz/2`, the full-cell normal cutoff leaves those charge pairs
unchanged. Bare integrals in the interior can then agree with the layered
evaluation. Finite tails, normal resolution and domain boundaries determine
the residual difference. The layered representation is chosen for its
open-normal integral and height-dependent operator; the full-cell form
remains an appropriate reference for translation-invariant interactions.

<a id="en-ws-choice"></a>

#### 3.7 Choosing the representation and setting the input

Choose the representation from the physical operator and charge geometry:

| Property | Full-cell representation: BVK_WS | Layer-resolved representation: SLAB_WS |
|---|---|---|
| Finite domain | Full three-dimensional BvK WS cell in the relative coordinate | In-plane BvK WS window; each density occupies the normal interval I |
| Normal electrostatics | Translation-invariant relative-coordinate kernel | Open or gate-bounded Green function retaining both heights |
| Fourier operator | Scalar at each three-dimensional momentum | Matrix in normal Fourier indices at each planar momentum |
| Implemented models | Bare Coulomb, Yukawa | Bare Coulomb, Yukawa, finite-height RK, RK2D, dual gate |
| Lattice geometry | General three-dimensional metric | Perpendicular third vector and one normal k point; oblique in-plane cells supported |
| Reference plane | No independent height reference | Density-interval center; also the sheet position or gate midpoint when applicable |
| Typical use | Bulk interactions or a localized layer with a suitable normal cutoff | Layer-resolved environmental screening and a matching finite-height reference |

For the present implementation, the sheet, RK2D and gate selectors use
`SLAB_WS`; BARE and YUKAWA support either representation. RK2D integrates
out height and consequently has no Z0 dependence. Converge the planar range,
normal density containment and Fourier resolution for the chosen model.

The shared VASP entry is `HFRCUT=-2`, used with `ALGO=2E4WA`. An explicit
`WINT_KERNEL` selects the physical interaction and `WINT_BOUNDARY` selects
its representation. Both full-cell and layered calculations use this entry.
With all `WINT_*` tags omitted, it retains the original bare full-cell
evaluation; explicitly selecting BARE uses the reconstructed-PAW density
route explained in Section 5. Native VASP interprets an unpatched `-2` through
its automatic-sphere branch.

| Input | Definition |
|---|---|
| `HFRCUT=-2` | Enable this patch's WS representation for the Wannier-integral workflow |
| `WINT_BOUNDARY=BVK_WS` | Full BvK WS domain; default boundary for explicit BARE/YUKAWA |
| `WS_NGRID` | Three-component kernel FFT grid; automatic by default |
| `WS_TOL` | Internal kernel construction/quadrature control; default `1e-8` |

| Input / condition | Meaning |
|---|---|
| `WINT_BOUNDARY=SLAB_WS` | In-plane BvK WS window plus open-normal density interval |
| `WINT_Z0` | Required interval center for finite-height models; also the sheet/gate reference plane when applicable |
| `WINT_NORMAL` | Optional Cartesian normal, default along the third lattice vector |
| Lattice and k mesh | Third vector perpendicular to the plane; one normal k point; oblique in-plane cells are supported |
| `WS_NGRID` | First two components refine the plane kernel FFT; third component is unused by this boundary |

Both representations require a complete, regular, automatic Gamma-centered
k mesh with no Fock mesh reductions. The layered form additionally uses one
normal k point and a third lattice vector perpendicular to the plane.
Oblique in-plane vectors are supported. A supplied `WINT_NORMAL` must be
parallel to the third vector; its default is along that vector.

For a finalized Wannier basis and otherwise unchanged interaction input,
an explicit full-cell bare calculation uses

```text
HFRCUT = -2
WINT_KERNEL = BARE
WINT_BOUNDARY = BVK_WS
```

For a layer centered at Cartesian z=0 with the default normal, the matching
layer-resolved bare reference uses

```text
HFRCUT = -2
WINT_KERNEL = BARE
WINT_BOUNDARY = SLAB_WS
WINT_Z0 = 0.0
```

Set the reference coordinate from the actual structure. Screened calculations
then select the desired kernel and its parameters from Section 4 while
keeping the finalized Wannier basis fixed. The current dual-gate integration
requires both plates to enclose the complete interval, giving `d>Lz/2`.


<a id="en-kernels"></a>

### 4. Physical interaction models

The following static kernels describe different electrostatic environments
for a fixed Wannier basis. Their formulas define the physical interaction
before the finite-domain projection in Section 3. The dielectric parameter
`epsilon` is the relative permittivity, while `epsilon_0` is the vacuum
permittivity. Screening lengths set physical scales in the model; FFT
parameters control numerical resolution.

#### 4.1 Coulomb repulsion in vacuum

$$
K_0(\mathbf r,\mathbf r')=\frac{C}{|\mathbf r-\mathbf r'|},
\qquad C=\frac{e^2}{4\pi\epsilon_0}.
$$

The bare kernel provides the vacuum reference. A more concentrated
normalized orbital generally has a larger electrostatic self-interaction,
which sets the onsite scale. Well-separated localized densities approach
the point-charge `C/distance` law away from finite-domain boundaries. The
electronic-structure calculation determines the Wannier charge clouds;
dielectric screening enters through the choice of interaction kernel.

| Keyword | Meaning |
|---|---|
| `WINT_KERNEL=BARE` | Explicit full-PAW bare projection |
| `WINT_BOUNDARY` | BVK_WS by default, or SLAB_WS |
| `WINT_Z0` | Required only with SLAB_WS; sets its numerical interval |

For explicit BARE, omit `WINT_EPS` and the screening lengths. Use the same
SLAB_WS domain and Z0 when it serves as the finite-height reference for the
layered screened models.

#### 4.2 Homogeneous screening: the Yukawa interaction

In a homogeneous medium, a local linear screening response leads to the
modified Poisson equation and its Yukawa solution:

$$
K_Y(r)=\frac{C}{\epsilon}\frac{e^{-\kappa r}}r,
\qquad
(-\nabla^2+\kappa^2)K_Y
=\frac{4\pi C}{\epsilon}\delta(\mathbf r),
\qquad
\widetilde K_Y(\mathbf Q)=\frac{4\pi C}{\epsilon(Q^2+\kappa^2)}.
$$

The inverse length `kappa` sets the exponential decay length
`lambda=1/kappa`. The functional form originates with
[Yukawa (1935)](https://doi.org/10.11429/ppmsj1919.17.0_48); here it describes
the static electronic-screening model specified by the displayed equation.
Varying this length provides a simple way to study the effect of interaction
range. Metallic gates produce a spatially structured response through their
boundary conditions, as described below.

| Keyword | Meaning |
|---|---|
| `WINT_KERNEL=YUKAWA` | Exponentially screened radial interaction |
| `WINT_KAPPA` | Required inverse screening length, Å⁻¹; finite and nonnegative |
| `WINT_EPS` | Homogeneous relative permittivity, default 1; finite and positive |
| `WINT_BOUNDARY` | BVK_WS or SLAB_WS; the latter also requires Z0 |

At `kappa=0` the kernel becomes bare/epsilon on the same domain, recovering
explicit BARE when `epsilon=1`. As kappa increases, shorter length scales
contribute more strongly relative to the distant interaction and require
adequate density resolution.

#### 4.3 Screening by a polarizable sheet with finite-height charges

In a two-dimensional insulator, the induced polarization is concentrated
in the material while the electric field extends into the surrounding
three-dimensional dielectric. Screening therefore depends on the in-plane
wavelength. An ideal polarizable sheet describes this response while the
Wannier charge clouds retain their finite height distributions.

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

determine the induced sheet response. Solving this ideal-sheet problem gives
the finite-height expression above. When both charges lie on the sheet, it
reduces to the Rytova–Keldysh interaction derived by
[Cudazzo, Tokatly and Rubio (2011)](https://doi.org/10.1103/PhysRevB.84.085406)
and based on the earlier work of Rytova and Keldysh [7–9](#references--参考文献).

On the sheet, the factor `r0 q/(1+r0 q)` increasingly suppresses
short-wavelength fields. At finite height, the exponential weights the
induced response by each charge's distance from the sheet. For two charges
on the same side and away from the sheet, the relative correction vanishes
at sufficiently large q; the distant in-plane interaction retains its
`C/(epsilon rho)` tail. Coincident charges away from the sheet retain the
local three-dimensional Coulomb singularity, while charges on the sheet have
the logarithmic RK singularity. The approximation treats the polarizable
material as one ideal plane within a uniform background dielectric.

| Keyword | Meaning |
|---|---|
| `WINT_KERNEL=RK_SHEET` | Finite-height ideal polarizable-sheet Green function |
| `WINT_BOUNDARY=SLAB_WS` | Required boundary representation |
| `WINT_R0` | Required effective screening length, Å; finite and nonnegative |
| `WINT_Z0` | Required sheet position and integration-interval center, Å |
| `WINT_EPS` | Uniform background relative permittivity, default 1; finite and positive |

The displayed equation defines `r0`, including its background-dielectric
convention. Use that convention and consistent units when converting a
two-dimensional polarizability into r0. At `r0=0` the sheet response vanishes,
leaving the finite-height bare/epsilon interaction.

#### 4.4 The strict two-dimensional Rytova–Keldysh interaction

When the charge thickness is small compared with the relevant length
scales, integrate each pair density over height and place the resulting
planar distribution on the sheet:

$$
G_{2D}(q)=\frac{2\pi C}{\epsilon q(1+r_0q)},\qquad
K_{2D}(\rho)=\frac{\pi C}{2\epsilon r_0}
\left[H_0(\rho/r_0)-Y_0(\rho/r_0)\right],
$$

$$
\rho_{2D}(\boldsymbol\rho)=\int_I dz\,\rho_{3D}(\boldsymbol\rho,z).
$$

For r0>0, the real-space expression uses the Struve function `H0` and the
Bessel function of the second kind `Y0`, with the asymptotic behavior

$$
K_{2D}(\rho)\simeq
\begin{cases}
\displaystyle\frac{C}{\epsilon r_0}
\left[\ln\!\left(\frac{2r_0}{\rho}\right)-\gamma\right],&\rho\ll r_0,\\[6pt]
\displaystyle\frac{C}{\epsilon\rho},&\rho\gg r_0.
\end{cases}
$$

Sheet polarization produces the short-distance logarithm, while electric
field outside the layer gives the distant Coulomb tail. Applying this
on-sheet screening law to planar densities provides a direct comparison
with the finite-height model. The thickness approximation is most consequential
for onsite terms, which probe large momenta. Their interpretation also
depends on how well the long-wavelength sheet model describes the microscopic
scale of the chosen orbitals.

| Keyword | Meaning |
|---|---|
| `WINT_KERNEL=RK2D` | RK operator acting on height-integrated pair densities |
| `WINT_BOUNDARY=SLAB_WS` | Required in-plane finite-domain representation |
| `WINT_R0` | Required screening length, Å; finite and nonnegative; same convention as RK_SHEET |
| `WINT_EPS` | Background relative permittivity, default 1; finite and positive |
| `WINT_Z0` | Omit: the normal coordinate has been integrated out |

At `r0=0`, the limiting kernel is the planar `C/(epsilon rho)` interaction.
Recovering finite-height BARE additionally requires retaining the normal
density coordinates. The two-dimensional label refers to the charges and
polarization; their electric field propagates in three dimensions.

#### 4.5 Screening by two grounded metal plates

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

The plate Green function and its image construction are given by
[Takae and Onuki (2013)](https://doi.org/10.1063/1.4821085). A parallel-plate
interaction also enters the electronic model of
[Throckmorton and Vafek (2012)](https://doi.org/10.1103/PhysRevB.86.115447).
Here `d` is the distance from the midpoint to each plate, so their separation
is twice that distance.

For point charges in the midplane, the image series can also be expanded
in transverse modes:

$$
K_g(\rho;0,0)=\frac{C}{\epsilon}
\sum_{n=-\infty}^{\infty}
\frac{(-1)^n}{\sqrt{\rho^2+(2nd)^2}}
=\frac{2C}{\epsilon d}\sum_{j=0}^\infty
K_0\!\left(\frac{(2j+1)\pi\rho}{2d}\right).
$$

The lowest transverse mode controls the distant tail,
`exp[−pi rho/(2d)]/sqrt(rho)`, with decay scale `2d/pi`. The full gate kernel
contains all transverse modes and the two charge heights. It therefore
captures both the spatial decay and the changing response as charges
approach either plate.

| Keyword | Meaning |
|---|---|
| `WINT_KERNEL=DUAL_GATE` | Two grounded plates, retaining both charge heights |
| `WINT_BOUNDARY=SLAB_WS` | Required boundary representation |
| `WINT_GATE_D` | Required midpoint-to-each-plate distance d in Å; plate spacing is 2d |
| `WINT_Z0` | Required plate midpoint and integration-interval center, Å |
| `WINT_EPS` | Homogeneous spacer relative permittivity, default 1; finite and positive |

The current implementation places the entire density interval between the
plates, requiring `d>Lz/2`. This condition follows from the integration
domain used here. At large d, the operator approaches finite-height
bare/epsilon. The Green function enters each two-electron integral with the
normalization shown above; the many-body Hamiltonian supplies the factor
1/2 for pair counting. The one-electron image potential contributes separately
to the hopping and onsite Hamiltonian.

<a id="en-use"></a>

### 5. PAW projection and the fixed-Wannier workflow

#### The same orbital density must meet the chosen operator

PAW represents a wavefunction through a smooth plane-wave part and
atom-centered reconstruction terms. Applying a chosen interaction consistently
requires the complete pair density. The explicit WINT route reconstructs it as

$$
\rho_{ij}=\rho^{\mathrm{PW}}_{ij}
+\sum_a\left(\rho^{\mathrm{AE}}_{ij,a}
-\rho^{\mathrm{PS}}_{ij,a}\right).
$$

Contracting the operator with this density includes the cross terms between
the smooth and atom-centered contributions, as well as those between atoms.
The reconstruction uses the PAW partial waves and projectors supplied with
the calculation. The PAW formalism and its VASP implementation are described
in [Blöchl (1994)](https://doi.org/10.1103/PhysRevB.50.17953) and
[Kresse and Joubert (1999)](https://doi.org/10.1103/PhysRevB.59.1758).

A common reconstructed density lets every screened Green function act on
the same orbital charge distribution. The original bare route, selected by
omitting WINT tags, uses its native compensation and local-correction
construction. Comparing the two at finite cutoff checks their numerical
agreement; the explicit route supplies the common bare reference for the
screened models.

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

Set the reference plane from the actual structure. In this example, each
gate lies 100 Å, or 10 nm, from that plane. Keep the Wannier basis and PAW
datasets fixed when comparing kernels. VASP writes VIJKL and VRijkl using
its native index conventions and records model parameters in the file
headers and OUTCAR. Use these interactions with the one-electron SAWF
Hamiltonian of the same Wannier basis.

For fRG, epsilon, r0 and kappa define the environmental or partially screened
interaction supplied to the low-energy model. This input screening should
exclude the target-band processes generated by the subsequent fRG flow.
The scalar epsilon assumes a homogeneous isotropic medium; an anisotropic
substrate or multilayer device calls for a suitable effective model or a
more general Poisson Green function. The current sheet and gate selectors
solve their respective electrostatic problems separately. A device containing
both requires the corresponding combined Green function.

<a id="en-convergence"></a>

### 6. Convergence and physical interpretation

- **Orbital model:** band selection, Wannier gauge and the one-electron band
  representation determine which charge clouds interact.
- **Normal vacuum and placement:** contain the significant density tails in I.
  For bare, moving the boundaries within vacuum checks the finite-domain
  representation; moving a sheet or gate relative to a fixed layer changes
  the physical environment.
- **Spatial range:** enlarge the BvK k mesh to increase the WS interaction range.
- **Kernel resolution:** refine WS_NGRID at fixed range. WS_TOL sets the internal
  construction and quadrature tolerance; convergence of the resulting matrix
  elements determines their numerical accuracy.
- **Density resolution:** converge VCUTOFF together with NGX/NGY/NGZ. Short-range
  and onsite interactions probe high-momentum density components. Raising
  VCUTOFF alone can encounter a density-grid ceiling.

When a larger density grid is incompatible with the existing CHGCAR,
`ICHARG=12` can be used for the fixed-WAVECAR interaction projection, as checked
in this workflow. Continue to use the established SAWF Hamiltonian for the
one-electron model. Screening lengths and dielectric constants specify the
physical model throughout these numerical convergence checks.

The implemented finite-basis operator has been checked against independent
Green-function calculations, normal-matrix tests and native projection/file
output comparisons. The convergence checks above establish the accuracy for
a particular material. A positive pointwise kernel acting on nonnegative
exact orbital densities gives positive density-direct terms. General
transition-density matrix elements need not be positive.

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

Then install the new package, check and apply it, and rebuild. Complete the
reverse step before overwriting the old script/patch pair or running `git pull`
in that directory, because `--reverse` checks the source hashes for its own
version. Recognized clean vendor or SCDM-only sources can use the current installer
directly. The script retains source backups.

Keep `apply_patch.sh` beside the unified `.patch` file when deploying.
`SOURCE_IDENTITIES.json` records the supported source hashes, and
`manifest.sha256` checks package integrity. Despite its historical filename,
the patch contains all features described above. This repository distributes
the third-party patch; obtain VASP sources, binaries and PAW datasets through
the licensed VASP distribution.

## 中文说明

### 1. 物理目标：从 Bloch 能带到相互作用 Wannier 模型

低能 Wannier 模型用跃迁振幅描述电子运动，用轨道电荷分布之间的矩阵元描述
电子相互作用。对 moiré 材料，费米能附近的少数能带可能已包含主要低能物理，
而底层平面波计算仍需数千条能带。构造模型时，需要先为目标能带建立局域轨道
基底，再选择与材料环境相符的相互作用。

对于与自旋无关的标量核，相互作用投影可写成下面的双线性形式。其中，密度
直接项给出两个归一化 Wannier 电荷云之间的静电排斥：

$$
\mathcal B_K[\rho,\sigma]
=\int d^3r\,d^3r'\,\rho^*(\mathbf r)
K(\mathbf r,\mathbf r')\sigma(\mathbf r'),\qquad
D_{mn}(\mathbf R)=\mathcal B_K[|w_{m0}|^2,|w_{n\mathbf R}|^2].
$$

轨道的大小与内部结构决定 onsite 排斥，off-center 矩阵元则描述轨道间相互作用
的空间衰减。将这些能量尺度与带宽比较，可以判断模型所处的相互作用区间。
一般的跃迁密度矩阵元还描述交换和轨道变化过程。VASP 的指标及文件约定见
[VRijkl 说明](https://vasp.at/wiki/VRijkl)。

本补丁为模型构造提供三项功能：

1. **Band-compact SCDM：** 降低大能带计算中初始 Wannier 投影的内存需求。
2. **Wigner–Seitz 库仑截断：** 采用与有限 k 网格匹配的积分区域，提供完整超胞和
   保留层状高度信息的两种表示。
3. **物理相互作用模型：** 在同一轨道基底中投影真空库仑作用、均匀屏蔽、极化片和
   金属栅极所对应的相互作用。

Wannier 基底确定轨道电荷分布，物理核定义电荷之间的相互作用，有限积分区域
规定它在 k 网格上的表示。固定 WAVECAR 和最终 WANPROJ 后，便可在同一组轨道
中比较不同相互作用；这组轨道可以来自 PBE 或 PBE/SOC 计算。

<a id="zh-scdm"></a>

### 2. Band-compact SCDM：降低初始 Wannier 基底的内存需求

#### 动机和原理

SCDM 从密度矩阵或谱过滤投影算符中选取实空间列，构造局域轨道的初始基底。
带列主元的 QR 分解确定锚点，锚点对应的列随后用于建立投影。原始 SCDM 方法、
周期推广 SCDM-k 和纠缠能带的处理见参考文献 [1–3](#references--参考文献)。

谱过滤为参与子空间构造的态分配权重：

$$
P_f(\mathbf k)=\sum_{n\in\mathcal A}
f(E_{n\mathbf k})\,|u_{n\mathbf k}\rangle\langle u_{n\mathbf k}|,
\qquad
M_{\mathrm{SCDM}}\propto
N_{\mathrm{grid}}N_{\mathrm{spinor}}N_{\mathrm{active}}.
$$

例如，ERFC 过滤在选定能量附近平滑降低高能态的权重：

$$
f(E)=\frac12\mathrm{erfc}\!\left(\frac{E-\mu}{\sigma}\right).
$$

这里，`CUTOFF_MU` 设定过滤中心，`CUTOFF_SIGMA` 设定能量宽度，两者共同
控制哪些能态参与初始子空间的构造。

原生 VASP 的 SCDM 工作矩阵可能仍按完整能带数分配，即使 Wannier90 已排除
其中大部分能带。moiré 晶胞同时增大实空间网格和能带数，使这部分内存开销尤为
突出。本补丁将已有的 `exclude_bands` 选带结果用于 SCDM 的波函数、能量、过滤
权重及变换矩阵。保留态按原始能带编号升序排列，从而保持 AMN 行与 VASP 能带
之间的对应关系。

主要 SCDM 工作矩阵的能带维度由 `NBANDS` 缩小为 `N_active`。全局波函数、
WAVECAR 和原有 MMN 临时矩阵保持原来的维度，因此内存节省来自 VASP 整体计算
中的 SCDM 工作空间。

#### 怎样使用

通过常规 VASP/Wannier90 接口，在 Wannier90 输入中指定选带。若孤立能带组中
的保留带数恰好等于 Wannier 函数数目，单位权重即可定义子空间投影：

```text
LSCDM = .TRUE.
LWANNIER90 = .TRUE.
NUM_WANN = 8
CUTOFF_TYPE = none
```

本例中，配套的 `exclude_bands` 应保留八条合适的能带。对于纠缠能带，需保留
更大的能量窗口，并显式指定过滤方式。

| 输入 | 在紧凑子空间中的作用 |
|---|---|
| Wannier90 `exclude_bands` | 从 SCDM 工作空间排除的原始 VASP 能带 |
| `NUM_WANN` | 需要构造的 Wannier 函数数目 |
| `CUTOFF_TYPE=none` | 保留态均使用单位权重 |
| `CUTOFF_MU`、`CUTOFF_SIGMA` | 显式谱过滤参数；ERFC 的中心与宽度定义见上式 |

排除部分能带后，若保留带数仍多于 Wannier 函数数目，ERFC/Gaussian 需同时给出
两个过滤参数，Fermi 需给出 `CUTOFF_MU`，`num_wann` 需给出 `CUTOFF_SIGMA`。
各过滤形式沿用本版本 VASP 的参数约定，上面的中心与宽度公式定义 ERFC。
若保留带数恰好为 `NUM_WANN`，省略过滤设置也会得到单位权重。

SCDM 投影为后续 Wannier 局域化提供初始基底。外部工具可通过 SAWF 施加对称性
约束并导出最终 WANPROJ；相互作用计算读取这一固定的轨道变换。

<a id="zh-ws"></a>

### 3. 基于 Wigner–Seitz 单元的库仑截断

平面波方法将库仑积分写成动量空间求和。有限 k 网格下，零动量附近的贡献与
实空间重复镜像密切相关。Wigner–Seitz（WS）截断从采样超晶格的几何出发，
为库仑相互作用建立与该网格匹配的有限区域表示。

#### 3.1 有限 k 网格为何需要处理库仑奇点

裸库仑核的倒空间形式为 `4πC/Q²`，其中 `C=e²/(4πepsilon_0)`。动量转移包含
两个 Bloch 波矢之差及一个倒格矢。在 `Q=0` 处，核逐点发散，但三维体积元含有平方动量因子，因此零点
附近的积分有限。连续积分包含这一小体积的贡献；离散求和则需要用有限数值
表示它。直接删去零模会遗漏这部分贡献，使长程库仑积分对 k 网格收敛缓慢。

这一问题也有实空间解释。对于无额外约化的完整 Γ 中心均匀网格，Bloch 态的离散
Fourier 变换给出在 Born–von Karman（BvK）超晶格上重复的 Wannier 轨道及
跃迁密度。原晶格矢量、采样点数与该超晶格的关系为

$$
\mathbf A_i=N_i\mathbf a_i,\qquad
N_k=N_1N_2N_3,\qquad V_B=N_k\Omega.
$$

因此，k 网格同时决定实空间的重复周期和可分辨的相互作用范围。原胞体积为
Ω，BvK 体积为原胞的 N_k 倍。处理库仑奇点时，对有限动量求和所作的选择也会
反映在这些重复电荷分布之间的相互作用中。

#### 3.2 VASP 中的原生处理

原生 VASP 的默认设置 `HFRCUT=0` 在使用规则网格计算未屏蔽相互作用时采用
probe-charge Ewald 修正，其他支持的情况使用辅助函数处理。两者用模型电荷
估计零动量附近的贡献，将其计入
有限网格求和。具体定义见 [HFRCUT](https://vasp.at/wiki/HFRCUT) 和
[Coulomb singularity](https://vasp.at/wiki/Coulomb_singularity) 文档。

在本工作流使用的规则、未屏蔽情形中，周期 Green 函数满足

$$
-\nabla^2G_B(\mathbf r)
=4\pi\left[\sum_{\mathbf S\in\Lambda_B}
\delta(\mathbf r-\mathbf S)-\frac1{V_B}\right],
\qquad \langle G_B\rangle=0.
$$

右端包含周期重复的点电荷和均匀补偿背景，电势的晶胞平均值取零。相应的
局域 Madelung 修正及有限网格核为

$$
K_{\mathrm{per}}(\mathbf r)=C\,[G_B(\mathbf r)-v_M],
\qquad
v_M=\lim_{r\to0}\left[G_B(\mathbf r)-\frac1r\right].
$$

该修正抵消原点处镜像与背景的正则贡献，使短距离电势恢复孤立库仑势的参考。
在有限距离上，周期镜像和背景仍产生随间距变化的贡献。BvK 描述 k 采样对应的
周期性，Ewald 求和则用于计算这一周期体系上的静电作用。

VASP 还提供球形截断，将相互作用限制在给定半径内。正值 `HFRCUT` 直接指定
半径，`HFRCUT=-1` 自动选取与 BvK 超胞等体积的球：

$$
K_{\mathrm{sph}}(r)=\frac{C}{r}\Theta(R_c-r),\qquad
R_c=\left(\frac{3N_k\Omega}{4\pi}\right)^{1/3}.
$$

球形截断使零动量分量有限，并通过一个长度参数控制空间作用范围。上述自动
半径使用完整均匀网格的体积；其适用程度取决于球与实际采样超胞的形状是否
相近。

#### 3.3 二维材料与 moiré 超胞中的局限

为孤立二维材料构造扩展 Hubbard 模型，需要同一材料层内 Wannier 电荷云之间
的相互作用。平面波计算沿法向加入真空，真空厚度通常按抑制层间轨道重叠来
选择。库仑场的作用距离远长于轨道重叠范围，因此周期重复层仍会影响长距离
矩阵元，尤其是在法向周期短于面内 moiré 尺度的晶胞中。

本项目的 tSnS 计算提供了一个例子。采用原生 `HFRCUT=0` 时，密度直接项在
远距离出现负尾部。对同一 Wannier 基底中投影算符的分析表明，这一尾部来自
周期 Green 函数及零模修正。此处需要的是孤立材料层中的排斥相互作用，因而
要重新选择有限 k 网格下的库仑表示。

球形截断另有形状上的限制。相关长度是 **`N_i a_i`**，包含原晶胞与 k 网格
共同产生的各向异性。法向仅一个 k 点时，BvK 法向周期保持为 `Lz`，面内范围
却随网格增加而扩大。

例如，约 `41 × 39 × 29 Å` 的近似正交晶胞配合 `6 × 6 × 1` 网格，得到约
`246 × 233 × 29 Å` 的 BvK 尺寸。自动球形半径约为 `73.5 Å`，而法向 WS
半宽约为 `14.5 Å`，面内半宽约为 `117–123 Å`。球在法向跨越重复晶胞，
同时截掉一部分有用的面内长距离。对球形核的离散 Fourier 采样使其周期重复，
这种重叠可能保留与重复材料层的相互作用。随采样超晶格变化的截断形状，更适合
这一高度各向异性的几何。

#### 3.4 Wigner–Seitz 截断的原理

[Sundararaman–Arias (2013)](https://doi.org/10.1103/PhysRevB.87.165122) 从
Wannier 跃迁密度的实空间局域性出发，分析了有限 k 网格的库仑误差。其核心
做法是将核限制在采样超晶格的 WS 单元中，在该区域内部保留原来的库仑作用，
并抑制采样超胞之间的人工耦合。WS 单元由原点相对于其他超晶格点的最近区域
组成，因此截断形状自然遵循晶格度量。

完整 BvK 单元上的定义为

$$
W_B=\{\mathbf s:\ |\mathbf s|\leq|\mathbf s-\mathbf S|
\text{ for every }\mathbf S\in\Lambda_B\},
\qquad
K_{\mathrm{WS}}(\mathbf s)=\frac{C}{|\mathbf s|}\chi_{W_B}(\mathbf s).
$$

其 Fourier 分量由同一实空间核统一给出：

$$
\widetilde K_{\mathrm{WS}}(\mathbf Q)
=\int_{W_B}d^3s\,\frac{C\,e^{-i\mathbf Q\cdot\mathbf s}}{|\mathbf s|}.
$$

积分区域有限，零模也随之有限。截断改变所有动量分量，并以这一统一构造替代
原生周期核及其零模修正。区域内的物理定律仍是 `C/r`，有限范围由晶格和
k 网格决定；环境屏蔽则通过后文的物理核定义。

文中的 WS 超胞方法涵盖体材料及部分周期体系。占据数加权的 Wannier 型函数
具有指数局域性时，交换能误差随采样超胞的线性尺寸指数衰减；零温金属则保持
代数收敛。用于 Wannier 相互作用时，单个矩阵元的精度取决于参与积分的密度
尾部、可用空间范围和数值分辨率。

#### 3.5 完整 BvK 单元中的表示（BVK_WS）

`BVK_WS` 在完整三维 BvK WS 单元上表示相互作用。正交晶格对应各方向半宽为
`|A_i|/2` 的长方体；斜晶格按实际晶格度量确定最近镜像。裸库仑和 Yukawa 核
均只依赖相对位移，截断后的周期表示因此可在完整倒空间中用标量 Fourier
乘子作用于密度。

这一表示也适用于层状材料。第三晶格矢量垂直面内，且法向只有一个 k 点时，
法向 WS 边界位于原胞的半宽处，面内范围则随 k 网格扩大。若材料电荷集中在一层中，参与所关注积分
的主要法向相对距离小于 `Lz/2`，且密度尾部足够小，完整 BvK 构造可以保留
这些电荷之间的裸相互作用。接近 WS 面的长距离项则需要更充分的范围与尾部
收敛检查。

对裸库仑核，数值上采用短程与平滑部分的分离：用 erfc(alpha r)/r 的解析 Fourier 项处理
原点附近的短程贡献，在 WS 网格上对平滑的 erf(alpha r)/r 部分作 FFT。
分离参数 alpha 由 WS 内切半径和 `WS_TOL` 确定；两部分相加恢复所需的截断
核。`WS_NGRID` 控制平滑部分的 Fourier 采样及 WS 面附近的积分误差。这一
分离参数用于数值构造，物理屏蔽长度由所选模型另行给出。

<a id="zh-slab"></a>

#### 3.6 保留轨道厚度的层状表示（SLAB_WS）

##### 层状环境为何需要两个高度坐标

极化片和金属板在器件中选定了一个高度。将两团电荷一起向某块金属板移动，
即使保持它们的相对位移不变，诱导表面电荷产生的屏蔽也会改变。因此，这类
相互作用具有 `V(rho−rho'; z,z')` 的形式，需要分别保留源和观察点的高度。
Wannier 电荷分布本身也有有限厚度，这些高度应参与积分。

`SLAB_WS` 是这里为层状环境定义的表示：面内采用 BvK WS 窗口，法向在一个
连续区间内保留两个密度坐标，并使用所选物理环境的 Green 函数。它将 WS 有限
范围与开放法向静电问题结合起来，具体算符由下文定义。

##### 物理 Green 函数与有限区域算符

先选择无限平面相互作用 `V_infinity(s;z,z')`，例如自由空间、极化片或接地
金属板的 Green 函数。以 `W_parallel` 表示两根面内 BvK 矢量生成的 WS 区域，
每个密度在以 z0 为中心、长度为 Lz 的区间内表示：

$$
I=[z_0-L_z/2,z_0+L_z/2],\qquad
K_W(\boldsymbol\rho,z;\boldsymbol\rho',z')
=\chi_{W_\parallel}(\boldsymbol\rho-\boldsymbol\rho')
V_\infty(\boldsymbol\rho-\boldsymbol\rho';z,z'),
\qquad z,z'\in I.
$$

WS 窗口限制相对面内位移，两个密度则在同一 `I` 内各自保留高度坐标。区间
端点规定密度的表示范围，物理静电边界由 Green 函数给定：自由空间和极化片
模型在法向无界，双栅极模型在 `z0±d` 处接地。

[Ismail-Beigi (2006)](https://doi.org/10.1103/PhysRevB.73.233103) 的 slab 截断
`Theta(zc−|z−z'|)/|r−r'|` 限制相对法向位移。这里的构造则为两个密度各自
选择连续法向区间，并保留环境 Green 函数对独立高度的依赖。已有库仑截断和
层状静电学提供其物理基础，上式规定本实现采用的组合。

面内 WS 窗口还决定解析 Green 函数怎样进入有限网格。将无限平面 Green
函数记为 `G_infinity(p;z,z')`，带窗口算符的 Fourier 变换为

$$
\widetilde V_W(\mathbf q;z,z')
=\int\frac{d^2p}{(2\pi)^2}\,
\widetilde\chi_{W_\parallel}(\mathbf q-\mathbf p)
G_\infty(\mathbf p;z,z').
$$

实空间窗口产生动量空间卷积。每个有限 WS 分量都包含窗口选定的各动量贡献，
因此实现中积分完整的带窗口算符。

##### 法向 Fourier–Galerkin 投影

从 WAVECAR 得到的跃迁密度可表示为 Fourier 系数。为将其与上述 Green 函数
收缩，在 `I` 上用 Fourier 基底表示高度依赖。取中心化坐标：

$$
x=z-z_0,\quad g_m=\frac{2\pi m}{L_z},\quad
\phi_m(x)=\frac{e^{ig_mx}}{\sqrt{L_z}},
$$

法向矩阵元同时对源和观察点积分：

$$
T_{mn}(\mathbf q)=\frac1{L_z}
\int_{-L_z/2}^{L_z/2}\!dx\int_{-L_z/2}^{L_z/2}\!dx'\,
e^{-ig_mx}\,
\widetilde V_W(\mathbf q;z_0+x,z_0+x')\,e^{ig_nx'}.
$$

转换到原密度的全局 Fourier 约定后，矩阵带有参考面相位：

$$
K^{\mathrm{global}}_{mn}(\mathbf q)
=e^{i(g_n-g_m)z_0}T_{mn}(\mathbf q).
$$

即使自由核只依赖 `z−z'`，在有限区间上一般也会得到稠密矩阵。这里计算的是
开放 Green 函数在两个高度坐标上的积分；周期 Green 函数所对应的循环卷积才
由法向 Fourier 基底对角化。静电边界由 Green 函数决定，Fourier 基底负责表示
密度及算符。

计入 C 后，面内 Green 函数的量纲为 `eV Å²`，法向矩阵为 `eV Å³`。与密度
收缩后，VASP 的 Fourier 和体积归一化将其转换为能量。提取 C 后，几何矩阵的
量纲为 `Å²`。

空间积分可通过以下恒等式分离：

$$
\frac1{\sqrt{\rho^2+(x-x')^2}}
=\frac2{\sqrt\pi}\int_0^\infty dt\,
e^{-t^2\rho^2}e^{-t^2(x-x')^2}.
$$

面内因子在 WS 窗口内积分，法向矩阵元化为一维矩积分。极化片诱导项和栅极
镜像项也有相应表示。通过这些矩积分，计算保留轨道厚度，并以积分表示处理
库仑奇点。

##### 密度区间与参考面的选择

法向密度区间应包含材料和显著密度尾部，两端放在真空内。若材料层跨越周期
坐标原点，先将原子坐标展开为连续的一层，再确定其参考面。在足够真空内
移动 bare 模型的区间端点，可检查剩余的尾部误差和数值误差。

`WINT_Z0` 是沿法向的绝对笛卡尔坐标，单位为 Å。对 bare/Yukawa，它确定密度
区间中心；对 RK_SHEET，它同时确定极化片位置；对 DUAL_GATE，它同时确定
两板中面。因此，相对固定材料移动极化片或栅极会改变物理模型。几何层中心
可作为初始估计，有效极化面与器件位置应依据实际环境选取。全局 Fourier
相位使这个坐标在 bare 情形下也参与有限区间投影。

面内窗口确定有限相互作用范围。趋近无限面积极限时，需要扩大该范围，并
检查法向真空和 Fourier 分辨率。物理屏蔽参数定义在这一过程中保持固定的
环境。

<a id="zh-ws-choice"></a>

#### 3.7 如何选择表示及设置输入

对只依赖相对位移的裸库仑或 Yukawa 相互作用，可用完整 BvK WS 单元。对于
较薄的单层，只要主要法向间距落在 BvK 截断内，且密度尾部足够小，BVK_WS 与
开放法向表示可以在区域内部给出一致积分。两种有限区域的完整矩阵通常不同。

极化片与栅极需要保留相互作用对各自高度的依赖。当前实现用 SLAB_WS 表示
RK_SHEET 和 DUAL_GATE；它也为 bare/Yukawa 提供保留厚度的开放法向参考。
RK2D 先将跃迁密度沿高度积分，使用同一面内 WS 基础，因此同样选择 SLAB_WS，
法向积分消除了对 Z0 的依赖。

| 比较项 | BVK_WS | SLAB_WS |
|---|---|---|
| 有限区域 | 完整三维 BvK WS 单元 | 面内 BvK WS 窗口与每个密度的法向区间 |
| 法向静电作用 | 只依赖相对坐标的核 | 保留两个高度的开放或栅极边界 Green 函数 |
| Fourier 算符 | 完整倒空间中的标量乘子 | 面内动量分量与法向矩阵 |
| 物理核 | BARE、YUKAWA | BARE、YUKAWA、RK_SHEET、RK2D、DUAL_GATE |
| 几何 | 一般晶格度量 | 第三矢量垂直面内，法向一个 k 点 |
| 参考面 | 无需独立参考面 | 区间中心，并在相应模型中确定极化片或栅极位置 |
| 主要用途 | 体材料或法向截断范围合适的局域单层 | 层状环境屏蔽及与其配套的有限厚度参考 |

补丁通过 `HFRCUT=-2` 启用 WS 截断。所有显式 WINT 模型，包括采用 SLAB_WS
的模型，都需要这一设置；它用于固定 Wannier 基底的 `ALGO=2E4WA` 相互作用
计算。`WINT_BOUNDARY` 选择有限区域表示，`WINT_KERNEL` 选择物理核。

未打补丁的 VASP 6.6.1 源码将所有负值，包括 `-2`，送入自动球形分支。
本补丁为 `-2` 赋予 WS 含义。

| 输入 | 定义 |
|---|---|
| `HFRCUT=-2` | 在 Wannier 积分工作流中启用本补丁的 WS 表示 |
| `WINT_BOUNDARY=BVK_WS` | 完整 BvK WS 区域；显式 BARE/YUKAWA 的默认边界 |
| `WS_NGRID` | 三分量核 FFT 网格；默认自动选择 |
| `WS_TOL` | 核构造及内部积分控制，默认 `1e-8` |

默认采用自动的核 FFT 网格。以下设置使用显式重建 PAW 密度，给出完整 BvK
单元中的裸库仑参考；将这些行加入已有的固定轨道相互作用输入即可：

```text
HFRCUT = -2
WINT_KERNEL = BARE
WINT_BOUNDARY = BVK_WS
```

同一裸核的开放法向参考改用 SLAB_WS，并指定实际材料的区间中心。以下示例
假定参考面位于沿法向的笛卡尔坐标 0 Å：

```text
HFRCUT = -2
WINT_KERNEL = BARE
WINT_BOUNDARY = SLAB_WS
WINT_Z0 = 0.0
```

| 输入/条件 | 含义 |
|---|---|
| `WINT_BOUNDARY=SLAB_WS` | 面内 BvK WS 窗口，加上各密度的开放法向区间 |
| `WINT_Z0` | 有限厚度模型必需的区间中心；相应模型中同时是极化片或两板参考面 |
| `WINT_NORMAL` | 可选法向，默认沿第三晶格矢量 |
| 晶格和 k 网格 | 第三矢量垂直面内；法向一个 k 点；面内可以为斜晶格 |
| `WS_NGRID` | 前两项细化面内核 FFT；本边界不使用第三项 |

`WINT_NORMAL` 必须平行于第三晶格矢量，默认取该矢量方向。面内晶格可以倾斜。
WS 构造使用完整、自动生成的 Γ 中心均匀网格，当前路径保持 Fock 网格完整，
要求关闭网格约化、平移及奇偶子网格选项。对 SLAB_WS，前两项 `WS_NGRID`
控制面内核分辨率，法向密度精度通过密度 Fourier 截断与网格收敛。

省略全部 `WINT_*` 标签时，`HFRCUT=-2` 沿用原有 bare-WS 的 PAW 补偿及局域
修正构造；显式 `WINT_KERNEL=BARE` 使用完整重建密度。两者计算同一裸库仑律，
以不同有限基组密度构造求值。比较不同物理核时，显式路径提供共同的密度参考，
具体 PAW 投影见后文。

屏蔽计算按第 4 节选择物理核及其参数，并保持最终 Wannier 基底固定。当前双栅极
积分要求两块板包围整个密度区间，因此还需满足 `d>Lz/2`。

<a id="zh-kernels"></a>

### 4. 物理相互作用模型

以下静态核描述同一 Wannier 基底所处的不同静电环境。各公式先定义物理相互作用，
再通过第 3 节的有限区域表示进行投影。`epsilon` 是相对介电常数，`epsilon_0` 是
真空介电常数。屏蔽长度确定模型中的物理尺度，FFT 参数控制计算分辨率。

#### 4.1 真空中的库仑排斥

$$
K_0(\mathbf r,\mathbf r')=\frac{C}{|\mathbf r-\mathbf r'|},
\qquad C=\frac{e^2}{4\pi\epsilon_0}.
$$

裸核提供真空中的相互作用参考。归一化轨道越集中，其静电自相互作用通常越大，
从而提高 onsite 排斥。远离有限区域边界时，相隔足够远的局域电荷云趋于点电荷
的 `C/距离` 相互作用。电子结构计算确定 Wannier 电荷云的形状，介电屏蔽则由
相互作用核引入。

| 关键词 | 定义 |
|---|---|
| `WINT_KERNEL=BARE` | 显式完整 PAW bare 投影 |
| `WINT_BOUNDARY` | 默认 BVK_WS，也可选 SLAB_WS |
| `WINT_Z0` | 仅 SLAB_WS 必需，用于确定数值积分区间 |

显式 BARE 应省略 `WINT_EPS` 和屏蔽长度。作为层状屏蔽模型的有限厚度参考时，
使用与之相同的 SLAB_WS 区域和 Z0。

#### 4.2 均匀介质中的 Yukawa 屏蔽

在均匀介质中，局域线性屏蔽响应给出修正 Poisson 方程及其 Yukawa 解：

$$
K_Y(r)=\frac{C}{\epsilon}\frac{e^{-\kappa r}}r,
\qquad
(-\nabla^2+\kappa^2)K_Y
=\frac{4\pi C}{\epsilon}\delta(\mathbf r),
\qquad
\widetilde K_Y(\mathbf Q)=\frac{4\pi C}{\epsilon(Q^2+\kappa^2)}.
$$

逆长度 kappa 决定指数衰减长度 `lambda=1/kappa`。这一函数形式源于
[Yukawa (1935)](https://doi.org/10.11429/ppmsj1919.17.0_48)，这里用它描述上述
方程定义的静态电子屏蔽。改变这一长度，可以研究相互作用范围对低能物理的
影响。下文的金属栅极则通过空间边界条件产生位置相关的屏蔽响应。

| 关键词 | 定义 |
|---|---|
| `WINT_KERNEL=YUKAWA` | 指数屏蔽的径向相互作用 |
| `WINT_KAPPA` | 必需，逆屏蔽长度 Å⁻¹，有限非负 |
| `WINT_EPS` | 均匀背景相对介电常数，默认 1，有限正值 |
| `WINT_BOUNDARY` | BVK_WS 或 SLAB_WS；后者还需 Z0 |

kappa=0 时，核退化为同一区域上的 bare/epsilon；再取 epsilon=1，即恢复显式
BARE。随着 kappa 增大，短尺度相对远程作用的贡献增加，需要相应的密度分辨率。

#### 4.3 极化片屏蔽与有限厚度轨道

二维绝缘体的诱导极化集中在材料附近，电场则延伸到周围三维介质，因此屏蔽
随面内波长变化。理想极化片描述这一响应，Wannier 电荷云同时保留有限的高度
分布。

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

这些条件确定极化片的诱导响应，解出上面的有限高度形式。源与观察点均在片上
时，它退化为 [Cudazzo、Tokatly 与 Rubio (2011)](https://doi.org/10.1103/PhysRevB.84.085406)
推导的 Rytova–Keldysh 相互作用，其基础是 Rytova 与 Keldysh 的早期工作
[7–9](#references--参考文献)。

在极化片上，`r0 q/(1+r0 q)` 对短波长场的抑制更强。电荷位于片外时，指数
因子还按它们到极化片的距离调节诱导响应。对于同侧且远离片的两团电荷，相对
修正在足够大的 q 下消失；远距离面内作用保留 `C/(epsilon rho)` 尾部。
片外重合电荷具有局域三维库仑奇点，片上电荷则具有对数 RK 奇点。这一近似把
极化材料视为均匀背景介质中的单个理想平面。

| 关键词 | 定义 |
|---|---|
| `WINT_KERNEL=RK_SHEET` | 保留两个高度的理想极化片 Green 函数 |
| `WINT_BOUNDARY=SLAB_WS` | 所需的边界表示 |
| `WINT_R0` | 必需，有效屏蔽长度 Å，有限非负 |
| `WINT_Z0` | 必需，极化片位置和积分区间中心，Å |
| `WINT_EPS` | 均匀背景相对介电常数，默认 1，有限正值 |

上述方程定义 r0 及其背景介电约定。从二维极化率换算时，需采用一致的单位和
同一介电约定。r0=0 时，极化片响应消失，剩下有限高度的 bare/epsilon 相互作用。

#### 4.4 严格二维的 Rytova–Keldysh 相互作用

当电荷厚度远小于所关注的长度尺度时，先将每个跃迁密度沿高度积分，再把得到
的二维分布放在极化片上：

$$
G_{2D}(q)=\frac{2\pi C}{\epsilon q(1+r_0q)},\qquad
K_{2D}(\rho)=\frac{\pi C}{2\epsilon r_0}
\left[H_0(\rho/r_0)-Y_0(\rho/r_0)\right],
$$

$$
\rho_{2D}(\boldsymbol\rho)=\int_I dz\,\rho_{3D}(\boldsymbol\rho,z).
$$

r0>0 时，实空间表达式包含 Struve 函数 H0 和第二类 Bessel 函数 Y0，
其渐近形式为

$$
K_{2D}(\rho)\simeq
\begin{cases}
\displaystyle\frac{C}{\epsilon r_0}
\left[\ln\!\left(\frac{2r_0}{\rho}\right)-\gamma\right],&\rho\ll r_0,\\[6pt]
\displaystyle\frac{C}{\epsilon\rho},&\rho\gg r_0.
\end{cases}
$$

片极化产生短距离的对数行为，层外电场则给出远距离库仑尾部。把片上的屏蔽律
作用于二维密度，可直接与有限高度模型比较。厚度近似对涉及大动量的 onsite
矩阵元影响更大；解释这些矩阵元时，还需考察长波长极化片模型在轨道微观尺度
上的适用程度。

| 关键词 | 定义 |
|---|---|
| `WINT_KERNEL=RK2D` | 作用于法向积分后密度的 RK 核 |
| `WINT_BOUNDARY=SLAB_WS` | 所需的面内有限表示 |
| `WINT_R0` | 必需，屏蔽长度 Å，有限非负，约定同 RK_SHEET |
| `WINT_EPS` | 背景相对介电常数，默认 1，有限正值 |
| `WINT_Z0` | 不设置；法向坐标已被积分掉 |

r0=0 时，极限核是平面的 `C/(epsilon rho)` 相互作用。恢复有限厚度 BARE 还需
保留密度的法向坐标。这里的二维指电荷及极化模型，电场在三维空间中传播。

#### 4.5 两块接地金属板的镜像屏蔽

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

平板 Green 函数及其镜像构造见
[Takae–Onuki (2013)](https://doi.org/10.1063/1.4821085)。
[Throckmorton–Vafek (2012)](https://doi.org/10.1103/PhysRevB.86.115447) 的电子
模型也采用平行板相互作用。这里 d 表示中面到每一块板的距离，两板间距为 2d。

对于中面上的点电荷，镜像级数也可展开为横向模：

$$
K_g(\rho;0,0)=\frac{C}{\epsilon}
\sum_{n=-\infty}^{\infty}
\frac{(-1)^n}{\sqrt{\rho^2+(2nd)^2}}
=\frac{2C}{\epsilon d}\sum_{j=0}^\infty
K_0\!\left(\frac{(2j+1)\pi\rho}{2d}\right).
$$

最低横向模控制远距离尾部 `exp[−pi rho/(2d)]/sqrt(rho)`，衰减尺度为 `2d/pi`。
完整的 gate 核包含全部横向模和两个电荷高度，因此同时描述空间衰减，以及电荷
接近上板或下板时响应的变化。

| 关键词 | 定义 |
|---|---|
| `WINT_KERNEL=DUAL_GATE` | 保留电荷高度的双接地板模型 |
| `WINT_BOUNDARY=SLAB_WS` | 所需的边界表示 |
| `WINT_GATE_D` | 必需，中面到每侧板的距离 d，Å；两板间距为 2d |
| `WINT_Z0` | 必需，两板中面及积分区间中心，Å |
| `WINT_EPS` | 均匀间隔介质的相对介电常数，默认 1，有限正值 |

当前实现将整个密度区间放在两板之间，因此要求 `d>Lz/2`。这一条件来自这里
采用的积分区域。d 很大时，算符趋于有限高度的 bare/epsilon。Green 函数按
上式的归一化进入每个两电子积分，多体 Hamiltonian 中的 1/2 用于成对计数。
单电子镜像势另行贡献于跃迁和 onsite Hamiltonian。

<a id="zh-use"></a>

### 5. PAW 投影与固定 Wannier 基底工作流

#### 选定的算符需一致作用于同一轨道密度

PAW 用平滑平面波部分和原子中心重建项表示波函数。为使所选相互作用一致地
作用于轨道，显式 WINT 路径重建完整跃迁密度：

$$
\rho_{ij}=\rho^{\mathrm{PW}}_{ij}
+\sum_a\left(\rho^{\mathrm{AE}}_{ij,a}
-\rho^{\mathrm{PS}}_{ij,a}\right).
$$

将算符与这一密度收缩，便包含平滑部分与原子中心项之间、以及不同原子之间的
交叉项。重建使用计算中已有的 PAW 分波和投影子。PAW 形式及其 VASP 实现见
[Blöchl (1994)](https://doi.org/10.1103/PhysRevB.50.17953) 和
[Kresse–Joubert (1999)](https://doi.org/10.1103/PhysRevB.59.1758)。

统一重建密度后，各屏蔽 Green 函数作用于同一轨道电荷分布。省略 WINT 标签
所选的原有 bare 路径采用其补偿和局域修正构造。在有限截断下比较两者，可检查
数值上的一致性；显式路径为新屏蔽模型提供共同的 bare 参考。

#### 从最终 Wannier 基底到相互作用文件

1. 完成电子结构计算，得到一致的 WAVECAR 与 k 网格。
2. 使用紧凑 SCDM 得到初始 Wannier 投影。
3. 固定最终 Wannier gauge；如需 SAWF，先由外部工具完成并导出兼容 VASP 的 WANPROJ。
4. 运行固定轨道的 `ALGO=2E4WA`，投影所选相互作用。

相互作用计算支持 CPU NCL、`ISPIN=1`、`NCORE=1`，包括非共线 SOC spinor。
若材料参考面位于笛卡尔 z=0，可在已有输入中加入以下双栅极设置：

```text
ALGO = 2E4WA
HFRCUT = -2
NCORE = 1
WINT_KERNEL = DUAL_GATE
WINT_BOUNDARY = SLAB_WS
WINT_Z0 = 0.0
WINT_GATE_D = 100.0
```

参考面按实际结构设定。本例中，每块板距离该面 100 Å，即 10 nm。比较不同
核时，保持 Wannier 基底和 PAW 数据集固定。VASP 按原有指标约定写出 VIJKL 和
VRijkl，并在文件头与 OUTCAR 中记录模型参数。后续将这些相互作用与同一
Wannier 基底中的 SAWF 单粒子 Hamiltonian 配套使用。

用于 fRG 时，epsilon、r0 和 kappa 定义输入模型的环境屏蔽或部分屏蔽。
后续 fRG 流将生成的目标低能带屏蔽过程应从输入屏蔽中排除。标量 epsilon 假定
均匀各向同性介质，各向异性基底或多层器件则需要合适的有效模型，或更一般的
Poisson Green 函数。当前极化片和栅极选项分别求解各自的静电问题；同时包含
二者的器件需要相应的组合 Green 函数。

<a id="zh-convergence"></a>

### 6. 收敛与物理解释

- **轨道模型：** 选带、Wannier gauge 和单粒子能带表示决定参与相互作用的电荷云。
- **法向真空与参考面：** I 应包含显著密度尾部。对 bare，在真空中移动区间端点可
  检查有限区域表示；相对固定材料移动极化片或栅极，则改变了物理环境。
- **空间范围：** 增大 BvK k 网格，扩大 WS 相互作用范围。
- **核分辨率：** 在固定范围内细化 WS_NGRID。WS_TOL 设定内部构造与积分容差，
  所得矩阵元的收敛情况决定其数值精度。
- **密度分辨率：** 联合收敛 VCUTOFF 与 NGX/NGY/NGZ。短程及 onsite 作用涉及
  大动量密度分量，单独增大 VCUTOFF 可能达到现有 FFT 网格的分辨率上限。

扩大密度网格后，若现有 CHGCAR 不再兼容，可在固定 WAVECAR 的相互作用投影中
采用 `ICHARG=12`，本工作流已检查这一用法。单粒子模型继续使用既定的 SAWF
Hamiltonian。在上述数值收敛检查中，屏蔽长度和介电常数始终定义所研究的物理
模型。

实现的有限基组算符已与独立 Green 函数计算、法向矩阵检验及原生投影和文件
输出比较。具体材料的精度由上述收敛检查确定。逐点为正的核作用于非负精确
轨道密度时，得到正的密度直接项；一般跃迁密度矩阵元不要求为正。

<a id="zh-install"></a>

### 安装与升级

当前支持的源码版本为 VASP 6.6.1 和 6.6.0，本地原生计算使用 6.6.1 CPU NCL。
使用者需持有 VASP 许可和自己的源码，并沿用已成功编译 VASP 的编译器、MPI
及 Wannier90 配置。对干净原版或可识别的仅 SCDM 源码执行：

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

随后换用新版安装包，检查、应用并重新编译。应在覆盖旧脚本与补丁，或在原目录
执行 `git pull` 之前完成回退，因为 `--reverse` 按其配套版本核对源码哈希。
可识别的干净原版或仅 SCDM 源码可直接使用当前安装脚本，脚本会保留源码备份。

部署时，将 `apply_patch.sh` 与统一的 `.patch` 文件放在同一目录。
`SOURCE_IDENTITIES.json` 记录支持的源码哈希，`manifest.sha256` 用于核对安装包
完整性。补丁沿用历史文件名，内容包含上文全部功能。本仓库分发第三方补丁；
VASP 源码、可执行文件和 PAW 数据集由使用者通过授权渠道获取。

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
