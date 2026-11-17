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
2. **`HFRCUT=-2`:** represent the bare interaction on the lattice-dependent
   Wigner–Seitz domain of the finite k mesh.
3. **Layered interactions:** retain the normal extent of the Wannier densities
   while introducing an ideal polarizable sheet or metallic gates.

The Wannier basis fixes the orbital charge distributions. The kernel specifies
their physical interaction, and the finite integration domain determines its
representation on the k mesh. Keeping WAVECAR and the final WANPROJ fixed
allows different interactions to be compared in the same basis, including
bases obtained from PBE or PBE/SOC calculations.

<a id="en-scdm"></a>

### 2. Band-compact SCDM: selecting the relevant Hilbert space without storing all bands

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

### 3. `HFRCUT=-2`: why a lattice-shaped Coulomb cutoff is useful

For an isolated two-dimensional material, the desired extended Hubbard
interaction couples Wannier charge clouds within one layer. A plane-wave
calculation represents that layer in a three-dimensional periodic cell with
vacuum along the normal. In a moiré cell, the normal period may be much shorter
than the in-plane interaction range. The Coulomb representation must account
for this geometry when the orbitals are projected onto a finite k mesh.

#### 3.1 What the original VASP treatment represents

The reciprocal bare kernel, `4πC/Q²`, has an integrable singularity at
`Q=0` whose contribution must be represented on a finite k mesh. Native
VASP uses a probe-charge Ewald correction for a regular unscreened mesh and
a Gaussian auxiliary-function treatment in other supported cases. The
`HFRCUT` parameter also provides spherical truncation. The native methods are
described in the official [HFRCUT](https://vasp.at/wiki/HFRCUT) and
[Coulomb singularity](https://vasp.at/wiki/Coulomb_singularity) documentation.

For a full, unreduced uniform mesh, the relevant real-space repetition is the
Born–von Karman (BvK) superlattice:

$$
\mathbf A_i=N_i\mathbf a_i,\qquad
N_k=N_1N_2N_3,\qquad V_B=N_k\Omega.
$$

The k mesh thus fixes a BvK repetition lattice. Ewald summation evaluates
electrostatics on that periodic lattice; the intended boundary conditions
determine which interaction the projected model represents.

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

The uniform term in the Poisson equation supplies a compensating background.
The Madelung correction cancels the regular image and background contribution
at the origin, while their separation-dependent contributions remain. In an
anisotropic cell, the resulting periodic potential can differ substantially
from the isolated `C/r` law at the distances relevant to an extended model.

This occurred in the tSnS calculations that motivated the patch: with
`HFRCUT=0`, the density-direct interactions developed a negative tail at large
separation. Analysis of the projected operator traced this tail to the
periodic Green function and its zero-mode correction. Recovering the
isolated-layer interaction therefore requires changing that finite-mesh
representation.

In native VASP 6.6.1, `HFRCUT=0` selects the periodic treatment above.
Positive values instead specify a spherical cutoff radius, and the documented
`HFRCUT=-1` chooses a sphere with the BvK volume:

$$
K_{\mathrm{sph}}(r)=\frac{C}{r}\Theta(R_c-r),\qquad
R_c=\left(\frac{3N_k\Omega}{4\pi}\right)^{1/3}.
$$

This radius uses the full uniform mesh without additional mesh reductions.
All negative values in the native code enter the automatic-sphere branch,
including `-2` in an unpatched source tree. The patch assigns `-2` to the WS
construction described below.

#### 3.2 Why a two-dimensional moiré geometry is demanding

A spherical cutoff introduces a geometric difficulty in a strongly
anisotropic cell. The vacuum needed to suppress orbital overlap may be short
compared with the in-plane moiré length, and the k mesh increases that contrast.
The relevant lengths are **`N_i a_i`**. With one normal k point, the normal BvK
period remains `Lz` as the in-plane range grows with the mesh.

For example, an approximately orthogonal `41 × 39 × 29 Å` cell with a
`6 × 6 × 1` mesh has BvK lengths near `246 × 233 × 29 Å`. Its automatic sphere
has radius about `73.5 Å`, compared with a normal WS half-width of about
`14.5 Å` and in-plane half-widths of about `117–123 Å`. It spans repeated
normal cells while cutting off useful in-plane separations. Sampling its
Fourier components periodizes the spherical kernel, so the overlap can retain
interactions with repeated layers.

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

For orthogonal vectors the domain is a box with half-widths `|A_i|/2`;
for oblique cells it follows the nearest-image geometry of the lattice.
The truncated kernel supplies every Fourier component, including its finite
head, in place of the periodic Ewald/FSG correction.

This WS-supercell construction follows
[Sundararaman and Arias (2013)](https://doi.org/10.1103/PhysRevB.87.165122),
whose analysis includes partially periodic systems. Matching the cutoff to
the repetition lattice makes it useful for both bulk and layered materials.
Their exchange-energy convergence analysis assumes appropriate localization;
for individual Wannier matrix elements, the density tails, available spatial
range and numerical resolution determine the remaining error.

| Input | Definition |
|---|---|
| `HFRCUT=-2` | Enable this patch's WS representation for the Wannier-integral workflow |
| `WINT_BOUNDARY=BVK_WS` | Full BvK WS domain; default boundary for explicit BARE/YUKAWA |
| `WS_NGRID` | Three-component kernel FFT grid; automatic by default |
| `WS_TOL` | Internal kernel construction/quadrature control; default `1e-8` |

With all `WINT_*` tags omitted, `HFRCUT=-2` uses the original bare-WS
evaluation. Setting `WINT_KERNEL=BARE` selects the reconstructed-PAW density
introduced below. Both evaluate the bare Coulomb law, using different
finite-basis density constructions.

<a id="en-slab"></a>

### 4. `SLAB_WS`: a finite-height representation of layered electrostatics

#### 4.1 Why keep two normal coordinates?

A metal plate or polarizable sheet singles out a height in the device.
Two charge pairs with the same relative displacement can then experience
different screening because their distances from that plane differ. The
interaction takes the form `V(rho−rho'; z,z')`, retaining both the source and
observation heights. Its Fourier representation consequently requires a
matrix in the normal coordinate.

`SLAB_WS` combines a finite in-plane WS window with a height-resolved Green
function. The construction below defines this combination, drawing on
Coulomb truncation and the Green functions of layered electrostatics.

#### 4.2 Specify the physical Green function, then its finite representation

Choose the infinite-plane interaction `V_infinity(s;z,z')` for the physical
environment: free space, a polarizable sheet or grounded plates. With
`W_parallel` denoting the WS cell of the two in-plane BvK vectors, the
finite-domain operator is

$$
I=[z_0-L_z/2,z_0+L_z/2],\qquad
K_W(\boldsymbol\rho,z;\boldsymbol\rho',z')
=\chi_{W_\parallel}(\boldsymbol\rho-\boldsymbol\rho')
V_\infty(\boldsymbol\rho-\boldsymbol\rho';z,z'),
\qquad z,z'\in I.
$$

The WS window limits the relative in-plane separation, while each density
has its own height coordinate within `I`. This interval defines the density
representation. The chosen Green function specifies the electrostatic
boundaries: unbounded in the normal direction for free space and the sheet,
or grounded at `z0±d` for the gate model.

The slab cutoff `Theta(zc−|z−z'|)/|r−r'|` of
[Ismail-Beigi (2006)](https://doi.org/10.1103/PhysRevB.73.233103) truncates the
relative normal displacement. Here the two densities instead occupy continuous
normal intervals, allowing the Green function to depend on their individual
heights as required by a sheet or gate.

Writing the analytic infinite-plane Green function as
`G_infinity(p;z,z')`, the planar Fourier transform of the windowed operator is

$$
\widetilde V_W(\mathbf q;z,z')
=\int\frac{d^2p}{(2\pi)^2}\,
\widetilde\chi_{W_\parallel}(\mathbf q-\mathbf p)
G_\infty(\mathbf p;z,z').
$$

The real-space window produces a momentum-space convolution. Each finite-WS
component at `q` therefore involves the infinite-plane Green function over
the momenta selected by the window. This convolution is included in the
operator evaluated by the implementation.

#### 4.3 Normal Fourier–Galerkin projection

The pair densities obtained from WAVECAR are available as Fourier
coefficients. To contract them with the open-normal Green function, represent
their height dependence in a Fourier basis on `I`. In centered coordinates,

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

The last expression converts the centered matrix to the global Fourier
convention of the supplied density. The phase makes `WINT_Z0` part of the
finite-interval projection even for bare Coulomb.

Even the free kernel, which depends only on `z−z'`, generally gives a dense
matrix on this finite interval. Its matrix elements integrate the open Green
function over the two height coordinates. A periodic Green function would
instead give a circular convolution and a diagonal Fourier representation.
The electrostatic boundary conditions thus enter through the Green function.

Including `C`, the planar Green function has units `eV Å²` and the normal
matrix has units `eV Å³`. VASP's Fourier and volume normalization converts
their contraction with the densities into energies. Factoring out `C` leaves
a geometric matrix with units `Å²`.

The spatial integrals separate through the identity

$$
\frac1{\sqrt{\rho^2+(x-x')^2}}
=\frac2{\sqrt\pi}\int_0^\infty dt\,
e^{-t^2\rho^2}e^{-t^2(x-x')^2}.
$$

The planar factor is integrated over the WS window, while the normal matrix
elements reduce to one-dimensional moments. Corresponding representations
handle sheet-induced and gate-image terms. These moments retain the orbital
thickness and treat the singularity analytically in the normal integration.

#### 4.4 Physical assumptions and the choice of reference plane

`WINT_Z0` is an absolute Cartesian coordinate along the layer normal, in Å.
It centers the density interval for bare/Yukawa, places the polarizable plane
for RK_SHEET and sets the plate midpoint for DUAL_GATE. The geometric center
of the layer gives a starting estimate; the effective polarization plane and
gate positions should reflect the physical model.

Place the ends of `I` in vacuum so the interval contains the layer and its
significant density tails. Moving the endpoints within sufficient vacuum
tests the remaining tail and numerical errors in the bare integral. If the
layer crosses the periodic coordinate origin, unwrap its coordinates before
determining a continuous layer and its center.

The in-plane window fixes the available interaction range. Convergence
toward the infinite-area result requires increasing that range alongside
checks of normal vacuum and Fourier resolution. The physical screening
parameters specify the model being converged. `BVK_WS` also applies to layers:
its bare integrals can agree with SLAB_WS in the interior when the relevant
normal separations lie within the BvK cutoff and density tails are small.
The two domains generally produce different complete finite-basis matrices.

| Input / condition | Meaning |
|---|---|
| `WINT_BOUNDARY=SLAB_WS` | In-plane BvK WS window plus open-normal density interval |
| `WINT_Z0` | Required interval center for finite-height models; also the sheet/gate reference plane when applicable |
| `WINT_NORMAL` | Optional Cartesian normal, default along the third lattice vector |
| Lattice and k mesh | Third vector perpendicular to the plane; one normal k point; oblique in-plane cells are supported |
| `WS_NGRID` | First two components refine the plane kernel FFT; third component is unused by this boundary |

RK_SHEET, RK2D and DUAL_GATE use this boundary explicitly. BARE and YUKAWA
support both boundaries. For RK2D, integration over height removes the
dependence on Z0.

<a id="en-kernels"></a>

### 5. Choosing a physical interaction

The following static kernels describe different electrostatic environments
for a fixed Wannier basis. The dielectric parameter `epsilon` is the relative
permittivity, while `epsilon_0` is the vacuum permittivity. Their screening
lengths set physical scales in the model; FFT parameters control numerical
resolution.

#### 5.1 `BARE`: the repulsion of the Wannier charge clouds

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

#### 5.2 `YUKAWA`: homogeneous three-dimensional screening

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

#### 5.3 `RK_SHEET`: a polarizable sheet with finite-height Wannier charges

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

#### 5.4 `RK2D`: the strict zero-thickness approximation

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

### 6. PAW projection and practical use

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

### 7. What must be converged for a physical interpretation?

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
2. **`HFRCUT=-2`：** 用与有限 k 网格对应超晶格匹配的 WS 区域表示裸库仑相互作用。
3. **层状相互作用：** 保留 Wannier 密度的法向厚度，并引入理想极化片或金属栅极。

Wannier 基底确定轨道电荷分布，物理核定义电荷之间的相互作用，有限积分区域
规定它在 k 网格上的表示。固定 WAVECAR 和最终 WANPROJ 后，便可在同一组轨道
中比较不同相互作用；这组轨道可以来自 PBE 或 PBE/SOC 计算。

<a id="zh-scdm"></a>

### 2. Band-compact SCDM：让工作空间对应真正保留的能带

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

### 3. `HFRCUT=-2`：为什么需要随晶格形状变化的库仑截断

为孤立二维材料构造扩展 Hubbard 模型时，我们需要同一材料层内 Wannier 电荷云
之间的相互作用。平面波计算用三维周期晶胞表示这一层，并在法向加入真空。
对 moiré 晶胞，法向周期可能远短于面内相互作用范围。把轨道投影到有限 k 网格
后，库仑核的表示需要与这一几何相适应。

#### 3.1 原生 VASP 的处理对应什么物理对象

裸库仑核的倒空间形式 `4πC/Q²` 在 `Q=0` 具有可积奇点，需要在有限 k 网格上
表示其贡献。原生 VASP 对规则、未屏蔽网格采用 probe-charge Ewald 修正，其他
支持的情况使用高斯辅助函数处理。`HFRCUT` 还提供球形截断选项。原生方法见
[HFRCUT](https://vasp.at/wiki/HFRCUT) 和
[Coulomb singularity](https://vasp.at/wiki/Coulomb_singularity) 文档。

对于没有额外约化的完整均匀网格，实空间周期性对应 Born–von Karman 超晶格：

$$
\mathbf A_i=N_i\mathbf a_i,\qquad
N_k=N_1N_2N_3,\qquad V_B=N_k\Omega.
$$

k 网格由此确定 BvK 重复晶格。Ewald 求和用于计算该周期晶格上的静电作用，
而所选边界条件决定投影后的低能模型对应哪一种相互作用。

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

Poisson 方程中的均匀项提供补偿背景。Madelung 修正抵消原点处镜像与背景的
正则贡献，有限距离上随间距变化的贡献则仍然存在。在高度各向异性的晶胞中，
由此得到的周期电势可能在扩展模型关心的距离上明显偏离孤立的 `C/r` 相互作用。

本项目的 tSnS 计算便出现了这一现象：使用 `HFRCUT=0` 时，密度直接项在远距离
出现负尾部。对投影算符的分析表明，负尾部来自周期 Green 函数及零模修正。
为了得到孤立材料层中的相互作用，需要调整这一步的有限网格表示。

原生 VASP 6.6.1 中，`HFRCUT=0` 选择上述周期处理。正值则指定球形截断半径，
文档中的 `HFRCUT=-1` 选择与 BvK 体积相同的球：

$$
K_{\mathrm{sph}}(r)=\frac{C}{r}\Theta(R_c-r),\qquad
R_c=\left(\frac{3N_k\Omega}{4\pi}\right)^{1/3}.
$$

该半径使用无额外约化的完整均匀网格。原生源码中的所有负值均进入自动球形
截断分支，包括未打补丁时的 `-2`。本补丁将 `-2` 指定为下文的 WS 构造。

#### 3.2 二维 moiré 几何为什么特别敏感

在高度各向异性的晶胞中，球形截断还面临几何上的困难。抑制层间轨道重叠所需的
真空可能短于面内 moiré 尺度，而 k 网格会进一步放大这一差异。相关长度为
**`N_i a_i`**：法向仅一个 k 点时，BvK 法向周期始终为 `Lz`，面内范围则随网格
增加而扩大。

例如，约 `41 × 39 × 29 Å` 的近似正交晶胞配合 `6 × 6 × 1` 网格，得到约
`246 × 233 × 29 Å` 的 BvK 尺寸。自动球形半径约为 `73.5 Å`，相比之下，
法向 WS 半宽约为 `14.5 Å`，面内半宽约为 `117–123 Å`。球在法向跨越重复晶胞，
同时截掉一部分有用的面内长距离。离散 Fourier 采样使球形核周期重复，这种
重叠因而可能保留与重复材料层的相互作用。

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

正交晶格的 WS 区域是各方向半宽为 `|A_i|/2` 的长方体；斜晶格则由最近镜像
几何确定。截断核统一给出所有 Fourier 分量及有限零模，替代原来的 Ewald/FSG
零模修正。

这一 WS 超胞构造来自
[Sundararaman–Arias (2013)](https://doi.org/10.1103/PhysRevB.87.165122)，其分析
包含部分周期体系。截断形状与重复晶格匹配，因此适用于体材料和层状材料。
文中的交换能收敛分析以适当的局域化为前提；对单个 Wannier 矩阵元，密度尾部、
可用空间范围和数值分辨率决定其剩余误差。

| 输入 | 定义 |
|---|---|
| `HFRCUT=-2` | 在 Wannier 积分工作流中启用本补丁的 WS 表示 |
| `WINT_BOUNDARY=BVK_WS` | 完整 BvK WS 区域；显式 BARE/YUKAWA 的默认边界 |
| `WS_NGRID` | 三分量核 FFT 网格；默认自动选择 |
| `WS_TOL` | 核构造及内部积分控制，默认 `1e-8` |

省略全部 `WINT_*` 标签时，`HFRCUT=-2` 使用原有 bare-WS 求值方式。设置
`WINT_KERNEL=BARE` 则选择下文的重建 PAW 密度。两者计算同一裸库仑律，采用
不同的有限基组密度构造。

<a id="zh-slab"></a>

### 4. `SLAB_WS`：保留层状静电问题中的两个高度坐标

#### 4.1 为什么需要分别保留 z 和 z'？

金属板或极化片在器件中选定了一个高度。两对电荷即使相对位移相同，也可能因
距离该平面的远近不同而受到不同屏蔽。因此，相互作用具有 `V(rho−rho'; z,z')`
的形式，需要同时保留源与观察点的高度；相应的 Fourier 表示在法向是一个矩阵。

`SLAB_WS` 将有限面内 WS 窗口与高度分辨的 Green 函数结合起来。下文定义这一
构造，其中的库仑截断和层状 Green 函数分别以已有方法为基础。

#### 4.2 先定义物理环境，再定义有限表示

先根据物理环境选择无限平面相互作用 `V_infinity(s;z,z')`，例如自由空间、
极化片或接地金属板。以 `W_parallel` 表示两根面内 BvK 矢量的 WS 区域，有限
算符定义为

$$
I=[z_0-L_z/2,z_0+L_z/2],\qquad
K_W(\boldsymbol\rho,z;\boldsymbol\rho',z')
=\chi_{W_\parallel}(\boldsymbol\rho-\boldsymbol\rho')
V_\infty(\boldsymbol\rho-\boldsymbol\rho';z,z'),
\qquad z,z'\in I.
$$

WS 窗口限制相对面内位移，两个密度则在 `I` 中各自保留高度坐标。这一区间
规定密度的表示范围，静电边界由所选 Green 函数确定：自由空间和极化片模型在
法向无界，双栅极模型在 `z0±d` 处接地。

[Ismail-Beigi (2006)](https://doi.org/10.1103/PhysRevB.73.233103) 的 slab 截断
`Theta(zc−|z−z'|)/|r−r'|` 限制相对法向位移。这里则把两个密度分别放在连续
法向区间内，使 Green 函数能够按极化片或栅极的要求依赖各自高度。

将无限平面的解析 Green 函数记为 `G_infinity(p;z,z')`，带窗口算符的面内
Fourier 变换为

$$
\widetilde V_W(\mathbf q;z,z')
=\int\frac{d^2p}{(2\pi)^2}\,
\widetilde\chi_{W_\parallel}(\mathbf q-\mathbf p)
G_\infty(\mathbf p;z,z').
$$

实空间窗口在动量空间产生卷积。每个有限 WS 分量都包含窗口选定的各动量贡献，
实现中对这一完整的带窗口算符进行积分。

#### 4.3 法向 Fourier–Galerkin 投影

从 WAVECAR 得到的跃迁密度可用 Fourier 系数表示。为使其与开放法向 Green
函数收缩，在 `I` 上用 Fourier 基底表示高度依赖。取中心化坐标：

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

最后一式将中心化矩阵转换到原密度的全局 Fourier 约定。这个相位使 `WINT_Z0`
在 bare 情形下也参与有限区间投影。

即使自由核只依赖 `z−z'`，在这一有限区间上通常也会得到稠密矩阵。矩阵元对
两个高度坐标上的开放 Green 函数积分；周期 Green 函数则对应循环卷积及对角
Fourier 表示。因此，静电边界条件由 Green 函数本身给定。

计入 C 后，面内 Green 函数的量纲为 `eV Å²`，法向矩阵为 `eV Å³`。与密度
收缩后，VASP 的 Fourier 和体积归一化将其转换为能量。提取 C 后，几何矩阵的
量纲为 `Å²`。

空间积分可通过以下恒等式分离：

$$
\frac1{\sqrt{\rho^2+(x-x')^2}}
=\frac2{\sqrt\pi}\int_0^\infty dt\,
e^{-t^2\rho^2}e^{-t^2(x-x')^2}.
$$

面内因子在 WS 窗口内积分，法向矩阵元则化为一维矩积分。极化片诱导项和栅极
镜像项也有相应表示。通过这些矩积分，计算保留轨道厚度，并在法向积分中解析
处理奇点。

#### 4.4 参考面及物理假设

`WINT_Z0` 是沿法向的绝对笛卡尔坐标，单位为 Å。它为 bare/Yukawa 确定密度
区间中心，为 RK_SHEET 确定极化片位置，为 DUAL_GATE 确定两板中面。材料层的
几何中心可作为初始估计；有效极化面和栅极位置应依据所研究的物理模型选取。

将 `I` 的两端放在真空中，使区间包含材料及其显著密度尾部。在足够真空内移动
端点，可检验 bare 积分中剩余的尾部误差和数值误差。若材料层跨越周期坐标原点，
先将原子坐标展开为连续的一层，再确定其中心。

面内窗口确定可计算的相互作用范围。趋近无限面积极限时，需要扩大这一范围，
同时检查法向真空和 Fourier 分辨率；物理屏蔽参数则定义所收敛的模型。
BVK_WS 也适用于层状体系。当主要法向间距位于 BvK 截断内，且密度尾部足够小，
其 bare 积分可在内部与 SLAB_WS 一致；两个区域的完整有限基组矩阵通常不同。

| 输入/条件 | 含义 |
|---|---|
| `WINT_BOUNDARY=SLAB_WS` | 面内 BvK WS 窗口，加上各密度的开放法向区间 |
| `WINT_Z0` | 有限厚度模型必需的区间中心；相应模型中同时是极化片或两板参考面 |
| `WINT_NORMAL` | 可选法向，默认沿第三晶格矢量 |
| 晶格和 k 网格 | 第三矢量垂直面内；法向一个 k 点；面内可以为斜晶格 |
| `WS_NGRID` | 前两项细化面内核 FFT；本边界不使用第三项 |

RK_SHEET、RK2D 和 DUAL_GATE 显式使用此边界，BARE 和 YUKAWA 支持两种边界。
对 RK2D，法向积分消除了对 Z0 的依赖。

<a id="zh-kernels"></a>

### 5. 选择物理相互作用

以下静态核为同一 Wannier 基底描述不同的静电环境。`epsilon` 是相对介电常数，
`epsilon_0` 是真空介电常数。屏蔽长度确定模型中的物理尺度，FFT 参数控制计算
分辨率。

#### 5.1 `BARE`：Wannier 电荷云之间的裸排斥

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

#### 5.2 `YUKAWA`：均匀三维介质中的指数屏蔽

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

#### 5.3 `RK_SHEET`：极化局限于平面，轨道电荷保留厚度

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

#### 5.4 `RK2D`：严格零厚度近似

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

### 6. PAW 投影与实际使用

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

### 7. 哪些收敛问题决定物理解释？

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
