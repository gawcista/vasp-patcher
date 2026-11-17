# VASP patcher: compact-band SCDM and Wannier interactions

[English](#english) · [中文说明](#中文说明) · [References / 参考文献](#references--参考文献)

## English

This third-party patch for VASP 6.6.1 combines compact-band SCDM working arrays,
Wigner–Seitz (WS) Coulomb truncation and fixed-orbital Wannier interaction
projection with bare, Yukawa, Rytova–Keldysh and dual-gate models. Exact source
profiles also support deployment to VASP 6.6.0; local native execution has been
validated with the 6.6.1 CPU NCL executable.

A valid VASP license and your own VASP source distribution are required. This
repository supplies a patch with the context needed to apply it, rather than a
VASP source distribution. Binaries, PAW datasets and simulation inputs are
obtained separately under their applicable licenses.

### Package contents

| File | Purpose |
|---|---|
| `apply_patch.sh` | Check source identities, apply the unified patch, report status and reverse it |
| `vasp-6.6.1-scdm-compact-bands.patch` | Combined source changes; the historical filename covers all features below |
| `SOURCE_IDENTITIES.json` | Exact vendor and patched hashes for both supported source profiles |
| `manifest.sha256` | Integrity hashes for the deployment package |
| `rollback/` | Matching earlier installer/patch pairs for upgrading existing installations |

The unified payload changes eight paths:

```text
src/mlwf.F           src/fock_glb.F       src/fock.F
src/dmft.F           src/.objects        src/ws_coulomb.F
src/wint_paw.F       src/wint_layered.F
```

### Installation and upgrades

Use the compiler, MPI, numerical libraries and Wannier90 configuration already
working for your VASP installation. From this repository:

```bash
sha256sum -c manifest.sha256
VASP_WINT_SRC=/absolute/path/to/vasp-source
bash apply_patch.sh --status "$VASP_WINT_SRC"
bash apply_patch.sh --check "$VASP_WINT_SRC"
bash apply_patch.sh --apply "$VASP_WINT_SRC"
make -C "$VASP_WINT_SRC" ncl DEPS=1 -j4
bash apply_patch.sh --status "$VASP_WINT_SRC"
```

The `--check` operation previews the source change. The recognized states are:

| Status | Meaning and next action |
|---|---|
| `clean-base` | Recognized vendor source: check and apply |
| `scdm-only` | Accepted compact-SCDM source: check and apply upgrades the remaining paths |
| `applied` | Complete current payload: application is idempotent; rebuild if required |
| `modified/unknown` | Source hashes do not match a recognized complete state |

For an earlier combined WS or Yukawa installation, use its matching archived
installer to reverse the old payload first. The [rollback guide](rollback/README.md)
identifies the available packages. For example, **after its status reports
`applied`**, the original WS package is reversed with:

```bash
bash rollback/ws-original-20261110/apply_patch.sh --reverse "$VASP_WINT_SRC"
bash apply_patch.sh --check "$VASP_WINT_SRC"
bash apply_patch.sh --apply "$VASP_WINT_SRC"
make -C "$VASP_WINT_SRC" ncl DEPS=1 -j4
```

Copying new patch files and running `make` does not apply their source changes.
The installer verifies the payload and each affected source file, uses
`git apply` or zero-fuzz `patch`, and keeps timestamped source backups.
Unrelated files remain outside its affected-path list.

To restore the recognized vendor sources after the current installation:

```bash
bash apply_patch.sh --reverse "$VASP_WINT_SRC"
```

Reversal restores all eight vendor-path states, including removal of the three
new source modules. Preserve the old executable if it is needed for previous
calculations. The unified payload includes the compact-SCDM changes, so no
additional SCDM patch is applied after installation.

### 1. Compact-band SCDM

The SCDM method selects columns of a density matrix to construct localized
orbitals; its periodic and entangled-band formulations are described in
Refs. [1–3](#references--参考文献). This patch changes the storage and indexing of
VASP's existing SCDM implementation.

The `exclude_bands` mask from Wannier90 setup defines a compact working space.
For `b` retained bands and `R` real-space/spinor samples, the dominant SCDM
buffers use `(R,b)` and `(b,R)`, and its transformation uses
`(b,NUM_WANN,NKPTS,ISPIN)`. Wavefunctions, eigenvalues, cutoff weights and final
AMN rows share ascending original-band order. Global wavefunction/WAVECAR band
dimensions retain their existing sizes; the MMN temporary matrix retains its
original full-band allocation.

Use the normal VASP/Wannier90 SCDM inputs (`LSCDM`, `NUM_WANN` and the
Wannier90 `exclude_bands` selection). Cutoff behavior in the compact space is:

- `CUTOFF_TYPE=none` assigns unit weights to retained bands.
- With exclusions and exactly `NUM_WANN` retained bands, omitting a cutoff also
  assigns unit weights.
- A larger retained space requires an explicit cutoff choice: ERFC/Gaussian
  require `CUTOFF_MU` and `CUTOFF_SIGMA`; Fermi requires `CUTOFF_MU`; `num_wann`
  requires `CUTOFF_SIGMA`. `none` explicitly selects unit weights.
- Explicit `CUTOFF_MU` and `CUTOFF_SIGMA` without `CUTOFF_TYPE` preserve the
  historical ERFC interpretation. Inputs without exclusions retain the native
  estimation behavior.

This optimization reduces selected-band SCDM storage; total VASP memory also
contains the other native allocations. Near-tie rank-revealing QR decisions
can depend on floating-point/parallel layout, so compare downstream data in a
verified common Wannier basis.

### 2. Fixed-orbital interaction workflow

A typical workflow is:

1. Prepare a converged electronic structure and matching WAVECAR/k mesh.
2. Generate the initial Wannier interface with compact-band SCDM.
3. Construct the final symmetry-adapted Wannier functions externally and export
   their transformation as a VASP-compatible `WANPROJ`.
4. Run VASP with `ALGO=2E4WA` to write `VIJKL` and `VRijkl` in that frozen basis.

SAWF localization and WANPROJ conversion are separate tools. The patch consumes
the resulting WANPROJ; it does not implement SAWF. The native tensor data
columns, index conventions and q-to-R phases are preserved; output headers
identify the chosen interaction. The unpatched tensor definition is documented
by the [VASP VRijkl reference](https://vasp.at/wiki/VRijkl).

With all `WINT_*` tags omitted, `HFRCUT=-2` retains the earlier bare-WS
evaluator. An explicit `WINT_KERNEL` selects a common full-PAW-density evaluator:

$$
\rho_{ij}=\rho^{\mathrm{PW}}_{ij}
 +\sum_a\left(\rho^{\mathrm{AE}}_{ij,a}-\rho^{\mathrm{PS}}_{ij,a}\right).
$$

The selected operator acts on this fixed density, including smooth–augmentation
and interatomic cross terms. Here full PAW refers to the supplied frozen-core,
finite-partial-wave PAW representation. No extra all-electron solid calculation
is needed. Explicit BARE and omitted-tag bare use the same physical Coulomb
law through different numerical evaluations and require their own convergence
assessment. Changing the interaction reuses the fixed orbitals and matching
WANPROJ. This interface is independent of the native hybrid/Thomas–Fermi and
cRPA screening branches.

### 3. Integration boundaries

`HFRCUT=-2` selects WS treatment; `WINT_KERNEL` selects the physical model;
`WINT_BOUNDARY` selects its finite-domain representation.

| Boundary | Definition | Supported models |
|---|---|---|
| `BVK_WS` (default) | Full real-space WS cell of the k-mesh Born–von Karman superlattice | BARE, YUKAWA |
| `SLAB_WS` | In-plane BvK WS window and an open normal integration interval for each density | All five models |

The `BVK_WS` construction follows Sundararaman–Arias [4], including its
application to partially periodic systems. It therefore also applies to layers.
The Cartesian metric defines the WS domain, including oblique cells.

`SLAB_WS` is this implementation's name for a specific combination of the
in-plane WS window and normal Fourier-Galerkin matrix. Its two normal
coordinates range independently over
`[z0−Lz/2,z0+Lz/2]`, where `Lz` is the DFT normal cell length. Periodic densities
are unwrapped into that interval. The finite-height operator couples normal
Fourier components at fixed in-plane momentum. This construction uses the
low-dimensional Coulomb-cutoff and mixed-space Green-function principles of
Refs. [5,8], but differs from the relative-height cutoff
`Theta(zc−|z−z'|)/|r−r'|` in Ref. [5].

`WINT_Z0` is a Cartesian coordinate in Å along the chosen normal. For slab
BARE/YUKAWA it centers the numerical integration interval. For RK_SHEET it also
locates the polarizable sheet; for DUAL_GATE it is the midpoint between plates.
Choose an interval containing the material and its density tails, with the
boundaries in vacuum. A geometric atomic-layer center is a possible starting
choice; the physical polarization center and actual device geometry may
require another model choice. RK2D integrates out height and rejects this tag.

### 4. Physical kernels

Let `C=e²/(4π epsilon_0)`, `epsilon=WINT_EPS`,
`q=|q_parallel|`, `rho=|r_parallel−r'_parallel|` and `x=z−z0`.
VASP reports the projected interaction energies in eV. The following mixed-space
expressions define the infinite-plane models before the common finite in-plane
WS window is applied.

| Selector | Physical picture |
|---|---|
| `BARE` | Vacuum Coulomb interaction between finite-height charge distributions |
| `YUKAWA` | Homogeneous exponential screening with length `1/kappa` |
| `RK_SHEET` | One polarizable insulating sheet; retain the Wannier charge heights |
| `RK2D` | Strictly planar RK interaction on height-integrated pair densities |
| `DUAL_GATE` | Finite-height charges between two symmetric grounded conducting plates |

**Bare and Yukawa** use

$$
K_{\rm bare}(r)=\frac{C}{r},\qquad
K_Y(r)=\frac{C\,e^{-\kappa r}}{\epsilon r}.
$$

The Yukawa functional form has the historical origin in Ref. [6]; its use here
is a phenomenological screening model. It is also the Green function of
`(-nabla²+kappa²)K=4πC delta(r)/epsilon`. Its decay differs from the field
screening of metal plates.

**Finite-height sheet** uses the independently derived Green function of an
ideal polarizable plane at `z0`:

$$
G_s(q;z,z')=\frac{2\pi C}{\epsilon q}
\left[e^{-q|z-z'|}
-\frac{r_0q}{1+r_0q}
 e^{-q(|z-z_0|+|z'-z_0|)}\right].
$$

At `z=z'=z0`, it reduces to the RK model in Ref. [7]. The finite-height
extension follows the corresponding sheet-polarization boundary condition;
it is not an isotropic radial substitution.

**Strict RK2D** uses

$$
G_{2D}(q)=\frac{2\pi C}{\epsilon q(1+r_0q)},\qquad
K_{2D}(\rho)=\frac{\pi C}{2\epsilon r_0}
\left[H_0(\rho/r_0)-Y_0(\rho/r_0)\right].
$$

Each pair density is first integrated over height. This potential has a
logarithmic short-distance form and a `C/(epsilon rho)` distant tail [7].
`r0` is defined by the displayed Fourier convention.

**Dual gate** places plates at `x=−d,+d` and uses the finite-height Dirichlet
Green function [8]:

$$
G_g(q;x,x')=\frac{4\pi C}{\epsilon q}
\frac{\sinh[q(x_<+d)]\sinh[q(d-x_>)]}{\sinh(2qd)}.
$$

Here `x_<`/`x_>` are the smaller/larger of `x,x'`. At the midplane,
`G_g(q;0,0)=2πC tanh(qd)/(epsilon q)`, with the finite limit
`G_g(0;0,0)=2πCd/epsilon`. The related planar gated model is discussed in Ref. [9].
`WINT_GATE_D` is the distance to **each** plate, so the plate separation is `2d`.
The present integration domain requires `d>Lz/2`; this is a numerical-domain
restriction. The operator supplies the two-electron interaction; one-electron
image potentials are a separate Hamiltonian contribution.

At `epsilon=1` and the same domain, `YUKAWA(kappa=0)`, `RK_SHEET(r0=0)` and
`DUAL_GATE(d→infinity)` recover finite-height BARE. For general epsilon these
limits are BARE/epsilon. `RK2D(r0=0)` instead gives planar
`C/(epsilon rho)` on the height-integrated densities. The finite-WS head and
all nonzero components are constructed from the same truncated operator.

The layered evaluator separates plane and normal coordinates through a
Gaussian representation. Smooth plane terms use the WS FFT, short-range plane
factors are analytic, and normal matrix elements use one-dimensional moment
quadrature. The same reconstructed PAW densities feed the native tensor
contraction and writer.

### 5. Inputs and supported calculation path

The interaction implementation supports CPU NCL, `ALGO=2E4WA`, `HFRCUT=-2`,
`ISPIN=1` and `NCORE=1`, including noncollinear SOC spinors. Layered calculations
require a single normal k point and the third lattice vector perpendicular to
the in-plane lattice. Oblique in-plane cells are supported.

| Tag | Definition / default |
|---|---|
| `WINT_KERNEL` | BARE, YUKAWA, RK_SHEET, RK2D or DUAL_GATE; physical default BARE; omit all WINT tags to retain the legacy evaluator |
| `WINT_BOUNDARY` | BVK_WS by default; explicitly set SLAB_WS for RK_SHEET/RK2D/DUAL_GATE |
| `WINT_KAPPA` | Required for YUKAWA, Å⁻¹, finite and nonnegative |
| `WINT_EPS` | Relative dielectric constant, default 1 for screened models, finite and positive; explicit BARE rejects it |
| `WINT_R0` | Required for RK_SHEET/RK2D, Å, finite and nonnegative |
| `WINT_GATE_D` | Required for DUAL_GATE, Å, distance to each plate; must exceed Lz/2 |
| `WINT_Z0` | Required for finite-height SLAB_WS models, Å; omitted for RK2D |
| `WINT_NORMAL` | Optional Cartesian normal, defaults to the third lattice-vector direction; supplied direction must be parallel to it |
| `WS_NGRID` | Automatic kernel FFT grid by default; three components, with only the first two used by SLAB_WS |
| `WS_TOL` | Internal kernel construction/quadrature tolerance, default `1e-8` |

For example, add the following to an existing fixed-orbital 2E4WA INCAR if the
material reference plane is at Cartesian `z=0`:

```text
WINT_KERNEL = DUAL_GATE
WINT_BOUNDARY = SLAB_WS
WINT_Z0 = 0.0
WINT_GATE_D = 100.0
```

Here `d=100 Å` is an illustrative 10 nm distance to each gate. Set the reference
plane and dielectric parameters for the actual device. Unsupported or
inconsistent native modes and irrelevant model parameters stop explicitly.

### 6. Numerical resolution and validation

Converge three different numerical choices: the kernel grid (`WS_NGRID`), the
reconstructed-density response basis (`VCUTOFF` with `NGX/NGY/NGZ`) and the
BvK range/k mesh (in-plane for SLAB_WS). `WS_TOL` controls the kernel construction and internal
quadrature; it is not a final tensor error bound. A larger response cutoff on an
unchanged density grid can encounter a grid ceiling. Check the effective grids
and diagnostics in OUTCAR alongside changes in the complete interaction tensor.

When a changed density FFT grid is incompatible with the old CHGCAR, the tested
`ICHARG=12` option is for fixed-WAVECAR interaction projection. Continue using
the established SAWF Hamiltonian for the one-electron model; this option does
not establish equivalence of companion one-electron/FCIDUMP quantities.

The retained implementation checks include 32 numerical tests under each of
GNU Fortran and Intel ifx, 76 source/input/audit tests, 21 final-binary native
calculations, and 12 upgrade/reversal workflows across 6.6.0/6.6.1 and Git/non-Git
targets. They check independent physical references, complex pair densities,
zero-screening limits, dielectric scaling, MPI/KPAR and legacy output
regression. Material-specific response/grid/range convergence and production
CPU/memory use remain application-level checks. The PAW and layered modules
each cap their additional owned workspace at 512 MiB; native caller-owned
wavefunction and density arrays contribute separately to whole-process memory.
Screening parameters for an
fRG model should be chosen consistently with the low-energy screening that
will be treated later.

## 中文说明

本补丁面向 VASP 6.6.1，将选带紧凑存储的 SCDM、Wigner–Seitz（WS）库仑截断，
以及固定 Wannier 轨道下的裸库仑、Yukawa、Rytova–Keldysh 和双栅极相互作用投影
集成在一个部署包中。VASP 6.6.0 具有精确的源码部署兼容配置；本地原生运行验证
使用的是 6.6.1 CPU NCL 可执行文件。

使用者需要有效的 VASP 许可及自行取得的 VASP 源码。本仓库提供第三方补丁和应用
补丁所需的上下文，不提供完整 VASP 源码、可执行文件、PAW 数据库或计算输入。

### 文件和安装

`apply_patch.sh` 与 `vasp-6.6.1-scdm-compact-bands.patch` 是部署所需的两个文件。
统一补丁的历史文件名保留，但实际包含 SCDM、WS 与全部新相互作用模型；修改范围
是上面列出的八条源码路径。`SOURCE_IDENTITIES.json` 记录两个版本的精确源码哈希，
`manifest.sha256` 用于校验部署包，`rollback/` 保存匹配的旧版安装脚本与补丁。

在本仓库目录中，使用原来已经能够编译 VASP 的编译器、MPI、数学库和 Wannier90
配置，逐步执行：

```bash
sha256sum -c manifest.sha256
VASP_WINT_SRC=/absolute/path/to/vasp-source
bash apply_patch.sh --status "$VASP_WINT_SRC"
bash apply_patch.sh --check "$VASP_WINT_SRC"
bash apply_patch.sh --apply "$VASP_WINT_SRC"
make -C "$VASP_WINT_SRC" ncl DEPS=1 -j4
bash apply_patch.sh --status "$VASP_WINT_SRC"
```

`clean-base` 表示原始源码，`scdm-only` 表示已安装受支持的紧凑 SCDM；二者均可
检查后应用新补丁。`applied` 表示当前完整补丁已安装，重复应用保持源码不变。
`modified/unknown` 表示受影响源码不符合完整的已知配置。

已有旧 WS 或 Yukawa 补丁时，先按 [回退说明](rollback/README.md)找到状态为
`applied` 的旧安装包，再用它回退，之后应用新包。例如，若匹配的是原始 WS 包：

```bash
bash rollback/ws-original-20261110/apply_patch.sh --reverse "$VASP_WINT_SRC"
bash apply_patch.sh --check "$VASP_WINT_SRC"
bash apply_patch.sh --apply "$VASP_WINT_SRC"
make -C "$VASP_WINT_SRC" ncl DEPS=1 -j4
```

仅覆盖补丁文件再运行 `make` 不会修改源码。脚本会校验补丁与源码哈希，使用
`git apply` 或零模糊匹配的 `patch`，并保留带时间戳的源码备份。升级可以在原
源码目录完成；新统一补丁已经包含紧凑 SCDM，无需另外叠加 SCDM 补丁。

当前补丁的完整回退命令为：

```bash
bash apply_patch.sh --reverse "$VASP_WINT_SRC"
```

它恢复所有受影响的原始源码状态，并删除三个新增模块。需要复现旧计算时，另行
保留旧可执行文件。

### 1. 紧凑选带 SCDM

SCDM 通过密度矩阵选列构造局域轨道，其原始、周期及纠缠能带方法见参考文献
[1–3](#references--参考文献)。本补丁优化的是 VASP 已有 SCDM 的内存存储与选带
索引实现。

Wannier90 的 `exclude_bands` 决定保留的 `b` 条能带；若实空间与旋量采样总数为
`R`，主要 SCDM 工作数组采用 `(R,b)`、`(b,R)`，变换矩阵采用
`(b,NUM_WANN,NKPTS,ISPIN)`。波函数、本征值、截断权重和输出 AMN 行均按原始
能带编号的升序排列。全局波函数与 WAVECAR 的能带维数保持原有大小；MMN 临时
矩阵仍保留完整能带分配。

使用正常的 VASP/Wannier90 SCDM 输入，即 `LSCDM`、`NUM_WANN` 和
`exclude_bands`。紧凑子空间中的权重规则是：

- `CUTOFF_TYPE=none` 对所有保留态赋单位权重。
- 启用排除能带且恰好保留 `NUM_WANN` 个态时，省略截断参数也使用单位权重。
- 保留态多于 `NUM_WANN` 时需显式选择截断：ERFC/Gaussian 需要
  `CUTOFF_MU` 与 `CUTOFF_SIGMA`，Fermi 需要 `CUTOFF_MU`，`num_wann` 需要
  `CUTOFF_SIGMA`；选择 `none` 则使用单位权重。
- 同时给出 `CUTOFF_MU`、`CUTOFF_SIGMA` 而省略 `CUTOFF_TYPE` 时，保留原有
  ERFC 解释。不排除能带时保留原生自动估计行为。

实际总内存还包含 VASP 的其他数组。秩揭示 QR 中接近简并的选点可能受浮点与
并行布局影响，下游比较应使用已验证一致的 Wannier 基底。

### 2. 固定轨道相互作用工作流

通常的流程为：电子结构计算和匹配的 WAVECAR → 紧凑 SCDM 初始 Wannier 接口 →
外部 SAWF 构造最终对称适配轨道并输出 WANPROJ → VASP `ALGO=2E4WA` 产生
`VIJKL`、`VRijkl`。SAWF 和 WANPROJ 转换属于独立工具；本补丁读取它们产生的
WANPROJ。张量数据列、索引约定和 q→R 相位保留原有格式，输出头部记录所选模型。
原始输出定义见 [VASP VRijkl 文档](https://vasp.at/wiki/VRijkl)。

省略全部 `WINT_*` 标签时，`HFRCUT=-2` 保留此前的 bare-WS 计算路径。显式设置
`WINT_KERNEL` 时，各模型统一投影完整的固定 PAW 对密度：平滑平面波部分，加上
各原子的全电子减赝部分，包含平滑项与重建项、不同原子重建项之间的交叉贡献。
这里完整 PAW 指给定冻结芯和有限分波 PAW 表示，无需另做固体全电子计算。
显式 BARE 与省略标签的 bare 使用同一裸库仑物理核，但采用不同数值求值路径，
需分别检查收敛性。更换模型可复用已有、规范匹配的 WAVECAR 与 WANPROJ；新接口
独立于原生 hybrid/Thomas–Fermi 和 cRPA 屏蔽分支。

### 3. 边界表示

`HFRCUT=-2` 选择 WS 处理，`WINT_KERNEL` 选择物理相互作用，`WINT_BOUNDARY`
选择有限区域表示。

| 取值 | 积分区域 | 支持模型 |
|---|---|---|
| `BVK_WS`，默认 | k 网格对应 BvK 超晶格的完整三维实空间 WS 区域 | BARE、YUKAWA |
| `SLAB_WS` | 面内 BvK WS 窗口，加上每个密度各自的开放法向积分区间 | 五种模型均支持 |

`BVK_WS` 直接采用 Sundararaman–Arias 的方法 [4]，同样适用于二维等部分周期
体系。WS 区域按笛卡尔度量构造，可处理斜晶格。

`SLAB_WS` 是本项目对“面内 WS 窗口＋法向 Fourier–Galerkin 矩阵”这一具体组合的
命名。源和目标的高度分别在 `[z0−Lz/2,z0+Lz/2]` 中积分，`Lz` 是 DFT 晶胞法向
长度，周期密度在该区间展开。有限厚度算符在固定面内动量下耦合法向 Fourier
分量。它采用低维库仑截断和混合空间 Green 函数的理论基础 [5,8]，具体形式与
文献 [5] 对两点相对高度差进行截断的方法不同。

`WINT_Z0` 的单位是 Å，表示沿所选法向的笛卡尔坐标。对 slab BARE/YUKAWA，它
定位数值积分区间；对 RK_SHEET，它同时定位极化片；对 DUAL_GATE，它同时定位
两块金属板的中面。积分边界应落在真空中，并完整包含材料及密度尾部。原子层
几何中心可作为初始参考，但实际极化中心和器件环境仍需物理建模。RK2D 已沿
法向积分密度，因此不使用此参数。

### 4. 物理模型和来源

沿用英文公式中的 `C=e²/(4π epsilon_0)`、`epsilon=WINT_EPS`，`q` 为面内波矢
模长，`rho` 为面内距离，`x=z−z0`。

| 模型 | 定义及物理含义 |
|---|---|
| BARE | `C/r`，有限厚度电荷分布间的真空裸库仑相互作用 |
| YUKAWA | `C exp(−kappa r)/(epsilon r)`，均匀介质中的指数屏蔽，长度为 `1/kappa` |
| RK_SHEET | 一张位于 `z0` 的可极化绝缘片，保留 Wannier 电荷的法向位置与厚度 |
| RK2D | 对法向积分后的密度使用严格二维 RK 核 |
| DUAL_GATE | 两个对称接地金属板之间的有限厚度电荷相互作用 |

Yukawa 函数形式的历史来源见 [6]；本项目将其作为现象学屏蔽模型，也可由
`(-nabla²+kappa²)K=4πC delta(r)/epsilon` 定义。

RK_SHEET 的完整 Green 函数见英文公式，由理想极化片边界条件独立推导。
当两点均在片上时，它化为文献 [7] 的二维 RK 形式；本实现保留任意源、目标高度。
RK2D 使用

$$
G_{2D}(q)=\frac{2\pi C}{\epsilon q(1+r_0q)},
$$

先对每个 PAW 对密度沿法向积分；实空间为 Struve 与 Bessel 函数之差，近距离
呈对数形式，远距离保留 `C/(epsilon rho)` 尾部 [7]。`r0` 按上述 Fourier 约定定义。

DUAL_GATE 采用两接地平板间的有限高度 Dirichlet Green 函数 [8]，片上极限为
`2πC tanh(qd)/(epsilon q)`，与文献 [9] 讨论的二维门控形式相关。金属板位于
`z0±d`；`WINT_GATE_D` 表示到每侧板的距离，总板间距是 `2d`。当前积分域要求
`d>Lz/2`，这是实现的区域限制。这里计算两电子相互作用；单电子镜像势需在
单粒子 Hamiltonian 中另行考虑。

在相同区域、`epsilon=1` 时，Yukawa 的 `kappa=0`、RK_SHEET 的 `r0=0` 和
双栅极的 `d→∞` 均回到有限厚度 BARE；一般 epsilon 下为 BARE/epsilon。
RK2D 的 `r0=0` 则是作用于高度积分密度的平面 `C/(epsilon rho)`。
有限 WS 区域的零模与非零分量由同一截断算符构造。

层状求值器使用 Gaussian 表示分离面内与法向坐标：平滑面内部分采用 WS FFT，
短程面内因子解析计算，法向矩阵元使用一维矩积分。同一重建 PAW 密度进入原生
张量收缩和输出。

### 5. 输入及适用范围

相互作用投影路径支持 CPU NCL、`ALGO=2E4WA`、`HFRCUT=-2`、`ISPIN=1` 和
`NCORE=1`，包括非共线 SOC 旋量。层状模型要求法向只有一个 k 点，第三晶格矢量
垂直于面内晶格；面内可为斜晶格。

| 标签 | 含义与默认值 |
|---|---|
| `WINT_KERNEL` | BARE、YUKAWA、RK_SHEET、RK2D、DUAL_GATE；物理默认 BARE，省略全部 WINT 标签保留旧求值路径 |
| `WINT_BOUNDARY` | 默认 BVK_WS；RK_SHEET、RK2D、DUAL_GATE 必须显式指定 SLAB_WS |
| `WINT_KAPPA` | YUKAWA 必需，Å⁻¹，有限非负 |
| `WINT_EPS` | 屏蔽模型的相对介电常数，默认 1，有限正值；显式 BARE 不接受此标签 |
| `WINT_R0` | RK_SHEET、RK2D 必需，Å，有限非负 |
| `WINT_GATE_D` | DUAL_GATE 必需，Å，到每侧金属板的距离，需大于 Lz/2 |
| `WINT_Z0` | 有限厚度 SLAB_WS 必需，Å；RK2D 省略 |
| `WINT_NORMAL` | 可选笛卡尔法向，默认第三晶格矢量方向；指定方向须与它平行 |
| `WS_NGRID` | 默认自动确定核 FFT 网格；三分量，SLAB_WS 只使用前两项 |
| `WS_TOL` | 核构造与内部积分容差，默认 `1e-8` |

例如材料参考面位于笛卡尔 `z=0` 时，在已有固定轨道 2E4WA 输入中添加：

```text
WINT_KERNEL = DUAL_GATE
WINT_BOUNDARY = SLAB_WS
WINT_Z0 = 0.0
WINT_GATE_D = 100.0
```

该示例表示每侧金属板距参考面 100 Å，即 10 nm；实际器件需使用相应参考面和介电
参数。不支持的计算模式、矛盾参数组合和模型无关参数会明确停止。

### 6. 收敛和验证范围

核网格 `WS_NGRID`、重建密度响应基组 `VCUTOFF` 与 `NGX/NGY/NGZ`、BvK
范围与 k 网格（SLAB_WS 中为面内方向）是不同的收敛方向。`WS_TOL` 控制核构造和内部积分，不直接给出
最终张量误差。仅增加响应截断可能触及原密度 FFT 网格的上限，应结合 OUTCAR
实际参数、诊断和完整张量变化来判断。

若增大密度 FFT 网格后原 CHGCAR 不兼容，已测试的 `ICHARG=12` 用法针对固定
WAVECAR 的相互作用投影；下游仍使用已建立的 SAWF 单粒子 Hamiltonian。这项
检验并不覆盖附带单粒子或 FCIDUMP 量的等价性。

保留的实现验证包括：GNU Fortran 与 Intel ifx 各 32 项数值测试、76 项源码/输入/
审计测试、21 次最终程序原生计算，以及两个源码版本和 Git/非 Git 目标上的 12 种
升级回退流程。检验内容覆盖独立物理参考、复数对密度、零屏蔽极限、介电缩放、
MPI/KPAR 和旧输出回归。具体材料仍需检查响应截断、网格、空间范围及实际 CPU/
内存。PAW 模块与层状模块各自将新增自有工作区限制为 512 MiB；原生调用方的
波函数和密度数组另行计入总进程内存。fRG 模型的屏蔽参数还应与后续低能空间的
屏蔽处理保持一致。

## References / 参考文献

1. A. Damle, L. Lin and L. Ying, “Compressed Representation of Kohn–Sham Orbitals via Selected Columns of the Density Matrix,” *J. Chem. Theory Comput.* **11**, 1463–1469 (2015). [DOI](https://doi.org/10.1021/ct500985f).
2. A. Damle, L. Lin and L. Ying, “SCDM-k: Localized orbitals for solids via selected columns of the density matrix,” *J. Comput. Phys.* **334**, 1–15 (2017). [DOI](https://doi.org/10.1016/j.jcp.2016.12.053), [arXiv](https://arxiv.org/abs/1507.03354).
3. VASP's [SCDM interface](https://vasp.at/wiki/LSCDM) cites A. Damle and L. Lin, “Disentanglement via Entanglement: A Unified Method for Wannier Localization,” *Multiscale Model. Simul.* **16**, 1392–1410 (2018). [DOI](https://doi.org/10.1137/17M1129696), [arXiv](https://arxiv.org/abs/1703.06958).
4. R. Sundararaman and T. A. Arias, “Ideal regularization of the Coulomb singularity in exact exchange by Wigner–Seitz truncated interactions,” *Phys. Rev. B* **87**, 165122 (2013), Sec. II and Appendix A. [DOI](https://doi.org/10.1103/PhysRevB.87.165122), [arXiv](https://arxiv.org/abs/1302.6204).
5. S. Ismail-Beigi, “Truncation of periodic image interactions for confined systems,” *Phys. Rev. B* **73**, 233103 (2006). [DOI](https://doi.org/10.1103/PhysRevB.73.233103), [arXiv](https://arxiv.org/abs/cond-mat/0603448).
6. H. Yukawa, “On the Interaction of Elementary Particles. I,” *Proc. Phys.-Math. Soc. Japan* **17**, 48–57 (1935). [DOI](https://doi.org/10.11429/ppmsj1919.17.0_48). Historical source of the Yukawa form; the electronic screening model used here is specified by its displayed equation.
7. P. Cudazzo, I. V. Tokatly and A. Rubio, “Dielectric screening in two-dimensional insulators: Implications for excitonic and impurity states in graphane,” *Phys. Rev. B* **84**, 085406 (2011), Eqs. 5–7. [DOI](https://doi.org/10.1103/PhysRevB.84.085406), [arXiv](https://arxiv.org/abs/1104.3346).
8. K. Takae and A. Onuki, “Applying electric field to charged and polar particles between metallic plates: Extension of the Ewald method,” *J. Chem. Phys.* **139**, 124108 (2013), Eqs. 2.17–2.21. [DOI](https://doi.org/10.1063/1.4821085), [arXiv](https://arxiv.org/abs/1305.3503).
9. R. E. Throckmorton and O. Vafek, “Fermions on bilayer graphene: symmetry breaking for B=0 and nu=0,” *Phys. Rev. B* **86**, 115447 (2012). [DOI](https://doi.org/10.1103/PhysRevB.86.115447), [arXiv](https://arxiv.org/abs/1111.2076). The present implementation defines the gate geometry explicitly through the physical distance `d` to each plate.
