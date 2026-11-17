# Matching installers for older source states / 旧补丁回退工具

Use an old installer only to reverse the exact older source state it recognizes.
Each directory contains its matching installer and unified payload:

| Directory | Older package |
|---|---|
| `ws-original-20261110` | Original bare-WS package |
| `ws-stable-20261113` | Bare-WS package with updated diagnostics |
| `yukawa-20261116` | Seven-path full-PAW BARE/YUKAWA package |

From the repository root, set the actual licensed source directory and inspect
which old package reports `status: applied`:

```bash
VASP_WINT_SRC=/absolute/path/to/vasp-source
bash rollback/ws-original-20261110/apply_patch.sh --status "$VASP_WINT_SRC"
bash rollback/ws-stable-20261113/apply_patch.sh --status "$VASP_WINT_SRC"
bash rollback/yukawa-20261116/apply_patch.sh --status "$VASP_WINT_SRC"
```

Use the matching pair to reverse, then use the current root installer. For example:

```bash
bash rollback/ws-original-20261110/apply_patch.sh --reverse "$VASP_WINT_SRC"
bash apply_patch.sh --check "$VASP_WINT_SRC"
bash apply_patch.sh --apply "$VASP_WINT_SRC"
make -C "$VASP_WINT_SRC" ncl DEPS=1 -j4
```

These pairs restore recognized vendor bytes. They do not convert arbitrary
locally modified sources. The current installer handles clean vendor and
recognized SCDM-only sources directly; those do not require an old rollback.

---

每个子目录中的脚本与补丁严格配对，仅用于回退它能识别的旧源码状态。
先逐个执行 `--status`，选择显示 `status: applied` 的旧版本，再执行
`--reverse`。之后使用根目录中的新版脚本检查、应用并重新编译。

已识别的干净原版源码或仅安装 SCDM 的源码可以直接使用新版脚本。
回退会恢复原版源码；局部修改、混合版本和未知哈希需要先核对。
