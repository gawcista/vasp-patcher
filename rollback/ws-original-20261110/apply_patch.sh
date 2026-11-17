#!/usr/bin/env bash
set -u

# Private development source provenance; deployment trusts exact hashes, not Git ancestry.
# Regenerated with the unified source payload; no version-text bypass exists.
BASE_COMMIT="71970129be8fbdd8a7fe12aba8bcaea0395f8131"
PATCHED_HEAD="1aec886c040dfd0bccd276aac8657dc02a883886"
PATCH_SHA256="387d745897017f022665459dd47c6d061aed4eae7aa55edcad35f711572e2803"
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
PATCH_FILE="$SCRIPT_DIR/vasp-6.6.1-scdm-compact-bands.patch"

PATHS=(
  "src/mlwf.F"
  "src/fock_glb.F"
  "src/fock.F"
  "src/dmft.F"
  "src/.objects"
  "src/ws_coulomb.F"
)
PRE_HASHES_661=(
  "7de0801b8cb93bb8813b9a8753e07f890308d317abf6851ac90060b3a2cf8e42"
  "380a98f539a9935e87c626d7688324460e4f20d62f34bf154e9aa2021b1dd198"
  "12d68e397d198af2bb6235fe7be303499af0d577bb517234b70acbba2395ef1c"
  "f88b7ed00729f87c85b641871a481d06faf106657e0060c2bc1d4e7b908c03a6"
  "d683fb52ffcc6db69726f23bb6bc79390b9088d72691a01ae040ff8ed0fca0b8"
  "MISSING"
)
POST_HASHES_661=(
  "48a6d8351f3987d286221c88f268c925f8d5ed10ecbfaa855d1f3905a9df59ca"
  "9d720e38f325a36a7134755bcc39873066c6add6886679dfe7db938bcf229bd5"
  "9f9c837662a0d14d293e27f190dcb0c2aa6381f19d8f1af9cd59109a0fd46777"
  "0d8b3af23138c9e9a4eb07926ce060cbdcc47d4adadc8ba710f9bb732e3916dc"
  "e04af0330cda2c7b14d38b8b4d133ead7298b0c265f837553776cda009601a16"
  "d18f9acc8554a0d158a46a878172d32139d9bec64972efa7472c3a91573582ff"
)
PRE_HASHES_660=(
  "7de0801b8cb93bb8813b9a8753e07f890308d317abf6851ac90060b3a2cf8e42"
  "380a98f539a9935e87c626d7688324460e4f20d62f34bf154e9aa2021b1dd198"
  "f2c34313379e6b21f93b3777696c7e74f4f01cd1ef7d4cb7414624a4c584dc5d"
  "f88b7ed00729f87c85b641871a481d06faf106657e0060c2bc1d4e7b908c03a6"
  "d683fb52ffcc6db69726f23bb6bc79390b9088d72691a01ae040ff8ed0fca0b8"
  "MISSING"
)
POST_HASHES_660=(
  "48a6d8351f3987d286221c88f268c925f8d5ed10ecbfaa855d1f3905a9df59ca"
  "9d720e38f325a36a7134755bcc39873066c6add6886679dfe7db938bcf229bd5"
  "7f74691fac65350db6e3252533346f0fcadfd6afbdb37a05329254cdc2fb7c5c"
  "0d8b3af23138c9e9a4eb07926ce060cbdcc47d4adadc8ba710f9bb732e3916dc"
  "e04af0330cda2c7b14d38b8b4d133ead7298b0c265f837553776cda009601a16"
  "d18f9acc8554a0d158a46a878172d32139d9bec64972efa7472c3a91573582ff"
)

usage() {
  echo "Usage: $0 --status|--check|--apply|--reverse TARGET_ROOT" >&2
  exit 2
}
die() { echo "ERROR: $*" >&2; exit 2; }

hash_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

verify_package() {
  [ -f "$PATCH_FILE" ] && [ ! -L "$PATCH_FILE" ] || die "patch payload is missing or symlinked"
  local actual
  actual=$(hash_file "$PATCH_FILE") || die "cannot hash patch payload"
  [ "$actual" = "$PATCH_SHA256" ] || die "patch payload hash mismatch"
}

select_profile() {
  case "$1" in
    6.6.1) PRE_HASHES=("${PRE_HASHES_661[@]}"); POST_HASHES=("${POST_HASHES_661[@]}") ;;
    6.6.0) PRE_HASHES=("${PRE_HASHES_660[@]}"); POST_HASHES=("${POST_HASHES_660[@]}") ;;
  esac
}

tree_state() {
  local path actual index profile pre post scdm state found=""
  ACTUAL_HASHES=()
  [ ! -L "$ROOT/src" ] || die "refusing symlinked source directory"
  for path in "${PATHS[@]}"; do
    [ ! -L "$ROOT/$path" ] || die "refusing symlinked source: $path"
    if [ -f "$ROOT/$path" ]; then
      actual=$(hash_file "$ROOT/$path") || die "cannot hash $path"
    elif [ -e "$ROOT/$path" ]; then
      actual="NOT_A_FILE"
    else
      actual="MISSING"
    fi
    ACTUAL_HASHES+=("$actual")
  done
  PROFILE="unknown"; TREE_STATE="modified/unknown"
  for profile in 6.6.1 6.6.0; do
    select_profile "$profile"
    [ "${#PRE_HASHES[@]}" -eq "${#PATHS[@]}" ] &&
      [ "${#POST_HASHES[@]}" -eq "${#PATHS[@]}" ] || die "incomplete embedded profile"
    pre=0; post=0; scdm=0
    for index in "${!PATHS[@]}"; do
      actual=${ACTUAL_HASHES[$index]}
      [ "$actual" != "${PRE_HASHES[$index]}" ] || pre=$((pre + 1))
      [ "$actual" != "${POST_HASHES[$index]}" ] || post=$((post + 1))
      if [ "${PATHS[$index]}" = "src/mlwf.F" ]; then
        [ "$actual" != "${POST_HASHES[$index]}" ] || scdm=$((scdm + 1))
      else
        [ "$actual" != "${PRE_HASHES[$index]}" ] || scdm=$((scdm + 1))
      fi
    done
    state=""
    if [ "$pre" -eq "${#PATHS[@]}" ]; then state="clean-base"
    elif [ "$post" -eq "${#PATHS[@]}" ]; then state="applied"
    elif [ "$scdm" -eq "${#PATHS[@]}" ]; then state="scdm-only"
    fi
    if [ -n "$state" ]; then
      [ -z "$found" ] || die "ambiguous embedded source profiles"
      found=$profile; PROFILE=$profile; TREE_STATE=$state
    fi
  done
  if [ "$PROFILE" != "unknown" ]; then select_profile "$PROFILE"; fi
}

print_status() {
  local index state
  tree_state
  echo "VASP unified compact-SCDM + WS source package"
  echo "  profile: $PROFILE"
  echo "  root: $ROOT"
  for index in "${!PATHS[@]}"; do
    state="unknown"
    if [ "$PROFILE" != "unknown" ]; then
      if [ "${ACTUAL_HASHES[$index]}" = "${PRE_HASHES[$index]}" ]; then state="pre"
      elif [ "${ACTUAL_HASHES[$index]}" = "${POST_HASHES[$index]}" ]; then state="post"; fi
    fi
    printf "  %-8s %s\n" "$state" "${PATHS[$index]}"
  done
  echo "status: $TREE_STATE"
  if [ "$TREE_STATE" = "applied" ]; then echo "configuration: combined SCDM+WS"; fi
}

is_target_git_root() {
  command -v git >/dev/null 2>&1 || return 1
  local top
  top=$(git -C "$ROOT" rev-parse --show-toplevel 2>/dev/null) || return 1
  [ "$(CDPATH= cd -- "$top" && pwd -P)" = "$ROOT" ]
}

prepare_patch() {
  local omit=$1 expected
  expected=$(printf '%s\n' "${PATHS[@]}")
  SELECTED_PATCH=$(mktemp "${TMPDIR:-/tmp}/vasp-unified-patch.XXXXXX") || die "cannot create temporary patch"
  # Keep complete diff blocks; validate the exact path set before using any of it.
  awk -v paths="$expected" -v omit="$omit" '
    BEGIN { n=split(paths,p,"\n"); for(i=1;i<=n;i++) allowed[p[i]]=1; emit=0; bad=0 }
    /^diff --git / {
      path=substr($3,3)
      if(NF!=4 || $3!="a/"path || $4!="b/"path || !allowed[path] || seen[path]++) {
        print "Invalid or duplicate unified-patch path" > "/dev/stderr"; bad=1; exit 2
      }
      emit=(path!=omit)
    }
    { if(emit) print }
    END {
      for(i=1;i<=n;i++) if(!seen[p[i]]) { print "Missing unified-patch path: "p[i] > "/dev/stderr"; bad=1 }
      if(bad) exit 2
    }
  ' "$PATCH_FILE" > "$SELECTED_PATCH" || die "invalid unified payload"
}

run_patch() {
  local direction=$1 dry=$2
  if is_target_git_root; then
    local args=(apply)
    [ "$dry" != "yes" ] || args+=(--check)
    [ "$direction" != "reverse" ] || args+=(--reverse)
    git -C "$ROOT" "${args[@]}" "$SELECTED_PATCH"
  else
    command -v patch >/dev/null 2>&1 || die "non-Git mode requires patch"
    local args=(--batch --fuzz=0 --no-backup-if-mismatch --reject-file=- -p1)
    [ "$dry" != "yes" ] || args+=(--dry-run)
    if [ "$direction" = "reverse" ]; then args+=(--reverse); else args+=(--forward); fi
    (cd "$ROOT" && patch "${args[@]}" < "$SELECTED_PATCH")
  fi
}

require_known() {
  [ "$PROFILE" != "unknown" ] || die "refusing mixed, missing or modified source files; run --status"
}

make_backup() {
  local label=$1 path destination stamp
  stamp=$(date -u +%Y%m%dT%H%M%SZ)
  BACKUP_DIR="$ROOT/.scdm-compact-bands-backups/${label}-${stamp}-$$"
  mkdir -p "$BACKUP_DIR/files" || die "cannot create backup directory"
  : > "$BACKUP_DIR/missing.list"
  for path in "${PATHS[@]}"; do
    destination="$BACKUP_DIR/files/$path"
    if [ -f "$ROOT/$path" ]; then
      mkdir -p "$(dirname -- "$destination")" || die "cannot create backup path"
      cp -p "$ROOT/$path" "$destination" || die "cannot back up $path"
    else
      printf '%s\n' "$path" >> "$BACKUP_DIR/missing.list"
    fi
  done
  printf '%s\n' "$TREE_STATE" > "$BACKUP_DIR/source-state"
  printf '%s\n' "$PROFILE" > "$BACKUP_DIR/source-profile"
  echo "source backup: $BACKUP_DIR"
}

restore_backup() {
  local path source
  for path in "${PATHS[@]}"; do
    source="$BACKUP_DIR/files/$path"
    if [ -f "$source" ]; then
      cp -p "$source" "$ROOT/$path" || die "cannot restore $path from $BACKUP_DIR"
    else
      rm -f -- "$ROOT/$path" || die "cannot remove newly added $path"
    fi
  done
}

change_source() {
  local direction=$1 wanted old_profile=$PROFILE
  if [ "$direction" = "forward" ]; then wanted="applied"; else wanted="clean-base"; fi
  run_patch "$direction" yes || die "patch dry-run failed; no source files changed"
  make_backup "before-$direction"
  if ! run_patch "$direction" no; then
    restore_backup
    die "patch failed; source files restored from $BACKUP_DIR"
  fi
  tree_state
  if [ "$PROFILE" != "$old_profile" ] || [ "$TREE_STATE" != "$wanted" ]; then
    restore_backup
    die "source hash verification failed; source files restored from $BACKUP_DIR"
  fi
  echo "Source state $TREE_STATE verified for VASP $PROFILE."
}

[ "$#" -eq 2 ] || usage
MODE=$1
case "$MODE" in --status|--check|--apply|--reverse) ;; *) usage ;; esac
[ -d "$2" ] || die "target root is not a directory: $2"
ROOT=$(CDPATH= cd -- "$2" && pwd -P)
command -v awk >/dev/null 2>&1 || die "awk is required"
command -v sha256sum >/dev/null 2>&1 || command -v shasum >/dev/null 2>&1 || die "SHA256 utility required"
SELECTED_PATCH=""
trap '[ -z "$SELECTED_PATCH" ] || rm -f -- "$SELECTED_PATCH"' EXIT
verify_package
tree_state
if [ "$MODE" = "--status" ]; then print_status; exit 0; fi
require_known
if [ "$MODE" = "--apply" ] && [ "$TREE_STATE" = "applied" ]; then
  echo "Patch is already applied; no files were changed."; exit 0
fi
if [ "$MODE" = "--reverse" ]; then
  [ "$TREE_STATE" = "applied" ] || die "reverse requires the complete combined SCDM+WS state"
  prepare_patch ""
  change_source reverse
else
  [ "$TREE_STATE" != "applied" ] || die "package is already applied; no files were changed"
  if [ "$TREE_STATE" = "scdm-only" ]; then prepare_patch "src/mlwf.F"; else prepare_patch ""; fi
  if [ "$MODE" = "--check" ]; then
    run_patch forward yes || die "patch dry-run failed"
    echo "Preimage hashes and unified patch dry-run verified for VASP $PROFILE ($TREE_STATE)."
  else
    change_source forward
    echo "Rebuild VASP using the site's normal private build procedure."
  fi
fi
