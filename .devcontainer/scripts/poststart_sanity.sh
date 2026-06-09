#!/usr/bin/env bash
set -euo pipefail

ok()  { echo "OK"; }
fail(){ echo "NOT OK"; }

sep(){ echo "==== $* ===="; }

sep "Devcontainer sanity checks"

printf "Python: "
if python -c "import sys; print(sys.version.split()[0])" >/dev/null 2>&1; then ok; else fail; fi

printf "R:      "
if R -q --vanilla -e "invisible(TRUE)" >/dev/null 2>&1; then ok; else fail; fi

printf "httpgd: "
if R -q --vanilla -e "suppressMessages(library(httpgd)); cat(as.character(packageVersion('httpgd')), '\n')" >/dev/null 2>&1; then ok; else 
  fail; 
  echo "Hint: install with CRAN or GitHub fallback:"; 
  echo "  R -q --vanilla -e 'install.packages(\"httpgd\")'";
  echo "  R -q --vanilla -e 'if(!requireNamespace(\"httpgd\",quietly=TRUE)){if(!requireNamespace(\"remotes\",quietly=TRUE)) install.packages(\"remotes\"); remotes::install_github(\"nx10/httpgd\")}'";
fi

printf "scanpy import: "
if python -c "import scanpy" >/dev/null 2>&1; then ok; else fail; fi

printf "default venv: "
if which python | grep -q "/opt/venvs/base/bin/python"; then echo "/opt/venvs/base" OK; else fail; fi

# radian and wrappers
printf "radian: "
if command -v radian >/dev/null 2>&1; then ok; else fail; fi

printf "r-base wrapper: "
if command -v r-base >/dev/null 2>&1; then ok; else echo "SKIP (wrapper not present)"; fi

printf "r-archr wrapper: "
if command -v r-archr >/dev/null 2>&1; then 
  if [ "${IMAGE_FLAVOR:-dev-core}" = "dev-archr" ]; then echo OK; else echo "OK (note: ArchR not installed in dev-core)"; fi
else 
  echo "SKIP (wrapper not present)"; 
fi

# Python env switchers
printf "usepy function: "
if grep -q "usepy(){" /etc/bash.bashrc; then ok; else fail; fi

printf "py-* wrappers: "
missing=0; for w in py-base py-squid py-atac py-comms; do command -v "$w" >/dev/null 2>&1 || missing=1; done; if [ $missing -eq 0 ]; then ok; else fail; fi

printf "py-squid path: "
if py-squid which python 2>/dev/null | grep -q "/opt/venvs/squid/bin/python"; then ok; else fail; fi

printf "py-atac path: "
if py-atac which python 2>/dev/null | grep -q "/opt/venvs/atac/bin/python"; then ok; else fail; fi

printf "py-comms path: "
if py-comms which python 2>/dev/null | grep -q "/opt/venvs/comms/bin/python"; then ok; else fail; fi

# R library write test (user library should be writable)
printf "R lib writable: "
if R -q --vanilla -e "lib <- .libPaths()[1]; ok <- file.access(lib, 2) == 0; if (!ok) stop('not writable: ', lib)" >/dev/null 2>&1; then ok; else 
  fail;
  echo "Hint: .libPaths()[1] must be writable by container user."
  echo "id: $(id -u):$(id -g) group '$(id -gn)'  pwd: $(pwd)"
fi

printf "R lib first: "
R_FIRST_LIB=$(R -q -e 'cat(.libPaths()[1])' 2>/dev/null || true)
if [ -n "$R_FIRST_LIB" ]; then echo "$R_FIRST_LIB"; else echo "UNKNOWN"; fi

# Workspace write test
printf "Workspace write: "
tmpdir="._write_sanity_$$"
if mkdir -p "$tmpdir" 2>/dev/null && echo ok > "$tmpdir/.touch" 2>/dev/null; then 
  ok; 
  rm -rf "$tmpdir" 2>/dev/null || true;
else 
  fail; 
  echo "Check mount ownership/ACLs for: $(pwd)";
  echo "ls -ld . =>"; ls -ld . || true;
  echo "id: $(id -u):$(id -g) group '$(id -gn)'";
fi

sep "Tips"
echo "usepy base|squid|atac|comms  # switch Python env"
if command -v r-base >/dev/null 2>&1; then
  echo "r-base / r-archr            # regular vs ArchR R sessions"
else
  echo "Run 'radian' for R; ArchR variant available in ArchR image"
fi

echo "Current image flavor: ${IMAGE_FLAVOR:-unknown}"
echo "Active Python: $(which python)"
echo "To switch Python env: usepy base|squid|atac|comms"
if [ "${IMAGE_FLAVOR:-dev-core}" = "dev-archr" ]; then
  echo "R sessions: r-base (no ArchR), r-archr (ArchR-enabled via USE_ARCHR=1)"
else
  echo "R sessions: radian only (no ArchR). For ArchR, use dev-archr service."
fi

