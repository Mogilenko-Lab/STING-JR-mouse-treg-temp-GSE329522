#!/usr/bin/env bash
# configure_sciagent_alias.sh - Install the si -> sciagent alias in shell rc.
# CWD-relative (not a PATH symlink) so it resolves to this project's own
# vendored 01_modules/SciAgent-toolkit wherever you cd -- a fixed symlink bakes
# one checkout as the resolution target and misfires against a different
# project's toolkit copy (SciAgent-toolkit README.md "Install",
# docs/architecture.md §5/§13). Re-applied on every start since the home dir
# does not persist across container recreations.
set -euo pipefail

MARKER="_sciagent_si_alias"
TARGET="${HOME}/.bashrc"

grep -q "$MARKER" "$TARGET" 2>/dev/null && exit 0

cat >> "$TARGET" << 'EOF'

# _sciagent_si_alias: CWD-relative, safe across multiple vendored toolkit copies
alias si="./01_modules/SciAgent-toolkit/bin/sciagent"
EOF
