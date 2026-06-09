#!/usr/bin/env bash
# source_env.sh - Inject .env auto-sourcing into shell profile
# Called by postCreateCommand to survive container rebuilds.
set -euo pipefail

MARKER="_source_project_env"
TARGET="${HOME}/.bashrc"

# Idempotent: skip if already present
if grep -q "$MARKER" "$TARGET" 2>/dev/null; then
    exit 0
fi

cat >> "$TARGET" << 'EOF'

# Auto-source project .env files (API keys for Gemini, PAL, etc.)
_source_project_env() {
    local project_dir=""
    if [ -d "/workspaces" ]; then
        for d in /workspaces/*/; do
            [ -f "${d}.devcontainer/.env" ] && project_dir="${d%/}" && break
        done
    fi
    [ -z "$project_dir" ] && project_dir="${PWD}"
    [ -f "${project_dir}/.devcontainer/.env" ] && { set -a; source "${project_dir}/.devcontainer/.env"; set +a; }
    [ -f "${project_dir}/.env" ] && { set -a; source "${project_dir}/.env"; set +a; }
}
_source_project_env
EOF
