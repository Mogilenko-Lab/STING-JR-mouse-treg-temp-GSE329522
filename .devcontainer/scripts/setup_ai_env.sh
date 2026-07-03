#!/usr/bin/env bash
# setup_ai_env.sh — poststart: install AI coding CLIs + seed Claude power-user config.
#
# Shipped by scbio-docker's init-container.sh into every project's
# .devcontainer/scripts/ and invoked from devcontainer.json's postStartCommand.
# Runs on every container START. The container home dir is NOT persisted across
# rebuilds, so this re-establishes the AI toolchain + settings each time.
#
# Design contract:
#   - IDEMPOTENT: skips any CLI already on PATH; safe to run every start.
#   - NETWORK-TOLERANT: a failed/unavailable installer only WARNS; it must never
#     block the container from starting (hence `set -uo pipefail`, not `-e`).
#   - Image prereqs (baked into scdock-r-dev): node/npm, curl, jq. AI tooling is
#     intentionally NOT baked into the image — it lives here, at runtime.

set -uo pipefail

log(){ printf '[setup_ai_env] %s\n' "$*"; }

# Installers drop binaries into these; put them on PATH for this run + probes
# (claude/codex/agy → ~/.local/bin, opencode → ~/.opencode/bin, npm-global fallback).
export PATH="$HOME/.local/bin:$HOME/.opencode/bin:$HOME/.npm-global/bin:$PATH"

######################################################################
# 1. AI coding CLIs  (install once per fresh container; skip if present)
######################################################################
install_if_missing(){
  # $1 = probe command   $2 = human name   $3 = installer command line
  local probe="$1" name="$2" cmd="$3"
  if command -v "$probe" >/dev/null 2>&1; then
    log "$name present ($(command -v "$probe")) — skip"
    return 0
  fi
  log "installing $name ..."
  if eval "$cmd"; then log "$name installed"; else log "WARN: $name install failed (continuing)"; fi
}

install_if_missing claude   "Claude Code"  'curl -fsSL https://claude.ai/install.sh | bash'
install_if_missing codex    "OpenAI Codex" 'curl -fsSL https://chatgpt.com/codex/install.sh | sh'
install_if_missing opencode "opencode"     'curl -fsSL https://opencode.ai/install | bash'

# Pi agent — the official pi.dev installer is interactive AND wants Node >=22.19;
# the npm package installs cleanly headless on the image's Node 20 into a
# user-writable prefix (~/.local, already on PATH). configure_pi_models.sh writes
# the local-model roster separately.
if command -v pi >/dev/null 2>&1; then
  log "Pi agent present ($(command -v pi)) — skip"
else
  log "installing Pi agent (npm) ..."
  npm config set prefix "$HOME/.local" >/dev/null 2>&1
  if npm install -g @mariozechner/pi-coding-agent >/dev/null 2>&1; then
    log "Pi agent installed ($("$HOME/.local/bin/pi" --version 2>/dev/null))"
  else
    log "WARN: Pi agent install failed (continuing)"
  fi
fi

# Antigravity CLI ships as `antigravity` (short alias `agy`); probe both.
if command -v antigravity >/dev/null 2>&1 || command -v agy >/dev/null 2>&1; then
  log "Antigravity present — skip"
else
  log "installing Antigravity CLI ..."
  curl -fsSL https://antigravity.google/cli/install.sh | bash \
    || log "WARN: Antigravity install failed (continuing)"
fi

######################################################################
# 2. Claude power-user settings + status line   (USER scope: ~/.claude)
######################################################################
# User scope applies to EVERY project opened in this container. Home is not
# persisted across rebuilds, so we re-seed each start. Any pre-existing keys are
# preserved via deep-merge (our values win on conflict), so a hand-edit made
# during a running session survives a same-container restart.
claude_dir="$HOME/.claude"
mkdir -p "$claude_dir"

# 2a. Status line — fields taken from the REAL statusline payload
#     (.model.display_name, .context_window.*), never the invented ones.
cat > "$claude_dir/statusline.sh" <<'SL'
#!/bin/bash
input=$(cat)
model=$(echo "$input" | jq -r '.model.display_name')
used=$(echo "$input" | jq -r '.context_window.total_input_tokens // 0')
free=$(echo "$input" | jq -r '.context_window.remaining_percentage // 100' | cut -d. -f1)
bar_size=10
filled=$(( (100 - free) * bar_size / 100 ))
empty=$(( bar_size - filled ))
bar=$(printf "%${filled}s" | tr ' ' '#')$(printf "%${empty}s" | tr ' ' '.')
echo "[$model] [$bar] ${used} tokens (${free}% free)"
SL
chmod +x "$claude_dir/statusline.sh"

# 2b. settings.json — schema-valid power-user defaults. Telemetry is left ON on
#     purpose: disabling it trips the feature-flag layer that gates agent teams /
#     1M context, so we do NOT set DISABLE_TELEMETRY here.
desired="$(cat <<'JSON'
{
  "$schema": "https://json.schemastore.org/claude-code-settings.json",
  "editorMode": "vim",
  "alwaysThinkingEnabled": true,
  "effortLevel": "xhigh",
  "showThinkingSummaries": true,
  "autoMemoryEnabled": false,
  "teammateMode": "auto",
  "statusLine": { "type": "command", "command": "bash \"$HOME/.claude/statusline.sh\"", "padding": 2 }
}
JSON
)"

settings="$claude_dir/settings.json"
if command -v jq >/dev/null 2>&1 && [ -s "$settings" ]; then
  if merged="$(jq -s '.[0] * .[1]' "$settings" <(printf '%s' "$desired") 2>/dev/null)"; then
    printf '%s\n' "$merged" > "$settings"
  else
    printf '%s\n' "$desired" > "$settings"   # existing file unparseable → overwrite
  fi
else
  printf '%s\n' "$desired" > "$settings"
fi
log "Claude settings seeded → $settings"

log "done."
