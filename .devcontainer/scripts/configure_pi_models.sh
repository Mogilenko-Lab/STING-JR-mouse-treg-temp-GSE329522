#!/usr/bin/env bash
# Configure Pi Coding Agent local-model roster.
# Writes ~/.pi/agent/models.json on every container start so the config
# survives rebuilds (home dir is not persisted across container recreations).
# Pi installation is separate — this only writes the config.

set -euo pipefail

OLLAMA_HOST="${OLLAMA_HOST:-http://172.17.0.1:11434}"

mkdir -p "$HOME/.pi/agent"

cat > "$HOME/.pi/agent/models.json" << EOF
{
  "providers": {
    "ollama": {
      "baseUrl": "${OLLAMA_HOST}/v1",
      "api": "openai-completions",
      "apiKey": "ollama",
      "compat": {
        "supportsDeveloperRole": false,
        "supportsReasoningEffort": false
      },
      "models": [
        { "id": "batiai/qwen3.6-27b:iq4",                                    "name": "Qwen 3.6 27B Dense (Primary Agent)",  "reasoning": false, "contextWindow": 131072 },
        { "id": "qwen3-coder:30b",                                            "name": "Qwen 3 Coder 30B (Code)",             "reasoning": false, "contextWindow": 131072 },
        { "id": "qwen3.6:35b",                                               "name": "Qwen 3.6 35B MoE (Fast Agent)",       "reasoning": false, "contextWindow": 131072 },
        { "id": "gemma4:31b",                                                "name": "Gemma 4 31B (Multimodal Agent)",      "reasoning": false, "contextWindow": 131072 },
        { "id": "Sub01/FuseO1-DeepSeekR1-QwQ-SkyT1-Flash-32B-Preview:q4_K_M", "name": "FuseO1 32B (Reasoning)",           "reasoning": true,  "contextWindow": 32768 },
        { "id": "hengwen/DeepSeek-R1-Distill-Qwen-32B:q4_K_M",              "name": "DeepSeek-R1 Qwen 32B (Reasoning)",   "reasoning": true,  "contextWindow": 32768 }
      ]
    }
  }
}
EOF

echo "Pi Agent: models configured (Ollama @ ${OLLAMA_HOST})"
