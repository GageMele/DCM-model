#!/bin/bash
set -euo pipefail
for name in onset_metrics.R predictive.R; do
  old_hash=$(awk -v n="$name" '$2==n {print $1}' "$ONSET_REPAIR_DIR/original_SHA256SUMS.txt")
  new_hash=$(awk -v n="$name" '$2==n {print $1}' "$ONSET_REPAIR_DIR/patched_SHA256SUMS.txt")
  actual=$(sha256sum "$EXPERIMENT_DIR/$name" | cut -d' ' -f1)
  if [[ "$actual" != "$old_hash" && "$actual" != "$new_hash" ]]; then
    echo "Unexpected scoring version: $name; no scoring files changed." >&2; exit 1
  fi
done
mkdir -p "$ONSET_REPAIR_DIR/scoring_backups"
for name in onset_metrics.R predictive.R; do
  new_hash=$(awk -v n="$name" '$2==n {print $1}' "$ONSET_REPAIR_DIR/patched_SHA256SUMS.txt")
  actual=$(sha256sum "$EXPERIMENT_DIR/$name" | cut -d' ' -f1)
  [[ "$actual" == "$new_hash" ]] && continue
  [[ -e "$ONSET_REPAIR_DIR/scoring_backups/$name" ]] || cp -- "$EXPERIMENT_DIR/$name" "$ONSET_REPAIR_DIR/scoring_backups/$name"
  cp -- "$ONSET_REPAIR_DIR/$name" "$EXPERIMENT_DIR/$name.repair_r1.tmp"
  mv -- "$EXPERIMENT_DIR/$name.repair_r1.tmp" "$EXPERIMENT_DIR/$name"
done
(cd "$EXPERIMENT_DIR"; sha256sum --check "$ONSET_REPAIR_DIR/patched_SHA256SUMS.txt")
