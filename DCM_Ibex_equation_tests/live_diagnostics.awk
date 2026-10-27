# Stream only sampler summaries; incomplete parameter rows are skipped.
BEGIN { phase = "warmup"; chain = "?" }
/^#[[:space:]]*id[[:space:]]*=/ {
  value = $0; sub(/.*=[[:space:]]*/, "", value); chain = value + 0
}
/^#[[:space:]]*max_depth[[:space:]]*=/ {
  value = $0; sub(/.*=[[:space:]]*/, "", value); cap = value + 0
}
/^# Adaptation terminated/ { phase = "sampling" }
/^#/ { next }
$1 == "lp__" {
  columns = NF
  for (i = 1; i <= NF; i++) {
    if ($i == "treedepth__") td = i
    if ($i == "divergent__") dv = i
    if ($i == "n_leapfrog__") lf = i
  }
  next
}
td && dv && lf && cap && NF == columns {
  n[phase]++
  depth = $(td) + 0
  hit = (depth >= cap)
  hits[phase] += hit
  divergences[phase] += $(dv)
  leaps[phase] += $(lf)
  if (depth > maximum[phase]) maximum[phase] = depth
  slot = (n[phase] - 1) % 50
  recent[phase] -= ring[phase, slot]
  ring[phase, slot] = hit
  recent[phase] += hit
}
END {
  for (j = 1; j <= 2; j++) {
    p = (j == 1 ? "warmup" : "sampling")
    if (!n[p]) continue
    window = (n[p] < 50 ? n[p] : 50)
    printf "%4d %5s %-8s %6d %3d %3d %6d %7.1f %8.1f %5d %10.1f\n",
      case_id, chain, p, n[p], cap, maximum[p], hits[p],
      100 * hits[p] / n[p], 100 * recent[p] / window,
      divergences[p], leaps[p] / n[p]
  }
  if (!n["warmup"] && !n["sampling"])
    printf "Case %s chain %s: no complete saved draws visible yet.\n", case_id, chain
}
