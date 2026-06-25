#!/usr/bin/env bash
# run_all_sims.sh — Run all simulation executables sequentially and save outputs
set -euo pipefail

# Arguments to pass to every simulation
args="0.0 0.4 30 100000"

# Create output folder (optional, keeps things tidy)
mkdir -p logs

# Loop through executables that match your naming pattern
for exe in simulations_*; do
    # Skip if not executable
    [[ -x "$exe" ]] || continue

    # Strip leading './' if present
    base="${exe#./}"

    # Construct output filename (replace 'simulations_' with '')
    name="${base#simulations_}"
    out="${name}_sims.txt"

    echo ">>> Running: $exe $args > $out"
    "./$exe" $args > "logs/$out"
done

echo "✅ All simulations complete. Results saved in ./logs/"

