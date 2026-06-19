#!/usr/bin/env bash
# Thin CLI wrapper for the Julia Baysor install (no PackageCompiler binary needed).
# Invokes Baysor.run_cli(ARGS) in the isolated project/depot built for headroom-linux.
set -euo pipefail
export JULIA_DEPOT_PATH="$HOME/baysor_depot"
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-4}"
exec "$HOME/julia_install/julia-1.10.11/bin/julia" --project="$HOME/baysor_proj" \
  -e 'using Baysor; Baysor.command_main(ARGS)' "$@"
