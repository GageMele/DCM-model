#!/bin/bash
# Monitor every registered DCM campaign, including previously started tests.
set -euo pipefail
root="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
python3 "$(dirname -- "${BASH_SOURCE[0]}")/monitor.py" "$root" "$@"
