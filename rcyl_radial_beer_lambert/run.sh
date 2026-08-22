#!/usr/bin/env bash
set -euo pipefail
test_root=$(cd "$(dirname "$0")/.." && pwd)
exec "${PYTHON:-/home/tomzhu0225/venvs/warpx-viz/bin/python}" "$test_root/run_all.py" --case rcyl_radial_beer_lambert
