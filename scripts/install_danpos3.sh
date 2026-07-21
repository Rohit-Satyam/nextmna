#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="${ROOT}/bin/DANPOS3"
if [[ -f "${TARGET}/danpos.py" ]]; then
  echo "DANPOS3 already present: ${TARGET}/danpos.py"
  exit 0
fi
rm -rf "${TARGET}"
git clone https://github.com/sklasfeld/DANPOS3.git "${TARGET}"
echo "Installed DANPOS3 at ${TARGET}"
