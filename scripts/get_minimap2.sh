#!/usr/bin/env bash
# Builds minimap2 from source into third_party/minimap2 (for the accuracy comparison).
# Skip this if minimap2 is already on your PATH (e.g. `brew install minimap2`, conda).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p third_party
if [[ ! -d third_party/minimap2 ]]; then
  git clone --depth 1 https://github.com/lh3/minimap2.git third_party/minimap2
fi
cd third_party/minimap2
case "$(uname -m)" in
  arm64|aarch64) make -j4 arm_neon=1 aarch64=1 ;;
  *) make -j4 ;;
esac
echo "Built: third_party/minimap2/minimap2 ($(./minimap2 --version))"
