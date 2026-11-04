#!/usr/bin/env bash
# Isolated user-local Mesa runtime; no permanent system or driver modifications.
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
mesa="$HOME/.local/share/ksi-mesa-opencl"
python="$HOME/.local/share/enz-kinetics-mm/bin/python"
for item in "$mesa/usr/lib/clc" "$mesa/usr/lib/llvm-20" "$mesa/etc/OpenCL/vendors" "$python"; do
  if [[ ! -e "$item" ]]; then echo "Missing dependency: $item" >&2; exit 1; fi
done
exec bwrap --ro-bind / / --dev-bind /dev /dev --tmpfs /usr/lib \
  --ro-bind /usr/lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu \
  --ro-bind "$mesa/usr/lib/clc" /usr/lib/clc \
  --ro-bind "$mesa/usr/lib/llvm-20" /usr/lib/llvm-20 \
  --bind /tmp /tmp --bind "$root/enzyme-kinetics/outputs" "$root/enzyme-kinetics/outputs" \
  --setenv OCL_ICD_VENDORS "$mesa/etc/OpenCL/vendors" \
  --setenv LD_LIBRARY_PATH "$mesa/usr/lib/x86_64-linux-gnu" \
  --setenv RUSTICL_ENABLE radeonsi --setenv MESA_SHADER_CACHE_DISABLE true \
  "$python" "$@"
