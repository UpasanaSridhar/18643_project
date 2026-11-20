#!/usr/bin/env bash
set -euo pipefail
# Run from this starter directory. Requires git, make and a supported Linux toolchain.
if [[ ! -d quickLD ]]; then git clone https://github.com/pephco/quickLD.git quickLD; fi
cd quickLD
./install_blis
./config install blis
./test_run full blis
