#!/usr/bin/env bash

set -evx

python test_charmm.py tests/01*/*.yaml --charmm --openmm-charmm --openmm-ffxml-fix-impropers --perturb-replicates 5 --perturb-seed 291700478 "$@"
