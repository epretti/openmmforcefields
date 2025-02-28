#!/usr/bin/env bash

python test_charmm.py tests/16*/*.yaml --charmm --openmm-charmm --openmm-ffxml-fix-impropers --perturb-replicates 5 --perturb-seed 291700478 "$@"
