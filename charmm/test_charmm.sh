#!/usr/bin/env bash

# TODO: want to test with --openmm-ffxml once OpenMM is fixed
python test_charmm.py tests/*/*.yaml --charmm --openmm-charmm --openmm-ffxml-fix-impropers --perturb-replicates 3 --perturb-seed 291700478 "$@"
