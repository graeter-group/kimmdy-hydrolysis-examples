#!/bin/env bash

ml purge
ml use /hits/fast/mbm/hartmaec/sw/easybuild/modules/all
ml GROMACS/2023.3-foss-2023a-CUDA-12.1.1-PLUMED-2.9.0

source ../../.venv-cascade/bin/activate

