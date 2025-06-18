#!/bin/env bash

if [[ "$(hostname)" == "cascade"* ]]; then
  ## cluster
  echo "cascade"
  ml purge
  ml use /hits/fast/mbm/hartmaec/sw/easybuild/modules/all
  ml GROMACS/2023.3-foss-2023a-CUDA-12.1.1-PLUMED-2.9.0
  gmx -version
elif [[ "$(hostname)" == "pop-desktop" ]]; then
  echo "local ws"
  # source /usr/local/gromacs/bin/GMXRC
  source /usr/local/gromacs-plumed/bin/GMXRC
  # gmx -version
  # local workstation
elif [[ "$(hostname)" == "pop-laptop" ]]; then
  # laptop
  echo "laptop"
  source /usr/local/gromacs/bin/GMXRC
else
  ## workstation
  echo "workstation"
  source /sw/mbm/riedmiki/plumed2/sourceme.sh
  source /sw/mbm/riedmiki/gromacs-2024.2_gpu_plumed_install/bin/GMXRC.bash
  gmx -version
fi

