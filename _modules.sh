#!/bin/env bash

if [[ "$(hostname)" == "cascade"* ]]; then
  ## cluster
  ml purge
  ml use /hits/fast/mbm/hartmaec/sw/easybuild/modules/all
  ml GROMACS/2023.3-foss-2023a-CUDA-12.1.1-PLUMED-2.9.0
  gmx -version
elif [[ "$(hostname)" == "pop-desktop" ]]; then
  # source /usr/local/gromacs/bin/GMXRC
  source /usr/local/gromacs-plumed/bin/GMXRC
  # gmx -version
  # local workstation
elif [[ "$(hostname)" == "pop-laptop" ]]; then
  # laptop
  source /usr/local/gromacs/bin/GMXRC
elif [[ "$(hostname)" == "pop-os" ]]; then
  # laptop
  :
else
  ## workstation
  source /sw/mbm/riedmiki/plumed2/sourceme.sh
  source /sw/mbm/riedmiki/gromacs-2024.2_gpu_plumed_install/bin/GMXRC.bash
  gmx -version
fi

