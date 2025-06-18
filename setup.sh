#!/bin/env bash

# set cwd to wherever this script is
current_dir=$(pwd)
dir=$(dirname "${BASH_SOURCE[0]}")
cd $dir

if [ ! -f examples/collagen-hydrolysis/assets/collagen.ndx ]; then 
  ./scripts/extract-assets.sh
fi

# load modules (gromacs, plumed, etc.)
if [ -f ./_modules.sh ]; then
  source ./_modules.sh
fi

# download dependencies
# choose the ones you need to work
# on and install the others from pypi
# see requirements.txt
if [ ! -d deps ]; then
  mkdir deps
  pushd deps
  git clone git@github.com:graeter-group/kimmdy.git
  git clone git@github.com:graeter-group/kimmdy-grappa.git
  git clone git@github.com:graeter-group/kimmdy-hydrolysis
  git clone git@github.com:graeter-group/kimmdy_paper_theme.git
  git clone git@github.com:caapontes/kimmdy-reactions-spec-binding.git kimmdy-binding
  popd
fi

# create virtual environment
if [[ "$(hostname)" == "cascade"* ]]; then
  if [ ! -d .venv-cascade ]; then
    python -m venv .venv-cascade
    source .venv-cascade/bin/activate
    python -m pip install --upgrade pip
    python -m pip install -r requirements.txt
  fi
else
  if [ ! -d .venv ]; then
    python -m venv .venv
    source .venv/bin/activate
    python -m pip install --upgrade pip
    python -m pip install -r requirements.txt
  fi
fi

export PS1='$ '

# activate virtual environment
if [[ "$(hostname)" == "cascade"* ]]; then
  source .venv-cascade/bin/activate
else
  source .venv/bin/activate
fi

cd $current_dir


