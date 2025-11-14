#!/bin/env bash

#%%

#%%
gro="./examples/triplehelix-hydrolysis/run_eq-triple_0nN/2_pull/pull.gro"
vmd $gro

#%%
vmd -e ./scripts/triplehelix-gly-pro.vmd

#%%
vmd -e ./scripts/triplehelix-gly-pro-pull.vmd


#%%
#%%
#%%
#%%
#%%
