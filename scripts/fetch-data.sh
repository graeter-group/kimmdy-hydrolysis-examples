
# %%
# download data from ao
rsync --exclude='*.trr' --exclude='*.xtc*' -av ao:/hits/fast/mbm/buhrjk/phd/kimmdy-examples/examples/triplehelix-hydrolysis/run_* ./examples/triplehelix-hydrolysis/


# %%
rsync -av ao:/hits/fast/mbm/buhrjk/phd/kimmdy-examples/examples/gly-hydrolysis/run ./examples/gly-hydrolysis/

# %%
rsync -av ao:/hits/fast/mbm/buhrjk/phd/kimmdy-examples/examples/triplehelix-hydrolysis/run_eq-triple_1nN/2_pull ./examples/triplehelix-hydrolysis/run_eq-triple_1nN/

# %%
rsync -av ao:/hits/fast/mbm/buhrjk/phd/kimmdy-examples/results/* ./results/
