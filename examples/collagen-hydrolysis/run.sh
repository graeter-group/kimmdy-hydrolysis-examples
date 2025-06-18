# %%
kimmdy --input ./kimmdy_eq-0nN.yml

# %%
kimmdy --input ./kimmdy_eq-1nN.yml

# %%
cd ./run_eq_collagen_0nN/2_pull/

# %%
gmx help sasa

# %%
gmx sasa -f pull.xtc -s pull.tpr -o sasa.xvg\
  -b 0\
  -dt 10\
  -oa\
  -surface '! group "Water_and_ions"'\
  -output '! group "Water_and_ions"'

# %%
gmx make_ndx -f pull.tpr

# %%
head -n40 sasa.xvg

# %%
head -n40 ./atomarea.xvg

# %%
vmd ./assets/nvt.gro

# %%
tail -n 50 ./run_eq_collagen_0.8nN_shear/2_pull/pull.log


