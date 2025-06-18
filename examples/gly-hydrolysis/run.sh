# %%
kimmdy --input ./kimmdy-eq.yml

# %%
kimmdy --input ./kimmdy.yml

# %%
kimmdy-analysis trjcat ./run_prod -f xtc -g0 -l -o

# %%
vmd -e vis-prod.vmd


# %%
cd render
ffmpeg -i untitled.%5d.ppm -c:v libx264 -pix_fmt yuv420p output.mp4
mv output.mp4 ..
cd ..


# groups for the chosen reaction to follow their distances
# (1-based top IDs) 12⚡14 f 71637➡12 71638➡14 71639➡14 

# %%
gmx distance -f ./run_prod/analysis/concat.xtc -s ./run_prod/analysis/concat.gro \
  -oall\
  -select -sf ./select.sel


