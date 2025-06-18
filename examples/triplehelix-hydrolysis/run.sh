# %% 
#
#  1  Bond             2  Morse            3  Angle            4  Proper-Dih.   
#  5  Per.-Imp.-Dih.   6  LJ-14            7  Coulomb-14       8  LJ-(SR)       
#  9  Disper.-corr.   10  Coulomb-(SR)    11  Coul.-recip.    12  Potential     
# 13  Kinetic-En.     14  Total-Energy    15  Conserved-En.   16  Temperature   
# 17  Pres.-DC        18  Pressure        19  dVremain/dl     20  Box-X         
# 21  Box-Y           22  Box-Z           23  Volume          24  Density       
# 25  pV              26  Enthalpy        27  Vir-XX          28  Vir-XY        
# 29  Vir-XZ          30  Vir-YX          31  Vir-YY          32  Vir-YZ        
# 33  Vir-ZX          34  Vir-ZY          35  Vir-ZZ          36  Pres-XX       
# 37  Pres-XY         38  Pres-XZ         39  Pres-YX         40  Pres-YY       
# 41  Pres-YZ         42  Pres-ZX         43  Pres-ZY         44  Pres-ZZ       
# 45  #Surf*SurfTen   46  Box-Vel-XX      47  Box-Vel-YY      48  Box-Vel-ZZ    
# 49  T-Protein                           50  T-non-Protein                     
# 51  Lamb-Protein                        52  Lamb-non-Protein                  

# %%
# with pairs, PR current state

# %%
kimmdy --input ./kimmdy-full-no-plumed.yml

# %%
kimmdy --input ./kimmdy-relax-test.yml

# %%
kimmdy --input kimmdy-just-reactions-single_1nN_7.4pH.yml

# %%
kimmdy --input ./kimmdy-just-reactions-single_0nN_7.4pH.yml

# %%
kimmdy-analysis trjcat ./run_single_0nN_7.4pH_001 -o -l -g0

# %%
vmd run_single_0nN_7.4pH_001/analysis/concat.gro run_single_0nN_7.4pH_001/analysis/concat.xtc

# %%
ls ./run_single_0nN_7.4pH_001

# %% 
kimmdy

# %% 
i="_001"

# %%
kimmdy-analysis energy run_prod$i -o --terms 'coulomb-(sr)' potential

# %% 
kimmdy-analysis energy run_prod$i -o --terms potential lj-14 'LJ-(SR)' coulomb-14 'coulomb-(sr)'

# %%
kimmdy-analysis trjcat run_prod$i -o -g 0 -l

# %% 
n=5
vmd ./run_00$n/4_relax/relax.gro ./run_00$n/4_relax/relax.xtc

# %% 
vmd ./run_00$1/0_setup/eq.gro ./run_00$n/4_relax/relax.xtc

# send the next to VMD terminal

# %%
mol color Name;
mol modselect 0 0 "not water and not ion or same residue as serial 74481";
mol modmaterial 0 0 AOChalky;
mol modstyle 0 0 CPK;
display resetview;
display update;
# %%
animate forward;
# %%
animate pause;

# %%
exit

# %%

# %%

# %%

# %

#
# %% 
cp ./run_001/3_apply_recipe/topol_after.top manual/
cp ./run_001/3_apply_recipe/topol_relax.top manual/
cp ./run_001/3_apply_recipe/topol_before_with_solvent_bonds.top manual/
cp ./run_001/0_setup/eq.gro manual/eq.gro
ln -sr ./relax.mdp manual/relax.mdp
ln -sr ./index.ndx manual/index.ndx


# %%
rm -r run_*

kimmdy

# %%
ls ./run_eq_single_*nN/2_pull/.kimmdy.bondstats

# %%
ls ./run_react_*
