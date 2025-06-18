#!/bin/env bash

# %%
gmx pdb2gmx -f ./gly.pdb -o ./gly.gro -water tip3p -ignh -ff amber99sb-star-ildnp
gmx editconf -f ./gly.gro -o ./gly_box.gro -c -d 5.0 -bt triclinic
gmx solvate -cp ./gly_box.gro  -o ./gly_solv.gro -p ./topol.top
gmx grompp -f ./ions.mdp -c ./gly_solv.gro -p ./topol.top -o ./ions.tpr

echo '13' | gmx genion -s ./ions.tpr -o ./gly_solv_ions.gro -p ./topol.top -pname NA -nname CL -neutral -conc 0.15

gmx grompp -f ./minim.mdp -c ./gly_solv_ions.gro -p ./topol.top -o ./em.tpr
gmx mdrun -v -deffnm em

gmx grompp -f ./nvt.mdp -c ./em.gro -p ./topol.top -o ./nvt.tpr
gmx mdrun -v -deffnm nvt

gmx grompp -f ./npt.mdp -c ./nvt.gro -p ./topol.top -o ./npt.tpr
gmx mdrun -v -deffnm npt

echo -e 'q\n' | gmx make_ndx -f ./npt.gro
rm \#*\#

# %%
gmx pdb2gmx -f ./two-gly.pdb -o ./two-gly.gro -water tip3p -ignh -ff amber99sb-star-ildnp -p two-topol.top
gmx editconf -f ./two-gly.gro -o ./two-gly_box.gro -c -d 5.0 -bt triclinic
gmx solvate -cp ./two-gly_box.gro  -o ./two-gly_solv.gro -p ./two-topol.top
gmx grompp -f ./ions.mdp -c ./two-gly_solv.gro -p ./two-topol.top -o ./two-ions.tpr
echo '13' | gmx genion -s ./two-ions.tpr -o ./two-gly_solv_ions.gro -p ./two-topol.top -pname NA -nname CL -neutral -conc 0.15

gmx grompp -f ./minim.mdp -c ./two-gly_solv_ions.gro -p ./two-topol.top -o ./two-em.tpr
gmx mdrun -v -deffnm two-em

gmx grompp -f ./nvt.mdp -c ./two-em.gro -p ./two-topol.top -o ./two-nvt.tpr
gmx mdrun -v -deffnm two-nvt

gmx grompp -f ./npt.mdp -c ./two-nvt.gro -p ./two-topol.top -o ./two-npt.tpr
gmx mdrun -v -deffnm two-npt

echo -e 'q\n' | gmx make_ndx -f ./two-npt.gro -n ./two-index.ndx
rm \#*\#

