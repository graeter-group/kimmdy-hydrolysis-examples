#!/bin/bash -l
#SBATCH -o collagen_%j.out
#SBATCH -e collagen_%j.err
#SBATCH -D ./
#SBATCH -J collagen
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=none
#SBATCH --time=24:00:00

module purge
module load gromacs/2024.3

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores

job=$1

MDRUN='gmx_mpi mdrun'

export CYCLE=23
START=$(date +"%s")

if [ -e ${job}.tpr ]; then
  srun $MDRUN -deffnm ${job} -px ${job}_pullx -pf ${job}_pullf -cpi ${job}.cpt -ntomp $OMP_NUM_THREADS -dlb auto -maxh 23.9
else
  srun --ntasks=1 gmx_mpi grompp -f ${job}.mdp -c npt.gro -t ${job}.trr -n collagen.ndx -p collagen.top -o ${job}.tpr -maxwarn 2
  srun $MDRUN -deffnm ${job} -v -ntomp $OMP_NUM_THREADS -dlb auto -maxh 23
fi

END=$(date +"%s")
echo "$(((END-START))) seconds ran"
echo "$(((END-START)/3600)) full hours ran"
let "CYCLE--"

if [ $(((END-START)/3600)) -lt $CYCLE ]; then
  echo "last cycle was just $(((END-START)/3600))h long and therefore finito"
else
  echo "script: resubmitting"
  sbatch --job-name=$SLURM_JOB_NAME viper-jobscript.sh $job
fi
