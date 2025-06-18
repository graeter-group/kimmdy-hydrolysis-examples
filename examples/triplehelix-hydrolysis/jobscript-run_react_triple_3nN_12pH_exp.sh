#!/bin/env bash
#SBATCH --job-name=run_react_triple_3nN_12pH_exp
#SBATCH --output=run_react_triple_3nN_12pH_exp.slurm.out
#SBATCH --error=run_react_triple_3nN_12pH_exp.slurm.err
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mincpus=40
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --gpus=2
#SBATCH --mail-type=ALL
#SBATCH --partition cascade.p
# # uncomment these to use:
# #SBATCH --mail-user=<your-email>


# Setup up your environment here
# modules.sh might load lmod modules, set environment variables, etc.
if [ -f ./_modules.sh ]; then
    source ./_modules.sh
fi

CYCLE=23
CYCLE_buffered=$(echo "scale=2; $CYCLE - 0.08" | bc)


START=$(date +"%s")

timeout ${CYCLE_buffered}h kimmdy -i kimmdy_just_reactions_triple_3nN_12pH_exp.yml --restart

END=$(date +"%s")

LEN=$((END-START))
HOURS=$((LEN/3600))

echo "$LEN seconds ran"
echo "$HOURS full hours ran"

let "CYCLE--"
if [ $HOURS -lt $CYCLE ]; then
  echo "last cycle was just $HOURS h long, KIMMDY is done."
  exit 3
else
  echo "jobscript resubmitting"
  sbatch ./jobscript-run_react_triple_3nN_12pH_exp.sh
  exit 2
fi