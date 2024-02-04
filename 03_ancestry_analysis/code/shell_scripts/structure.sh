#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition Lewis
#SBATCH --nodes=1
#SBATCH --ntasks=1  # used for MP#SBATCH -e error_%A_%a.err # Standard errorI codes, otherwise leav$
##SBATCH --cpus-per-task=12  # cores per task
#SBATCH --mem-per-cpu=16G  # memory per core (default is 1GB/core)
#SBATCH --time 2-00:00  # days-hours:minutes
#SBATCH --qos=normal
#SBATCH --array=1-260

## labels and outputs
#SBATCH --job-name=smb_fit_structure_jgunn
#
#SBATCH -o test_%A_%a.out # Standard output
#SBATCH -e error_%A_%a.err # Standard error

## notifications
#SBATCH --mail-user=jcg5g9@mail.missouri.edu  # email address for notifications
#SBATCH --mail-type=END,FAIL  # which type of notifications to send
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"


# load packages
module load rss/rss-2020
module load structure/structure-2.3.4

COMMANDA=`head -n ${SLURM_ARRAY_TASK_ID} ../batch_cmd_lists/structure_batch_cmd_list.txt | tail -n 1`
eval $COMMANDA


echo "### Ending at: $(date) ###"
