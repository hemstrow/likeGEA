#!/bin/bash -l

#SBATCH -t 2-24:00:00
#SBATCH -J likeGEA
#SBATCH --array=1-ARRAY_REPLACE
#SBATCH --mem=100G
#SBATCH -n 30
#SBATCH -N 1
#SBATCH --output=./slurm_%x_%A_%a.out
#SBATCH --error=./slurm_%x_%A_%a.err

##############################################################
# IMPORTANT NOTES:
# This script is an EXAMPLE and will need to be edited for use on any actual
# HPC, slurm-enabled clusters. At minimum, the following should be changed:
#   * Partitions, quality-of-service ("qos"), and account #SBATCH tags.
#     Depending on your HPC setup, these will probably need to set to enable the
#     job to actually run.
#   * #SBATCH --mem, -n, and -t tags. These are set according to your project's
#     memory use and parallel ("par" argument) needs.
#   * Module activation. Depending on your specific HPC set up, you may need to
#     activate some modules (like conda/mamba).
#
# As noted in the example_slurm_job() function documentation, this all assumes
# that a correct likeGEA conda environment is set-up and ready-to-use. The
# package's "setup_likeGEA_env" function is a good way to acheive this without
# versioning/compilation errors.
#
# Remember that likeGEA can automatically partition/split beagle files for you
# with the "get_chunks()" function. You can make a list of these files with
# "ls FILES*.beagle.gz > beagle_list.txt" where FILES is your file prefix.
#
# Please also confirm that all file paths are set appropriately below!
###############################################################

# Input files
beagle_list=BEAG_REPLACE
metadata=META_REPLACE
parameter_file=PARAM_REPLACE
dir_prefix=DIR_REPLACE
env=ENV_REPLACE
execution_Rscript=RSCRIPT_REPLACE

conda activate $env

string="sed -n ${SLURM_ARRAY_TASK_ID}p ${beagle_list}"
beagle=$($string)

# extract the beagle basename
bn=$(basename $beagle .gz)
bn=$(basename $bn .beagle)

# adjust file paths to absolute since we are moving dirs
beagle=$(realpath $beagle)
metadata=$(realpath $metadata)
parameter_file=$(realpath $parmeter_file)
execution_Rscript=$(realpath $execution_Rscript)

# make a working directory and move to it
working_dir="${dir_prefix}/partition_${SLURM_ARRAY_TASK_ID}_${bn}"

mkdir $working_dir

cd $working_dir

# Run the actual R script
Rscript $execution_Rscript $beagle $metadata $parameter_file $SLURM_ARRAY_TASK_ID

