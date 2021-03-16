#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Tue Nov 03 2020 00:33:23 GMT+1100 (Australian Eastern Daylight Time)

# Partition for the job:
#SBATCH --partition=physical

# Multithreaded (SMP) job: must run on one node
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="N30D2L2"

# The project ID which this job should run under:
#SBATCH --account="punim0609"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


# The amount of memory in megabytes per process in the job:
#SBATCH --mem=2048

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL

# Use this email address:
#SBATCH --mail-user=qiuyi.li1993@gmail.com

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=4-0:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The job command(s):
module load foss/2019b
module load python/3.7.4
python3 submit.py -n N30D2L2