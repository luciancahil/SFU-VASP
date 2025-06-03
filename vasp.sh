#!/bin/bash
# The session has one node with 24 cores and 72 GB memory, and a 4 hour time limit.
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=64GB
#SBATCH --account=def-samiras-ab
#SBATCH --mail-user=royhe@student.ubc.ca
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT,TIME_LIMIT_90,FAIL
#SBATCH --output=DFT.txt
#SBATCH --error=DFT_Error.txt



CRYSTAL=$1
SETTINGS=$2

# Start timing the entire job
start_time=$(date +%s)

# This loads the preliminary modules necessary for Vasp, then modules necessary for ASE environment, and finally loads Vasp
module purge
module load StdEnv/2023 intel/2023 intelmpi/2021.9.0
module load python/3.11.5 scipy-stack
module load vasp/6.4.2

# This checks if the modules are loaded correctly. If not, the script is prevented from further running.
if [[ $(module list | grep 'intel/2023') == ""  || $(module list | grep 'python/3.11.5') == "" || $(module list | grep 'vasp/6.4.2') == "" ]]; then
	module list
	echo "Your modules are not loaded correctly for VASP. Cancelling job... "

else
	module list
	echo "Your modules are loaded correctly for VASP. Proceeding to activate ASE..."
fi

# This loads the python virtualenv using ASE
if test -e "$HOME/.bashrc"; then
	source "$HOME/.bashrc"
fi

load_ase() {
	source $ASE_ENV
}
load_ase

echo "####################################
# We are now Running job ${SLURM_JOB_ID}! #
####################################"

# Get the job's cwd
cwd=$(pwd -P)


# Start the job
time python3 ./../../flow.py --crystal $CRYSTAL

# End timing the entire job
end_time=$(date +%s)

# Printing the task complete message allows the job restarter to find jobs which never completed correctly.
elapsed_time=$((end_time - start_time))
elapsed_time_hms=$(printf "%02d:%02d:%02d" $((elapsed_time / 3600)) $(( (elapsed_time % 3600) / 60 )) $((elapsed_time % 60)))
echo "Task ${SLURM_JOB_ID} completed in $elapsed_time_hms (HH:MM:SS)." 
echo $(pwd)