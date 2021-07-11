#!/bin/bash -e
#SBATCH --job-name=closedLoop_SpO2L60            # job name (shows up in the queue)
#SBATCH --account=uoa00596     # Project Account
#SBATCH --time=168:00:00         # Walltime (HH:MM:SS)
#SBATCH --mem=8000      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --cpus-per-task=8       # number of cores per task (e.g. OpenMP)
#SBATCH --profile=all		#Profiles job.
#SBATCH --array=1-16                # Array jobs

export model_name="closedLoop_SpO2L"
export para_name="parL2"
export root_path="/nesi/nobackup/uoa00596/lungPacemaker/" # Everything happens in here. Having this a variable makes it easy to move where you run your jobs. 
export script_path="${root_path}model/"  #This is where all the scripts are.
export data_path="${root_path}data/$SLURM_JOB_NAME"  #This is where all the scripts are.
mkdir -vp ${data_path}  #This will create the directory for the data
export working_path="${root_path}Working/$SLURM_JOB_NAME/run_${SLURM_ARRAY_TASK_ID}/" #This the jobs put their individual files that need to be kept seperate.
mkdir -vp ${working_path}  #This will create the directory and move you into it.
cp ${script_path}${model_name}.mdl  ${working_path}${model_name}${SLURM_ARRAY_TASK_ID}.mdl  		#This will take a copy of the simulink model into the running directory.
cp ${script_path}${para_name}.mat   ${working_path}${para_name}.mat  		#This will take a copy of the parameter into the running directory.
export TMPDIR=${working_path} #Stop matlab temp files clashing.
module load MATLAB/2019b
matlab -nojvm -nodisplay -r "run('${script_path}runClosedLoopSo2_60.m');exit;"
rm -r ${working_path}

