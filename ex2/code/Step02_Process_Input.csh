#!/bin/csh

#SBATCH --constraint=cpu
#SBATCH --job-name=step02         ## job_name
#SBATCH --account=m4267           ## project_name 
#SBATCH -q regular
#SBATCH --time=05:00:00           ## time_limit
#SBATCH --nodes=1                 ## number_of_nodes                                                                                              
#SBATCH --output=mat.stdout1      ## job_output_filename
#SBATCH --error=mat.stderr1       ## job_errors_filename

ulimit -s unlimited
module load matlab

matlab  -nodisplay -nosplash <Step02_Process_Input.m> Step02_Process_Input.m.log
