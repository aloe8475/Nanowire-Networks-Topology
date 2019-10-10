#/usr/bin/bash 
## Specify your shell
## Required PBS Directives ----------------------------------------------------
#PBS -t 68-92
#PBS -N asnseed
# use specific nodes
#PBS -q yossarian
#PBS -j oe
#PBS -l nodes=1:ppn=4
#PBS -l walltime=4:00:00
#PBS -l mem=16gb
#PBS -m abe
#PBS -M paula.sanz-leon@sydney.edu.au
# pass PYTHONPATH environment variable to the batch process
#PBS -v PATH
# use the submission environment
#PBS -V
#PBS -v PYTHONPATH

## Execution Block ----------------------------------------------------
#  Environment Setup 
cd $PBS_O_WORKDIR 

# create a job-specific subdirectory based on JOBID and cd to it
# JOBID=`echo ${PBS_JOBID} | cut -d '.' -f 1`
# mkdir -p ${JOBID}
# cd $JOBID

#  Launching  ----------------------------------------------------
LAUNCH_REPO=~/Code/AtomicSwitchNetworks/atomic-switch-network
LAUNCH_SUBDIR=asn/connectivity
THIS_SCRIPT=generate_nanowires_network.py


echo "####################################################"
echo "User: $PBS_O_LOGNAME"
echo "Batch job started on $PBS_O_HOST"
echo "PBS job id: $PBS_JOBID"
echo "PBS array id: $PBS_ARRAYID"
echo "PBS job name: $PBS_JOBNAME"
echo "PBS working directory: $PBS_O_WORKDIR"
echo "Job started on" `hostname` `date`
echo "Current directory:" `pwd`
echo "PBS environment: $PBS_ENVIRONMENT"
echo "####################################################"

echo "####################################################"
echo "Full Environment:"
printenv
echo "####################################################"

echo "####################################################"
echo "The following script is being executed: $THIS_SCRIPT"
echo "####################################################"


echo "Git ref is:"
cat ${LAUNCH_REPO}/.git/refs/heads/master
echo "##########################################################"

echo "The Job is using the following python:"
which python 
echo "##########################################################"

echo "The Job is being executed on the following node:"
cat ${PBS_NODEFILE}
echo "##########################################################"

# Execute script
python $LAUNCH_REPO/$LAUNCH_SUBDIR/$THIS_SCRIPT --nwires 16384 --seed $PBS_ARRAYID

echo "Job Finished: " `date`
echo "##########################################################"
exit 0
