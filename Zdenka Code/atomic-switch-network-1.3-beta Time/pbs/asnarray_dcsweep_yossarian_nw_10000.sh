#/usr/bin/bash
#PBS -t 1-1024
#PBS -N asndc
### specify the queue: 
#PBS -q yossarian
#PBS -j oe
#PBS -l nodes=1:ppn=2
#PBS -l walltime=72:00:00
#PBS -l mem=32gb
#PBS -m abe
#PBS -M p.sanz-leon@sydney.edu.au
#PBS -V
### cd to directory where the job was submitted:
cd $PBS_O_WORKDIR
date
echo "----------------"
echo "PBS job running on: `hostname`"
echo "in directory:       `pwd`"
echo "----------------"
echo "Full Environment:"
printenv
echo "####################################################"

### Pass it to matlab 2016
matlab -nodesktop -nosplash -nodisplay -singleCompThread -r  "addpath(genpath('/headnode1/paula123/Code/AtomicSwitchNetworks/atomic-switch-network/examples')); script_cluster_run_dcwait_simulation('10000'); exit" 

##### Execute Program #####
date

##### Exit successfully 
exit 0