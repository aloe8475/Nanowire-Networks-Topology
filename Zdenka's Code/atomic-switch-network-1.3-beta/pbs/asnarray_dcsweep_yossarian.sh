#/usr/bin/bash
#PBS -t 1-256
#PBS -N asn2dsweep
### specify the queue: 
#PBS -q yossarian
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -l mem=4gb
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
matlab -nodesktop -nosplash -nodisplay -singleCompThread -r  "addpath(genpath('/headnode2/paula123/Code/AtomicSwitchNetworks/atomic-switch-network/examples')); script_cluster_run_dcwait_duration_simulation('100', 'max_av', [], [])
; exit" 

##### Execute Program #####
date

##### Exit successfully 
exit 0
