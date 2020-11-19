**to run jupyter notebook:**

Load Nomachine connection to Gateway

open 2 terminals

terminal 1 sshHPC

qjupyter

terminal 2 sshL

enter password

go to localhost:9345 on browser

load + run dask_control.ipynb

port for dask dashboard: 18778



**Other Commands:**

qclear - deletes jupyter session

qls - shows Artemis queue

/project/NASN/ - path to project storage

**move file from silo to Artemis**:

- Storage: rsync -avz --progress /import/silo2/aloe8475/CODE/Data/Functional\ Connectivity/<file name>  aloe8475@hpc.sydney.edu.au:/project/NASN/Alon/<destination folder>
- Personal: rsync -avz --progress /import/silo2/aloe8475/Documents/CODE/Analysis/Functional\ Connectivity/Functional\ Tasks/"NonlinearTransformationAnalysis.ipynb" aloe8475@hpc.sydney.edu.au:/home/aloe8475/Documents/FunctionalConnectivity/



**to kill process:**

type: ps aux| grep aloe8475

type: kill <processID>