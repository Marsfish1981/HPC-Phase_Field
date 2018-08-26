#!/bin/bash
#BSUB -q hpc_linux
#BSUB -n 12
#BSUB -o output.%J
#BSUB -e error.%J

# Total number of computing process
export MPI_NUM_PROCESS=12
export I_MPI_DEBUG=5

rm -rf ./hosts
touch ./hosts

#construct the hosts file for the job

j=''
k=0
for i in `echo $LSB_HOSTS`
do
     echo $i>>./hosts
done

mpirun -machinefile ./hosts -np $MPI_NUM_PROCESS ./main.Linux.Intel.mpi.exe inputs_2d seed_def temp_def


