#!/bin/bash
#PBS -l select=1:ncpus=32
#PBS -N 3d_npt_394
cd $PBS_O_WORKDIR
#PBS -j oe
#PBS -k oed
#PBS -q workq

#module purge
#module load mpi/openmpi-4.1.4



export OMP_NUM_THREADS=1

for i in $(seq 125600 1 125799)
do
START=$i;
STOP=$(( ${START} + 1 ))

for ((  iter=${START}; iter < ${STOP}; iter++ ))
#for i in $(seq 1 1 20)
do
		
    	INIT=$( printf '%06d' ${iter} )

         ~/apps/bin/mpirun -np 32 ~/lammps-23Jun2022/src/lmp_mpi  \
        -in in.3d_KABLJM_PBC_zero_pressure \
        -var SEED `bash -c 'echo $(( $(od -An -tu4 -N4 /dev/urandom) >> 3 ))'` \
        -var INIT ${INIT} &
done
wait

done
