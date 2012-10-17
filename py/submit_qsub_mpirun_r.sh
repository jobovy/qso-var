#!/bin/sh
#Run qsub -l h_rt=1000:00:00
#$ -pe orte 9
#$ -cwd
#$ -V
#$ -R y
#$ -r y
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
export MPI=/usr/local/openmpi/gcc/x86_64
export PATH=${MPI}/bin:${PATH}
export LD_LIBRARY_PATH=${MPI}/lib
mpirun -x PYTHONPATH /home/bovy/local/bin/python skewQSOCluster.py ../skew/powerlawSF_constmean_r_skew_1000.sav -b r -t powerlawSF --mean=const -f ../fits/powerlawSF_constmean_r.sav -n 1000 --saveevery=25
