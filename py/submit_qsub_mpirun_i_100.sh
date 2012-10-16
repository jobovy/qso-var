#!/bin/sh
#Run qsub -l h_rt=10:00:00
#$ -pe orte 250
#$ -cwd
#$ -V
#$ -R y
#$ -r y
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
export MPI=/usr/local/openmpi/gcc/x86_64
export PATH=${MPI}/bin:${PATH}
export LD_LIBRARY_PATH=${MPI}/lib
mpirun -x PYTHONPATH /home/bovy/local/bin/python skewQSOCluster.py ../skew/powerlawSF_constmean_i_skew.sav -b i -t powerlawSF --mean=const -f ../fits/powerlawSF_constmean_i.sav -n 100 --saveevery=5 --split=250
