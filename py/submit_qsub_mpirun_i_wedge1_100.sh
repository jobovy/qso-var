#!/bin/sh
#Run qsub -l h_rt=4:00:00
#$ -pe orte 250
#$ -cwd
#$ -V
#$ -R y
#$ -r y
band=i
nsamples=100
wedgerate=1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
export MPI=/usr/local/openmpi/gcc/x86_64
export PATH=${MPI}/bin:${PATH}
export LD_LIBRARY_PATH=${MPI}/lib
mpirun -x PYTHONPATH /home/bovy/local/bin/python skewQSOCluster.py ../skew/powerlawSF_constmean_${band}_skew_wedge${wedgerate}.sav -b ${band} -t powerlawSF --mean=const -f ../fits/powerlawSF_constmean_${band}.sav -n ${nsamples} --saveevery=5 --split=250 --wedge --wedgerate=${wedgerate}
