#!/bin/sh
#Run qsub -t 1-251 -l h_rt=10:00:00
#$ -cwd
#$ -V
#$ -R y
#$ -r y
band=u
nsamples=100
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
export MPI=/usr/local/openmpi/gcc/x86_64
export PATH=${MPI}/bin:${PATH}
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MPI}/lib
/home/bovy/local/bin/python skewQSO.py ../skew-wpv/powerlawSF_constmean_${band}_skew_${SGE_TASK_ID}.sav -b ${band} -t powerlawSF --mean=const -f ../fits/powerlawSF_constmean_${band}.sav -n ${nsamples} --saveevery=5 --split=250 --ra=${SGE_TASK_ID}
