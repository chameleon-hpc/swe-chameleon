#!/usr/local_rwth/bin/zsh
#SBATCH --exclusive
# ask for memory (per node)
#SBATCH --mem=16G   #M is the default and can therefore be omitted, but could also be K(ilo)|G(iga)|T(era)

### Change to the work directory
source /home/sc427635/.zshrc
cd /home/sc427635/master/swe-benchmark

### Load modules

### Compile
#CHARM_PATH=/home/sc427635/sw/charm make charm

### Set environment variables
#export OMP_NUM_THREADS=1
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=auto
export I_MPI_DEBUG=5
export OMP_PROC_BIND=close

### Execute your application
echo /home/sc427635/sw/charm/bin/charmrun +p$NUM_TASKS_PER_NODE ++mpiexec ./batch/cpu_set_wrapper.sh ./build/SWE_intel_release_charm_omp_hybrid -t 1 -n 1 -x $SIZE -y $SIZE -o ./output/charm_batch
/home/sc427635/sw/charm/bin/charmrun +p$NUM_TASKS_PER_NODE ++mpiexec ./batch/cpu_set_wrapper.sh ./build/SWE_intel_release_charm_omp_hybrid -t 1 -n 1 -x $SIZE -y $SIZE -o ./output/charm_batch
