#!/usr/local_rwth/bin/zsh
# ask for number of nodes
#SBATCH --nodes=2
#SBATCH --exclusive
# ask for tasks (MPI ranks)
#SBATCH --ntasks=48
# ask for threads per task (which is 1 thread per core on one socket on CLAIX18)
#SBATCH --cpus-per-task=1
# ask for memory (per node)
#SBATCH --mem=16G   #M is the default and can therefore be omitted, but could also be K(ilo)|G(iga)|T(era)
# partition
#SBATCH --partition=c16m
# project
#SBATCH --account=jara0001
# name the job
#SBATCH --job-name=SWE_CHAMELEON
# declare the merged STDOUT/STDERR file
#SBATCH --output=batch/output/swe_charm++_output.%J.txt

### Change to the work directory
source /home/sc427635/.zshrc
cd /home/sc427635/master/swe-benchmark

### Load modules

### Compile
#CHARM_PATH=/home/sc427635/sw/charm make charm

### Set environment variables
export OMP_NUM_THREADS=1
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=auto
export I_MPI_DEBUG=5
export OMP_PROC_BIND=close

### Execute your application
echo /home/sc427635/sw/charm/bin/charmrun +p48 ./cpu_set_wrapper.sh ./build/SWE_intel_release_charm_omp_hybrid -t 1 -n 1 -x 1024 -y 1024 -o ./output/charm_batch
/home/sc427635/sw/charm/bin/charmrun +p48 ./cpu_set_wrapper.sh ./build/SWE_intel_release_charm_omp_hybrid -t 1 -n 1 -x 1024 -y 1024 -o ./output/charm_batch
