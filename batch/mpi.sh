#!/usr/local_rwth/bin/zsh
#SBATCH --exclusive
# ask for memory (per node)
#SBATCH --mem=16G   #M is the default and can therefore be omitted, but could also be K(ilo)|G(iga)|T(era)

### Change to the work directory
source /home/sc427635/.zshrc
cd /home/sc427635/master/swe-benchmark

### Load modules
module load chameleon

### Compile
#make chameleon

### Set environment variables
#export OMP_NUM_THREADS=12
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=auto
export I_MPI_DEBUG=5
export OMP_PROC_BIND=close

### Print execution statement
echo $MPIEXEC $FLAGS_MPI_BATCH ./batch/cpu_set_wrapper.sh ./build/SWE_intel_release_mpi_omp_hybrid -t 1 -n 1 -x $SIZE -y $SIZE -o ./output/mpi_batch

### Execute your application
$MPIEXEC $FLAGS_MPI_BATCH ./batch/cpu_set_wrapper.sh ./build/SWE_intel_release_mpi_omp_hybrid -t 1 -n 1 -x $SIZE -y $SIZE -o ./output/mpi_batch
