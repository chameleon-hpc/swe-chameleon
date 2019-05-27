#!/usr/local_rwth/bin/zsh
# ask for four tasks (which is 4 MPI ranks)
#SBATCH --ntasks=4
# ask for 24 threads per task (which is 1 thread per core on one socket on CLAIX18)
#SBATCH --cpus-per-task=24
#################
# ATTENTION !!! #
#################
# Divide the needed memory per task through the cpus-per-task, as slurm requests memory per cpu, not per task !
# Example:
# You need 24 GB memory per task, you have 24 cpus per task ordered
# order 24GB/24 -> 1G memory per cpu (i.e., per thread)
#SBATCH --mem-per-cpu=1G   #M is the default and can therefore be omitted, but could also be K(ilo)|G(iga)|T(era)
# name the job
#SBATCH --job-name=SWE_CHAMELEON
# declare the merged STDOUT/STDERR file
#SBATCH --output=swe_chameleon_output.%J.txt

### Change to the work directory
source /home/sc427635/.zshrc
cd /home/sc427635/master/swe-benchmark
make mpi_hybrid

### Execute your application
echo $FLAGS_MPI_BATCH
mpiexec $FLAGS_MPI_BATCH ./build/SWE_intel_release_mpi_omp_hybrid -t 1 -n 1 -x 1000 -y 1000 -o ./output/mpi_batch
