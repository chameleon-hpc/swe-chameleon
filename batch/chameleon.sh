#!/usr/local_rwth/bin/zsh
#SBATCH --exclusive

### Change to the work directory
source /home/sc427635/.zshrc
cd /home/sc427635/master/swe-benchmark

### Load modules
module load chameleon

### Compile
#make chameleon

### Set environment variables
#OMP_NUM_THREADS=$(($OMP_NUM_THREADS - 2))
OMP_NUM_THREADS=$(($OMP_NUM_THREADS - 1))
export OMP_NUM_THREADS=$OMP_NUM_THREADS
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=auto
export I_MPI_DEBUG=5
export OMP_PROC_BIND=close
export DRY_FRACTION=$DRY_FRACTION

### Print execution statement
COMMAND="mpiexec $FLAGS_MPI_BATCH ./batch/cpu_set_wrapper.sh ./build/SWE_intel_release_chameleon_omp_augrie "
COMMAND+="-t 1 -n 1 "
COMMAND+="-x $GRID_DIMENSION -y $GRID_DIMENSION "
COMMAND+="-i $BLOCK_COUNT -j $BLOCK_COUNT "
COMMAND+="-o ./output/chameleon_batch -i 200"

if [ "$EXTRA" != "none" ] 
then
	NODES=`scontrol show hostnames $SLURM_JOB_NODELIST`
	echo $NODES
	for NODE in `echo $NODES | head -n 1`
	do
		echo $NODE
		INTERFERENCE_COMMAND="ssh -f $NODE export OMP_NUM_THREADS=24 OMP_PLACES=cores OMP_PROC_BIND=close; ~/master/swe-benchmark/batch/interference/main $EXTRA"
		echo $INTERFERENCE_COMMAND
		$INTERFERENCE_COMMAND
	done
fi

### Execute your application
echo $COMMAND

for i in $( eval echo {1..$NUM_EXECUTIONS} )
do
	$COMMAND
	echo %
done

if [ "$EXTRA" != "none" ] 
then
	killall ssh
fi
