#!/usr/local_rwth/bin/zsh
#SBATCH --exclusive

### Change to the work directory
source /home/sc427635/.zshrc
cd /home/sc427635/master/swe-benchmark

### Load modules
#module load chameleon

### Compile
#CHARM_PATH=/home/sc427635/sw/charm make charm

### Set environment variables
#export OMP_NUM_THREADS=1
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=auto
export I_MPI_DEBUG=5
export OMP_PROC_BIND=close
export DRY_FRACTION=$DRY_FRACTION

### Execute your application
COMMAND="$MPIEXEC $FLAGS_MPI_BATCH ./batch/cpu_set_wrapper.sh ./build/SWE_intel_release_charm_omp_augrie "
COMMAND+="-t 1 -n 1 "  
COMMAND+="-x $GRID_DIMENSION -y $GRID_DIMENSION "
COMMAND+="-i $BLOCK_COUNT -j $BLOCK_COUNT "
COMMAND+="-o ./output/charm_batch -i 200"

echo $COMMAND

if [ $EXTRA != none ] 
	INTERFERENCE_COMMAND="$MPIEXEC $FLAGS_MPI_BATCH ./batch/cpu_set_wrapper.sh ./batch/interference/main $EXTRA"
	echo $INTERFERENCE_COMMAND
	$INTERFERENCE_COMMAND &
	PID=$!
fi

### Execute your application

for i in $( eval echo {1..$NUM_EXECUTIONS} )
do
	$COMMAND
	echo %
done

if [ $EXTRA != none ] 
	echo "Killing $PID"
	kill $PID
fi
