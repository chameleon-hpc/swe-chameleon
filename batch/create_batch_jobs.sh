#!/usr/local_rwth/bin/zsh

#EXPERIMENTS=(chameleon charm++ upcxx mpi)
EXPERIMENTS=(mpi chameleon charm++)
SIZES=(4096)
NODE_COUNTS=(1 2 4 8)
NUM_TASKS_PER_NODE=2
# Claix 16 settings
NUM_CPUS_PER_TASK=12
PROJECT=jara0001
PARTITION=c16m
CLUSTER=clx16

# Claix 18 settings
if [ "${CLX18}" = "1" ]; then
	NUM_CPUS_PER_TASK=24
	PROJECT=thes0570
	PARTITION=c18m
	CLUSTER=clx18
fi

for EXP in ${EXPERIMENTS[@]}
do
	for NODE_COUNT in ${NODE_COUNTS[@]}
	do
		for SIZE in ${SIZES[@]}
		do
			# Charm++ settings
			if [ "${EXP}" = "charm" ]; then
				NUM_TASKS_PER_NODE=$(($NUM_TASKS_PER_NODE * $NUM_CPUS_PER_TASK))
				NUM_CPUS_PER_TASK=1
			fi

			export OMP_NUM_THREADS=$NUM_CPUS_PER_TASK
			export SIZE=$SIZE
			export NODE_COUNT=$NODE_COUNT
			export EXP=$EXP
			export CLUSTER=$CLUSTER
			COMMAND="sbatch --partition=$PARTITION --nodes=${NODE_COUNT} --ntasks-per-node=$NUM_TASKS_PER_NODE --cpus-per-task=$NUM_CPUS_PER_TASK --job-name=swe_${CLUSTER}_${EXP}_${NODE_COUNT}_${SIZE} --output=output/swe_${CLUSTER}_${EXP}_${NODE_COUNT}_${SIZE}.txt --export=NODE_COUNT,NUM_TASKS_PER_NODE,OMP_NUM_THREADS,SIZE --account=$PROJECT ${EXP}.sh"
			echo $COMMAND
			$COMMAND
		done
	done
done
