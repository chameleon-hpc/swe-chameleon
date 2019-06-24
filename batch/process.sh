#!/usr/local_rwth/bin/zsh

#EXPERIMENTS=(chameleon charm++ upcxx mpi)
EXPERIMENTS=(chameleon)
NUM_NODES=(1 2)
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
	for NUM in ${NUM_NODES[@]}
	do
		COMMAND="cat output/swe_${CLUSTER}_${EXP}_${NUM}.txt | grep [RESULT] | sort"
		echo $COMMAND
		$COMMAND
	done
done
