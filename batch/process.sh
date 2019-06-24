#!/bin/bash

#EXPERIMENTS=(chameleon charm++ upcxx mpi)
EXPERIMENTS=(chameleon)
#SIZES=(2048 4096 8192)
SIZES=(2048)
NODE_COUNTS=(1 2)
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
	echo "config,walltime" > ${EXP}.csv
	for NODE_COUNT in ${NODE_COUNTS[@]}
	do
		for SIZE in ${SIZES[@]}
		do
			COMMAND="cat output/swe_${CLUSTER}_${EXP}_${NODE_COUNT}_${SIZE}.txt"
			#echo $COMMAND
			OUTPUT="${CLUSTER}_${EXP}_${NODE_COUNT}_${SIZE},"$($COMMAND | grep RESULT | sort | tail -n 1 | sed 's/[a-zA-Z:, ]//g')
			echo $OUTPUT >> ${EXP}.csv
		done
	done
done
