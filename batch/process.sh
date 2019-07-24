#!/bin/bash

#FRAMEWORKS=(chameleon charm++ upcxx mpi)
FRAMEWORKS=(mpi charm++ chameleon)
#SIZES=(2048 4096 8192)
SIZES=(4096)
NODE_COUNTS=(2)
IMBALANCES=(0.0 0.1 0.2 0.3 0.4 0.5)
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

for IMBALANCE in ${IMBALANCES[@]}
do
	echo "config,walltime" > ${IMBALANCE}.csv
	for EXP in ${FRAMEWORKS[@]}
	do
		for NODE_COUNT in ${NODE_COUNTS[@]}
		do
			for SIZE in ${SIZES[@]}
			do
				COMMAND="cat output/swe_${CLUSTER}_${EXP}_${NODE_COUNT}_${SIZE}_${IMBALANCE}.txt"
				#echo $COMMAND
				OUTPUT="${CLUSTER}_${EXP}_${NODE_COUNT}_${SIZE}_${IMBALANCE},"$($COMMAND | grep RESULT | sort | tail -n 1 | sed 's/[a-zA-Z:, ]//g')
				echo $OUTPUT >> ${IMBALANCE}.csv
			done
		done
	done
done
