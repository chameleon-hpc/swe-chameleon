#!/bin/bash

NUM_EXECUTIONS=5

#scaling
JOB_CATEGORY=scaling
CLUSTER=clx16
FRAMEWORKS=(mpi chameleon charm++)
NODE_COUNTS=(1 2 4 8)
GRID_DIMENSIONS=(4096)
BLOCK_COUNTS=(16)
DRY_FRACTIONS=(0.0)
EXTRAS=(none)

echo -n "Frameworks/#Nodes" > ${JOB_CATEGORY}.csv
for NODE_COUNT in ${NODE_COUNTS[@]}
do
	for i in $( eval echo {1..$NUM_EXECUTIONS} )
	do
		echo -n ",${NODE_COUNT}" >> ${JOB_CATEGORY}.csv
	done
done
echo "" >> ${JOB_CATEGORY}.csv
for FRAMEWORK in ${FRAMEWORKS[@]}
do
	echo -n "$FRAMEWORK" >> ${JOB_CATEGORY}.csv
	for NODE_COUNT in ${NODE_COUNTS[@]}
	do
		FILE="output/swe_${JOB_CATEGORY}_${CLUSTER}_${FRAMEWORK}_${NODE_COUNT}_4096_16_0.0_none.txt"
		INPUT=$(head -n -1 $FILE)
		# delimiter between experiments
		IFS='%'
				read -d '' -ra RUNS <<< "$INPUT"
		IFS=' '
		for RUN in "${RUNS[@]}"
		do
			OUTPUT=$(echo $RUN | grep RESULT | sort | tail -n 1 | sed 's/[a-zA-Z:, ]//g')
			echo -n ",$OUTPUT" >> ${JOB_CATEGORY}.csv
		done
	done
	echo "" >> ${JOB_CATEGORY}.csv
done

#imbalance
JOB_CATEGORY=imbalance
CLUSTER=clx16
FRAMEWORKS=(mpi chameleon charm++)
NODE_COUNTS=(4)
GRID_DIMENSIONS=(4096)
BLOCK_COUNTS=(16)
DRY_FRACTIONS=(0.0 0.1 0.2 0.3 0.4 0.5)
EXTRAS=(none)

echo -n "Frameworks/#Dry Fraction" > ${JOB_CATEGORY}.csv
for DRY_FRACTION in ${DRY_FRACTIONS[@]}
do
	for i in $( eval echo {1..$NUM_EXECUTIONS} )
	do
		echo -n ",${DRY_FRACTION}" >> ${JOB_CATEGORY}.csv
	done
done
echo "" >> ${JOB_CATEGORY}.csv
for FRAMEWORK in ${FRAMEWORKS[@]}
do
	echo -n "$FRAMEWORK" >> ${JOB_CATEGORY}.csv
	for DRY_FRACTION in ${DRY_FRACTIONS[@]}
	do
		FILE="output/swe_${JOB_CATEGORY}_${CLUSTER}_${FRAMEWORK}_4_4096_16_${DRY_FRACTION}_none.txt"
		INPUT=$(head -n -1 $FILE)
		# delimiter between experiments
		IFS='%'
				read -d '' -ra RUNS <<< "$INPUT"
		IFS=' '
		for RUN in "${RUNS[@]}"
		do
			OUTPUT=$(echo $RUN | grep RESULT | sort | tail -n 1 | sed 's/[a-zA-Z:, ]//g')
			echo -n ",$OUTPUT" >> ${JOB_CATEGORY}.csv
		done
	done
	echo "" >> ${JOB_CATEGORY}.csv
done

#granularity
JOB_CATEGORY=granularity
CLUSTER=clx16
FRAMEWORKS=(chameleon charm++)
NODE_COUNTS=(4)
GRID_DIMENSIONS=(4096)
BLOCK_COUNTS=(8 16 32 64 128)
DRY_FRACTIONS=(0.0)
EXTRAS=(none)

echo -n "Frameworks/#Blocks" > ${JOB_CATEGORY}.csv
for BLOCK_COUNT in ${BLOCK_COUNTS[@]}
do
	for i in $( eval echo {1..$NUM_EXECUTIONS} )
	do
		echo -n ",${BLOCK_COUNT}x${BLOCK_COUNT}" >> ${JOB_CATEGORY}.csv
	done
done
echo "" >> ${JOB_CATEGORY}.csv
for FRAMEWORK in ${FRAMEWORKS[@]}
do
	echo -n "$FRAMEWORK" >> ${JOB_CATEGORY}.csv
	for BLOCK_COUNT in ${BLOCK_COUNTS[@]}
	do
		FILE="output/swe_${JOB_CATEGORY}_${CLUSTER}_${FRAMEWORK}_4_4096_${BLOCK_COUNT}_0.0_none.txt"
		INPUT=$(head -n -1 $FILE)
		# delimiter between experiments
		IFS='%'
				read -d '' -ra RUNS <<< "$INPUT"
		IFS=' '
		for RUN in "${RUNS[@]}"
		do
			OUTPUT=$(echo $RUN | grep RESULT | sort | tail -n 1 | sed 's/[a-zA-Z:, ]//g')
			echo -n ",$OUTPUT" >> ${JOB_CATEGORY}.csv
		done
	done
	echo "" >> ${JOB_CATEGORY}.csv
done

# replace names
sed -i 's/mpi/MPI\/OpenMP/g' *.csv
sed -i 's/charm++/Charm++/g' *.csv
sed -i 's/chameleon/Chameleon/g' *.csv

# plot charts
python3 process.py
