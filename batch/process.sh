#!/bin/bash

NUM_EXECUTIONS=5

#scaling
JOB_CATEGORY=scaling
CLUSTER=clx16
FRAMEWORKS=(mpi chameleon charm++)
NODE_COUNTS=(1 2 4 8)
GRID_DIMENSIONS=(4096)
BLOCK_COUNTS=(32)
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
		FILE="output/swe_${JOB_CATEGORY}_${CLUSTER}_${FRAMEWORK}_${NODE_COUNT}_4096_32_0.0_none.txt"
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

#interference
JOB_CATEGORY=interference
CLUSTER=clx16
FRAMEWORKS=(mpi chameleon charm++)
NODE_COUNTS=(2)
GRID_DIMENSIONS=(4096)
BLOCK_COUNTS=(32)
DRY_FRACTIONS=(0.2)
EXTRAS=(0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5)

echo -n "Frameworks/Interference" > ${JOB_CATEGORY}.csv
for EXTRA in ${EXTRAS[@]}
do
	for i in $( eval echo {1..$NUM_EXECUTIONS} )
	do
		echo -n ",${EXTRA}" >> ${JOB_CATEGORY}.csv
	done
done
echo "" >> ${JOB_CATEGORY}.csv
for FRAMEWORK in ${FRAMEWORKS[@]}
do
	echo -n "$FRAMEWORK" >> ${JOB_CATEGORY}.csv
	for EXTRA in ${EXTRAS[@]}
	do
		FILE="output/swe_${JOB_CATEGORY}_${CLUSTER}_${FRAMEWORK}_2_4096_32_0.2_${EXTRA}.txt"
		INPUT=$(cat $FILE | grep -v Killed | head -n -1 )
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

#interference
JOB_CATEGORY=interference
CLUSTER=clx16
FRAMEWORKS=(mpi chameleon charm++)
NODE_COUNTS=(2)
GRID_DIMENSIONS=(4096)
BLOCK_COUNTS=(32)
DRY_FRACTIONS=(0.2)
EXTRAS=(1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5)

echo -n "Frameworks/Interference" > single_${JOB_CATEGORY}.csv
for EXTRA in ${EXTRAS[@]}
do
	for i in $( eval echo {1..$NUM_EXECUTIONS} )
	do
		echo -n ",${EXTRA}" >> single_${JOB_CATEGORY}.csv
	done
done
echo "" >> single_${JOB_CATEGORY}.csv
for FRAMEWORK in ${FRAMEWORKS[@]}
do
	echo -n "$FRAMEWORK" >> single_${JOB_CATEGORY}.csv
	for EXTRA in ${EXTRAS[@]}
	do
		FILE="output/swe_${JOB_CATEGORY}_${CLUSTER}_${FRAMEWORK}_2_4096_32_0.2_${EXTRA}.txt"
		INPUT=$(cat $FILE | grep -v Killed | head -n -1 )
		# delimiter between experiments
		IFS='%'
				read -d '' -ra RUNS <<< "$INPUT"
		IFS=' '
		for RUN in "${RUNS[@]}"
		do
			OUTPUT=$(echo $RUN | grep RESULT | sort | tail -n 1 | sed 's/[a-zA-Z:, ]//g')
			echo -n ",$OUTPUT" >> single_${JOB_CATEGORY}.csv
		done
	done
	echo "" >> single_${JOB_CATEGORY}.csv
done

#imbalance
JOB_CATEGORY=imbalance
CLUSTER=clx16
FRAMEWORKS=(mpi chameleon charm++)
NODE_COUNTS=(2)
GRID_DIMENSIONS=(4096)
BLOCK_COUNTS=(32)
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
		FILE="output/swe_${JOB_CATEGORY}_${CLUSTER}_${FRAMEWORK}_2_4096_32_${DRY_FRACTION}_none.txt"
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
NODE_COUNTS=(2)
GRID_DIMENSIONS=(4096)
BLOCK_COUNTS=(128 64 32 16 8)
DRY_FRACTIONS=(0.2)
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
		FILE="output/swe_${JOB_CATEGORY}_${CLUSTER}_${FRAMEWORK}_2_4096_${BLOCK_COUNT}_0.2_none.txt"
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
sed -i 's/mpi/MPI+OpenMP (24 Threads per Node)/g' *.csv
sed -i 's/charm++/Charm++ (24 Processors per Node)/g' *.csv
sed -i 's/chameleon/Chameleon (22 Threads per Node)/g' *.csv
sed -i 's/starpu/StarPU (24 Threads per Node)/g' *.csv

# plot charts
python3 process.py
