#!/bin/bash

create_batch_jobs() {
NUM_TASKS_PER_NODE=2
NUM_EXECUTIONS=5

# Claix 16 settings
NUM_CPUS_PER_TASK=12
PROJECT=jara0001
PARTITION=c16m
CLUSTER=clx16

# Claix 18 settings
if [ "${CLX18}" = "1" ]; then
  NUM_CPUS_PER_TASK=24
  PROJECT=thes0573
  PARTITION=c18m
  CLUSTER=clx18
fi

for FRAMEWORK in ${FRAMEWORKS[@]}
do
  for NODE_COUNT in ${NODE_COUNTS[@]}
  do
    for GRID_DIMENSION in ${GRID_DIMENSIONS[@]}
    do
      for BLOCK_COUNT in ${BLOCK_COUNTS[@]}
      do
        for DRY_FRACTION in ${DRY_FRACTIONS[@]}
        do
          for EXTRA in ${EXTRAS[@]}
          do
            # Charm++ settings
            if [ "${FRAMEWORK}" = "charm++" ]; then
              NUM_TASKS_PER_NODE=$(($NUM_TASKS_PER_NODE * $NUM_CPUS_PER_TASK))
              NUM_CPUS_PER_TASK=1
            fi
            # MPI settings
            if [ "${FRAMEWORK}" = "mpi" ]; then
              NUM_CPUS_PER_TASK=$(($NUM_CPUS_PER_TASK-1))
            fi
            # Chameleon settings
            if [ "${FRAMEWORK}" = "mpi" ]; then
              NUM_CPUS_PER_TASK=$(($NUM_CPUS_PER_TASK-1))
            fi

            export NUM_EXECUTIONS=$NUM_EXECUTIONS
            export CLUSTER=$CLUSTER
            export FRAMEWORK=$FRAMEWORK
            export NODE_COUNT=$NODE_COUNT
            export OMP_NUM_THREADS=$NUM_CPUS_PER_TASK
            export GRID_DIMENSION=$GRID_DIMENSION
            export BLOCK_COUNT=$BLOCK_COUNT
            export DRY_FRACTION=$DRY_FRACTION
            COMMAND="sbatch "
            COMMAND+="--partition=$PARTITION "
            COMMAND+="--nodes=${NODE_COUNT} "
            COMMAND+="--ntasks-per-node=$NUM_TASKS_PER_NODE "
            COMMAND+="--cpus-per-task=$NUM_CPUS_PER_TASK "
            COMMAND+="--job-name=swe_${JOB_CATEGORY}_${CLUSTER}_${FRAMEWORK}_${NODE_COUNT}_${GRID_DIMENSION}_${BLOCK_COUNT}_${DRY_FRACTION}_${EXTRA} "
            COMMAND+="--output=output/swe_${JOB_CATEGORY}_${CLUSTER}_${FRAMEWORK}_${NODE_COUNT}_${GRID_DIMENSION}_${BLOCK_COUNT}_${DRY_FRACTION}_${EXTRA}.txt "
            COMMAND+="--export=NUM_EXECUTIONS,NODE_COUNT,NUM_TASKS_PER_NODE,OMP_NUM_THREADS,GRID_DIMENSION,BLOCK_COUNT,DRY_FRACTION "
            COMMAND+="--account=$PROJECT "
            COMMAND+="${FRAMEWORK}.sh"
            echo $COMMAND
            $COMMAND
          done
        done
      done
    done
  done
done

JOB_CATEGORY=default
FRAMEWORKS=(mpi)
NODE_COUNTS=(2)
GRID_DIMENSIONS=(4096)
BLOCK_COUNTS=(32)
DRY_FRACTIONS=(0.0)
EXTRAS=(none)
}

JOB_CATEGORY=default
FRAMEWORKS=(mpi)
NODE_COUNTS=(2)
GRID_DIMENSIONS=(4096)
BLOCK_COUNTS=(32)
DRY_FRACTIONS=(0.0)
EXTRAS=(none)

create_scaling_jobs() {
  JOB_CATEGORY=scaling
  FRAMEWORKS=(mpi chameleon charm++)
  NODE_COUNTS=(1 2 4 8)
  create_batch_jobs
}

create_imbalance_jobs() {
  JOB_CATEGORY=imbalance
  FRAMEWORKS=(mpi chameleon charm++)
  NODE_COUNTS=(2)
  DRY_FRACTIONS=(0.0 0.1 0.2 0.3 0.4 0.5)
  create_batch_jobs
}

create_granularity_jobs() {
  JOB_CATEGORY=granularity
  FRAMEWORKS=(chameleon charm++)
  NODE_COUNTS=(2)
  DRY_FRACTIONS=(0.2)
  BLOCK_COUNTS=(8 16 32 64 128)
  create_batch_jobs
}

create_scaling_jobs
create_imbalance_jobs
create_granularity_jobs
