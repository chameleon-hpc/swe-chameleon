#!/bin/bash
 
#COMMAND="$MPIEXEC $FLAGS_MPI_BATCH ./batch/cpu_set_wrapper.sh ./interference/interference &"
export OMP_NUM_THREADS=2
COMMAND="mpiexec.hydra -np 2 ./interference/main 1"

$COMMAND &

PID=$!

sleep 100

echo $PID
kill $PID
