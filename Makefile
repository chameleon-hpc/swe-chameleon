#run_smp:
#	./build/SWE_gnu_release_none_omp_augrie -t 3600 -n 20 -x 1000 -y 1000 -o ~/storage/tsunami/simulation/tohu_1m_new -b ~/master/data/tohoku/bath.nc -d ~/master/data/tohoku/displ.nc

run_upcxx:
	${UPCXX_PATH}/bin/upcxx-run -n 4 ./build/SWE_gnu_release_upcxx_augrie -t 3600 -n 20 -x 1000 -y 1000 -o ~/storage/tsunami/simulation/upcxx -b ~/master/data/tohoku/bath.nc -d ~/master/data/tohoku/displ.nc
run_upcxx_asagi:
	${UPCXX_PATH}/bin/upcxx-run -n 4 ./build/SWE_gnu_release_upcxx_augrie -t 3600 -n 20 -x 1000 -y 1000 -o ~/storage/tsunami/simulation/upcxx -b ~/master/data/tohoku/bath.nc -d ~/master/data/tohoku/displ.nc
run_upcxx_test:
	${UPCXX_PATH}/bin/upcxx-run -n 4 ./build/SWE_gnu_release_upcxx_augrie -t 60 -n 10 -x 10 -y 10 -o ~/storage/tsunami/simulation/radial_upcxx

run_mpi_single:

run_mpi:
	I_MPI_PIN=1 I_MPI_PIN_DOMAIN=auto OMP_NUM_THREADS=12 OMP_PLACES=cores OMP_PROC_BIND=close mpiexec.hydra -np 2 ./build/SWE_intel_release_mpi_omp_augrie -t 0.1 -n 1 -x 4096 -y 4096 -o ./output/test -i 200 
run_mpi_asagi:
	I_MPI_PIN=1 I_MPI_PIN_DOMAIN=auto OMP_NUM_THREADS=12 OMP_PLACES=cores OMP_PROC_BIND=close mpiexec.hydra -np 2 ./build/SWE_intel_release_mpi_omp_augrie -t 0.1 -n 1 -x 4096 -y 4096 -o ./output/test -u 1 -v 1 -b ~/master/data/tohoku/bath.nc -d ~/master/data/tohoku/displ.nc

#run_ampi:
#	./charmrun +p4 ./build/SWE_gnu_release_mpi_augrie -t 3600 -n 20 -x 1000 -y 1000 -o ~/storage/tsunami/simulation/mpi -b ~/master/data/tohoku/bath.nc -d ~/master/data/tohoku/displ.nc

run_charm:
	I_MPI_PIN=1 I_MPI_PIN_DOMAIN=auto mpiexec.hydra -np 24 ./build/SWE_intel_release_charm_omp_augrie -t 0.1 -n 1 -x 4096 -y 4096 -o ./output/test
run_charm_asagi:
	I_MPI_PIN=1 I_MPI_PIN_DOMAIN=auto mpiexec.hydra -np 24 ./build/SWE_intel_release_charm_omp_augrie -t 0.1 -n 1 -x 4096 -y 4096 -o ./output/test -b ~/master/data/tohoku/bath.nc -d ~/master/data/tohoku/displ.nc
run_charm_test:
	I_MPI_PIN=1 I_MPI_PIN_DOMAIN=auto mpiexec.hydra -np 24 ./build/SWE_intel_release_charm_omp_augrie -t 0.01 -n 1 -x 320 -y 320 -o ./output/test -b ~/master/data/tohoku/bath.nc -d ~/master/data/tohoku/displ.nc

run_chameleon:
	I_MPI_PIN=1 I_MPI_PIN_DOMAIN=auto OMP_NUM_THREADS=11 OMP_PLACES=cores OMP_PROC_BIND=close mpiexec.hydra -np 2 ./build/SWE_intel_release_chameleon_omp_augrie -t 0.1 -n 1 -x 4096 -y 4096 -o ./output/test -i 200
run_chameleon_asagi:
	I_MPI_PIN=1 I_MPI_PIN_DOMAIN=auto OMP_NUM_THREADS=11 OMP_PLACES=cores OMP_PROC_BIND=close mpiexec.hydra -np 2 ./build/SWE_intel_release_chameleon_omp_augrie -t 0.1 -n 1 -x 4096 -y 4096 -o ./output/test -u 1 -v 1 -b ~/master/data/tohoku/bath.nc -d ~/master/data/tohoku/displ.nc
run_chameleon_test:
	I_MPI_PIN=1 I_MPI_PIN_DOMAIN=auto OMP_NUM_THREADS=11 OMP_PLACES=cores OMP_PROC_BIND=close mpiexec.hydra -np 2 ./build/SWE_intel_release_chameleon_omp_augrie -t 0.01 -n 1 -x 320 -y 320 -o ./output/test

#debug_charm_test:
#	/home/jurek/repository/tum/ccs_tools/bin/charmdebug +p4 ./build/SWE_gnu_release_charm_augrie -t 60 -n 10 -x 10 -y 10 -o ~/storage/tsunami/simulation/charm

#ampi:
#	scons writeNetCDF=True openmp=False parallelization=ampi

smp:
	scons writeNetCDF=True openmp=True parallelization=none compiler=intel
smp_asagi:
	scons writeNetCDF=True openmp=True parallelization=none compiler=intel asagi=true asagiDir=${ASAGI_PATH}

mpi:
	scons writeNetCDF=True openmp=True parallelization=mpi compiler=intel
mpi_asagi:
	scons writeNetCDF=True openmp=True parallelization=mpi compiler=intel asagi=true asagiDir=${ASAGI_PATH}
#mpi_noopenmp:
#	scons writeNetCDF=True openmp=false parallelization=mpi asagi=true asagiDir=${ASAGI_PATH}

upcxx:
	scons writeNetCDF=True openmp=True parallelization=upcxx compiler=intel
upcxx_asagi:
	scons writeNetCDF=True openmp=True parallelization=upcxx compiler=intel asagi=true asagiDir=${ASAGI_PATH}
#upcxx_noopenmp:
#	scons openmp=false parallelization=upcxx compiler=intel

charm:
	scons openmp=true parallelization=charm compiler=intel
charm_asagi:
	scons openmp=true parallelization=charm compiler=intel asagi=true asagiDir=${ASAGI_PATH}

chameleon:
	scons writeNetCDF=True compiler=intel openmp=true parallelization=chameleon
chameleon_asagi:
	scons writeNetCDF=True compiler=intel openmp=true parallelization=chameleon asagi=true asagiDir=${ASAGI_PATH}
chameleon_debug:
	scons writeNetCDF=True compiler=intel openmp=true parallelization=chameleon compileMode=debug


clean:
	rm -rf ./build/SWE_*
	rm -rf ./build/build_*
