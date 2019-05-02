/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *         Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * Basic setting of SWE, which uses a wave propagation solver and an artificial or ASAGI scenario on a single block.
 */

#include <mpi.h>
#include "chameleon.h"
#include <cassert>
#include <string>
#include <ctime>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <algorithm>

#include "tools/args.hh"

#ifdef WRITENETCDF
#include "writer/NetCdfWriter.hh"
#else
#include "writer/VtkWriter.hh"
#endif

#ifdef ASAGI
#include "scenarios/SWE_AsagiScenario.hh"
#else
#include "scenarios/SWE_simple_scenarios.hh"
#endif

#include "blocks/SWE_DimensionalSplittingChameleon.hh"

int main(int argc, char** argv) {



	/**************
	 * INIT INPUT *
	 **************/


	// Define command line arguments
	tools::Args args;

#ifdef ASAGI
	args.addOption("bathymetry-file", 'b', "File containing the bathymetry");
	args.addOption("displacement-file", 'd', "File containing the displacement");
#endif
	args.addOption("simulation-duration", 't', "Time in seconds to simulate");
	args.addOption("checkpoint-count", 'n', "Number of simulation snapshots to be written");
	args.addOption("resolution-horizontal", 'x', "Number of simulation cells in horizontal direction");
	args.addOption("resolution-vertical", 'y', "Number of simulated cells in y-direction");
	args.addOption("output-basepath", 'o', "Output base file name");


	// Declare the variables needed to hold command line input
	float simulationDuration;
	int numberOfCheckPoints;
	int nxRequested;
	int nyRequested;
	std::string outputBaseName;

	// Declare variables for the output and the simulation time
	std::string outputFileName;
	float t = 0.;

	// Parse command line arguments
	tools::Args::Result ret = args.parse(argc, argv);
	switch (ret)
	{
		case tools::Args::Error:
			return 1;
		case tools::Args::Help:
			return 0;
		case tools::Args::Success:
			break;
	}

	// Read in command line arguments
	simulationDuration = args.getArgument<float>("simulation-duration");
	numberOfCheckPoints = args.getArgument<int>("checkpoint-count");
	nxRequested = args.getArgument<int>("resolution-horizontal");
	nyRequested = args.getArgument<int>("resolution-vertical");
	outputBaseName = args.getArgument<std::string>("output-basepath");

	// Initialize Scenario
#ifdef ASAGI
	SWE_AsagiScenario scenario(args.getArgument<std::string>("bathymetry-file"), args.getArgument<std::string>("displacement-file"));
#else
	SWE_RadialDamBreakScenario scenario;
#endif

	// Compute when (w.r.t. to the simulation time in seconds) the checkpoints are reached
	float* checkpointInstantOfTime = new float[numberOfCheckPoints];
	// Time delta is the time between any two checkpoints
	float checkpointTimeDelta = simulationDuration / numberOfCheckPoints;
	// The first checkpoint is reached after 0 + delta t
	checkpointInstantOfTime[0] = checkpointTimeDelta;
	for(int i = 1; i < numberOfCheckPoints; i++) {
		checkpointInstantOfTime[i] = checkpointInstantOfTime[i - 1] + checkpointTimeDelta;
	}


	/***************
	 * INIT BLOCKS *
	 ***************/

	int myRank, numRanks;
	int provided;
	int requested = MPI_THREAD_MULTIPLE;
	MPI_Init_thread(&argc, &argv, requested, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	// Print status
	char hostname[HOST_NAME_MAX];
        gethostname(hostname, HOST_NAME_MAX);

	printf("%i Spawned at %s\n", myRank, hostname);

	/*
	 * Calculate the simulation grid layout.
	 * The cell count of the scenario as well as the scenario size is fixed, 
	 * Get the size of the actual domain and divide it by the requested resolution.
	 */
	int widthScenario = scenario.getBoundaryPos(BND_RIGHT) - scenario.getBoundaryPos(BND_LEFT);
	int heightScenario = scenario.getBoundaryPos(BND_TOP) - scenario.getBoundaryPos(BND_BOTTOM);
	float dxSimulation = (float) widthScenario / nxRequested;
	float dySimulation = (float) heightScenario / nyRequested;

	// hardcode just for testing
	int num_blocks_per_rank = 1;
	SWE_DimensionalSplittingChameleon* blocks[num_blocks_per_rank];

	int x_blocksize = nxRequested / num_blocks_per_rank;
	// y values are the same for all blocks on the same rank
	int y_blocksize = nyRequested / numRanks;
	int y_pos = std::min(myRank*y_blocksize, nyRequested);
	float originY = y_pos * dySimulation;
	int block_ny = std::min(y_blocksize, nyRequested - y_pos);
		
	for(int i = 0; i < num_blocks_per_rank; i++) {
		// TODO: extract into function in chameleonSplitting
		// divide horizontal slice vertically into blocks
		int x_pos = std::min(i*x_blocksize, nxRequested);
		float originX = x_pos * dxSimulation;
		int block_nx = std::min(x_blocksize, nxRequested - x_pos);

		BoundaryType boundaries[4];

		if(i == 0)
			boundaries[BND_LEFT] = scenario.getBoundaryType(BND_LEFT);
		else
			boundaries[BND_LEFT] = CONNECT;
		if(i == num_blocks_per_rank - 1)
			boundaries[BND_RIGHT] = scenario.getBoundaryType(BND_RIGHT);
		else
			boundaries[BND_RIGHT] = CONNECT;
		if(myRank == 0)
			boundaries[BND_BOTTOM] = scenario.getBoundaryType(BND_BOTTOM);
		else
			boundaries[BND_BOTTOM] = CONNECT;
		if(myRank == numRanks - 1)
			boundaries[BND_TOP] = scenario.getBoundaryType(BND_TOP);
		else
			boundaries[BND_TOP] = CONNECT;

		blocks[i] = new SWE_DimensionalSplittingChameleon(block_nx, block_ny, dxSimulation, dySimulation, originX, originY);
		blocks[i]->initScenario(scenario, boundaries);

		blocks[i]->myRank = myRank;
		if(myRank != 0)
			blocks[i]->neighbourRankId[BND_BOTTOM] = myRank-1;
		if(myRank != numRanks - 1)
			blocks[i]->neighbourRankId[BND_TOP] = myRank+1;

		//printf("Init block %d with originX:%d and originY:%d\n", i, blocks[i]->getOriginX(), blocks[i]->getOriginY());
	}

	for(int i = 0; i < num_blocks_per_rank; i++) {
		if(i != 0)
			blocks[i]->left = blocks[i-1];
		if(i != num_blocks_per_rank-1)
			blocks[i]->right = blocks[i+1];
	}


	/***************
	 * INIT OUTPUT *
	 ***************/

	// block used for writing (only used on rank 0)
	// all ranks write their blocks to this write block on rank 0 (using one-sided communication)
	// This block is then written resulting in a single file
	SWE_DimensionalSplittingChameleon writeBlock(nxRequested, nyRequested, dxSimulation, dySimulation, 0, 0);
	BoundaryType boundaries[4];
	boundaries[BND_LEFT] = scenario.getBoundaryType(BND_LEFT);
	boundaries[BND_RIGHT] = scenario.getBoundaryType(BND_RIGHT);
	boundaries[BND_TOP] = scenario.getBoundaryType(BND_TOP);
	boundaries[BND_BOTTOM] = scenario.getBoundaryType(BND_BOTTOM);
	writeBlock.initScenario(scenario, boundaries);
  	MPI_Win writeBlockWin_h;
	MPI_Win writeBlockWin_hu;
  	MPI_Win writeBlockWin_hv;

	MPI_Win_create(writeBlock.h.getRawPointer(), (nxRequested+2)*(nyRequested*2), sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &writeBlockWin_h);
	MPI_Win_create(writeBlock.hu.getRawPointer(), (nxRequested+2)*(nyRequested*2), sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &writeBlockWin_hu);
	MPI_Win_create(writeBlock.hv.getRawPointer(), (nxRequested+2)*(nyRequested*2), sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &writeBlockWin_hv);

	// Initialize boundary size of the ghost layers
	BoundarySize boundarySize = {{1, 1, 1, 1}};
#ifdef WRITENETCDF
	// Construct netCDF writers
	outputFileName = outputBaseName;
	NetCdfWriter writer(
		outputFileName,
		writeBlock.getBathymetry(),
		boundarySize,
		writeBlock.getCellCountHorizontal(),
		writeBlock.getCellCountVertical(),
		dxSimulation,
		dySimulation,
		writeBlock.getOriginX(),
		writeBlock.getOriginY());
	//printf("Init writer %d with originX:%d and originY:%d\n", i, blocks[i]->getOriginX(), blocks[i]->getOriginY());
#else
	// Construct a vtk writer
	outputFileName = outputBaseName;
	VtkWriter writer(
		outputFileName,
		writeBlock.getBathymetry(),
		boundarySize,
		writeBlock.getCellCountHorizontal(),
		writeBlock.getCellCountVertical(),
		dxSimulation,
		dySimulation);
#endif // WRITENETCDF

	// Write the output at t = 0
	writer.writeTimeStep(
		writeBlock.getWaterHeight(),
		writeBlock.getMomentumHorizontal(),
		writeBlock.getMomentumVertical(),
		(float) 0.);

	/********************
	 * START SIMULATION *
	 ********************/

    // chameleon_init();
    #pragma omp parallel
    {
        chameleon_thread_init();
    }
	
    // necessary to be aware of binary base addresses to calculate offset for target functions
    chameleon_determine_base_addresses((void *)&main);

    MPI_Barrier(MPI_COMM_WORLD);

	// Initialize timers
	std::clock_t computeClock;
	std::clock_t commClock;

	struct timespec startTime;
	struct timespec endTime;

	float computeTime = 0.;
	float commTime = 0.;
	float wallTime = 0.;

	t = 0.0;

	float timestep;
	unsigned int iterations = 0;

	// loop over the count of requested checkpoints
	for(int i = 0; i < numberOfCheckPoints; i++) {
		// Simulate until the checkpoint is reached
		while(t < checkpointInstantOfTime[i]) {

			// Start measurement
			clock_gettime(CLOCK_MONOTONIC, &startTime);
			commClock = clock();

			timestep = std::numeric_limits<float>::max();

			#pragma omp parallel for
			for(int i=0; i<num_blocks_per_rank; i++) {

				// set values in ghost cells.
				// we need to sync here since block boundaries get exchanged over ranks
				blocks[i]->setGhostLayer();
			}
			//chameleon_distributed_taskwait(0);

			// Accumulate comm time and start compute clock
			commClock = clock() - commClock;
			commTime += (float) commClock / CLOCKS_PER_SEC;
			computeClock = clock();

			#pragma omp parallel for
			for(int i=0; i<num_blocks_per_rank; i++) {
				// compute numerical flux on each edge
				blocks[i]->computeNumericalFluxes();
			}
			chameleon_distributed_taskwait(0);
			
			//#pragma omp parallel for
			for(int i=0; i<num_blocks_per_rank; i++) {
				// reduce timestep
				if(blocks[i]->getMaxTimestep() < timestep)
					timestep = blocks[i]->getMaxTimestep();
			}

			// reduce over all ranks
			float maxTimestepGlobal;
			MPI_Allreduce(&timestep, &maxTimestepGlobal, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
			timestep = maxTimestepGlobal;

			for(int i=0; i<num_blocks_per_rank; i++) {
				// update the cell values
				blocks[i]->maxTimestep = timestep;
				blocks[i]->updateUnknowns(timestep);
			}
			//chameleon_distributed_taskwait(0);

			// Accumulate compute time
			computeClock = clock() - computeClock;
			computeTime += (float) computeClock / CLOCKS_PER_SEC;

			// Accumulate wall time
			clock_gettime(CLOCK_MONOTONIC, &endTime);
			wallTime += (endTime.tv_sec - startTime.tv_sec);
			wallTime += (float) (endTime.tv_nsec - startTime.tv_nsec) / 1E9;
			// update simulation time with time step width.
			t += timestep;
			iterations++;

			MPI_Barrier(MPI_COMM_WORLD);
			printf("%d: Step, current time:%f\n", myRank, t);
		}

		printf("Write timestep (%fs)\n", t);

		// write output
		MPI_Win_fence(0, writeBlockWin_h);
		MPI_Win_fence(0, writeBlockWin_hu);
		MPI_Win_fence(0, writeBlockWin_hv);
		for(int i=0; i<num_blocks_per_rank; i++) {
			// Send all data to rank 0, which will write it to a single file
			// send each column separately
			for(int j=0; j<blocks[i]->nx; j++) {
				int x_pos = std::min(i*x_blocksize, nxRequested);
				int y_pos = std::min(i*y_blocksize, nyRequested);
				MPI_Put(blocks[i]->h.getRawPointer()+1+(blocks[i]->ny+2)*j, blocks[i]->ny, MPI_FLOAT,
					0, 1+(nyRequested+2)*(1+j+x_pos)+y_pos, blocks[i]->ny, MPI_FLOAT, writeBlockWin_h);
				MPI_Put(blocks[i]->hu.getRawPointer()+1+(blocks[i]->ny+2)*j, blocks[i]->ny, MPI_FLOAT,
					0, 1+(nyRequested+2)*(1+j+x_pos)+y_pos, blocks[i]->ny, MPI_FLOAT, writeBlockWin_hu);
				MPI_Put(blocks[i]->hv.getRawPointer()+1+(blocks[i]->ny+2)*j, blocks[i]->ny, MPI_FLOAT,
					0, 1+(nyRequested+2)*(1+j+x_pos)+y_pos, blocks[i]->ny, MPI_FLOAT, writeBlockWin_hv);
			}
		}
		MPI_Win_fence(0, writeBlockWin_h);
		MPI_Win_fence(0, writeBlockWin_hu);
		MPI_Win_fence(0, writeBlockWin_hv);

		if(myRank == 0) {
			writer.writeTimeStep(
				writeBlock.getWaterHeight(),
				writeBlock.getMomentumHorizontal(),
				writeBlock.getMomentumVertical(),
				t);
		}
	}


	/************
	 * FINALIZE *
	 ************/


	printf("SMP : Compute Time (CPU): %fs - (WALL): %fs | Total Time (Wall): %fs\n", blocks[0]->computeTime, blocks[0]->computeTimeWall, wallTime); 
	for(int i=0; i<num_blocks_per_rank; i++) {
		blocks[i]->freeMpiType();
	}
	MPI_Finalize();

	return 0;
}
