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

	//printf("%i Spawned at %s\n", myRank, hostname);

	/*
	 * Calculate the simulation grid layout.
	 * The cell count of the scenario as well as the scenario size is fixed, 
	 * Get the size of the actual domain and divide it by the requested resolution.
	 */
	int widthScenario = scenario.getBoundaryPos(BND_RIGHT) - scenario.getBoundaryPos(BND_LEFT);
	int heightScenario = scenario.getBoundaryPos(BND_TOP) - scenario.getBoundaryPos(BND_BOTTOM);
	float dxSimulation = (float) widthScenario / nxRequested;
	float dySimulation = (float) heightScenario / nyRequested;

	int xRankCount = 1;
	int yRankCount = numRanks;
	int xBlockCount = 2;
	int yBlockCount = 2;
	//int num_blocks_per_rank = 4;
	float xWeights[xRankCount];
	float yWeights[yRankCount];
	float xSum = 0.0;
	float ySum = 0.0;
	// read and normalize weights
	for(int i = 0; i < xRankCount; i++){
		xWeights[i] = 1.0;
		xSum += xWeights[i];
	}
	for(int i = 0; i < yRankCount; i++) {
		yWeights[i] = 1.0;
		ySum +=yWeights[i];
	}

	for(int i = 0; i < xRankCount; i++)
		xWeights[i] = xWeights[i] / xSum * xBlockCount;
	for(int i = 0; i < yRankCount; i++)
		yWeights[i] = yWeights[i] / ySum * yBlockCount;

	int xBounds[xRankCount+1];
	xBounds[0] = 0;
	int yBounds[yRankCount+1];
	yBounds[0] = 0;
	for(int i = 1; i < xRankCount; i++)
		xBounds[i] = xBounds[i-1] + xWeights[i];
	xBounds[xRankCount] = xBlockCount;
	for(int i = 1; i < yRankCount; i++)
		yBounds[i] = yBounds[i-1] + yWeights[i];
	yBounds[yRankCount] = yBlockCount;
	//printf("%d: xBlockCount:%d\n", myRank, xBlockCount);
	//printf("%d: xBounds:%d, %d\n", myRank, xBounds[0], xBounds[1]);
	//printf("%d: yBounds:%d, %d, %d\n", myRank, yBounds[0], yBounds[1], yBounds[2]);
	//printf("%d: xDim:%d, yDim:%d\n", myRank, xBounds[(myRank%xRankCount)+1]-xBounds[myRank%xRankCount], yBounds[(myRank/xRankCount)+1]-yBounds[myRank/xRankCount]);

	SWE_DimensionalSplittingChameleon* blocks[xBlockCount][yBlockCount];

	int x_blocksize = nxRequested / xBlockCount;
	int y_blocksize = nyRequested / yBlockCount;

	//printf("%d: loop bounds: %d, %d, %d, %d\n", myRank, xBounds[myRank%xRankCount], xBounds[(myRank%xRankCount)+1], yBounds[myRank/xRankCount], yBounds[(myRank/xRankCount)+1]);
	for(int x = xBounds[myRank%xRankCount]; x < xBounds[(myRank%xRankCount)+1]; x++) {
		for(int y = yBounds[myRank/xRankCount]; y < yBounds[(myRank/xRankCount)+1]; y++) {
			//printf("%d: x=%d, y=%d\n", myRank, x, y);
			int x_pos = x*x_blocksize;
			float originX = x_pos * dxSimulation;
			int block_nx = x_blocksize;
			if(x == xBlockCount-1)
				block_nx += nxRequested % x_blocksize;
			int y_pos = y*y_blocksize;
			float originY = y_pos * dySimulation;
			int block_ny = y_blocksize;
			if(y == yBlockCount-1)
				block_ny += nyRequested % y_blocksize;

			BoundaryType boundaries[4];			

			if(x == 0)
				boundaries[BND_LEFT] = scenario.getBoundaryType(BND_LEFT);
			else if(x == xBounds[myRank%xRankCount])
				boundaries[BND_LEFT] = CONNECT;
			else
				boundaries[BND_LEFT] = CONNECT_WITHIN_RANK;
			
			if(x_pos+block_nx == nxRequested)
				boundaries[BND_RIGHT] = scenario.getBoundaryType(BND_RIGHT);
			else if(x == xBounds[(myRank%xRankCount)+1]-1)
				boundaries[BND_RIGHT] = CONNECT;
			else
				boundaries[BND_RIGHT] = CONNECT_WITHIN_RANK;

			if(y == 0)
				boundaries[BND_BOTTOM] = scenario.getBoundaryType(BND_BOTTOM);
			else if(y == yBounds[myRank/xRankCount])
				boundaries[BND_BOTTOM] = CONNECT;
			else
				boundaries[BND_BOTTOM] = CONNECT_WITHIN_RANK;

			//printf("%d:First condition: %d==%d\n", myRank, y_pos + block_ny, nyRequested-1);
			//printf("%d:Second condition: %d==%d\n", myRank, y, yBounds[myRank+1]-1);
			if(y_pos + block_ny == nyRequested)
				boundaries[BND_TOP] = scenario.getBoundaryType(BND_TOP);
			else if(y == yBounds[(myRank/yRankCount)+1]-1)
				boundaries[BND_TOP] = CONNECT;
			else
				boundaries[BND_TOP] = CONNECT_WITHIN_RANK;

			blocks[x][y] = new SWE_DimensionalSplittingChameleon(block_nx, block_ny, dxSimulation, dySimulation, originX, originY);
			blocks[x][y]->initScenario(scenario, boundaries);

			blocks[x][y]->myRank = myRank;
			if(myRank != 0)
				blocks[x][y]->neighbourRankId[BND_LEFT] = myRank-1;
			if(myRank != numRanks - 1)
				blocks[x][y]->neighbourRankId[BND_RIGHT] = myRank+1;
			if(myRank >= xRankCount)
				blocks[x][y]->neighbourRankId[BND_BOTTOM] = myRank-xRankCount;
			if(myRank < numRanks - xRankCount)
				blocks[x][y]->neighbourRankId[BND_TOP] = myRank+xRankCount;
		}
	}
	for(int x = xBounds[myRank%xRankCount]; x < xBounds[(myRank%xRankCount)+1]; x++) {
		for(int y = yBounds[myRank/xRankCount]; y < yBounds[(myRank/xRankCount)+1]; y++) {
			if(x != 0)
				blocks[x][y]->left = blocks[x-1][y];
			if(x != xBounds[myRank+1]-1)
				blocks[x][y]->right = blocks[x+1][y];
			if(y != 0)
				blocks[x][y]->bottom = blocks[x][y-1];
			if(y != yBounds[myRank+1]-1)
				blocks[x][y]->top = blocks[x][y+1];
			//printf("%d: Init blocks[%d,%d] with  block_nx:%d, block_ny:%d, originX:%d and originY:%d\n", myRank, x, y, blocks[x][y]->nx, blocks[x][y]->nx, blocks[x][y]->getOriginX(), blocks[x][y]->getOriginY());
		}
	}

	/***************
	 * INIT OUTPUT *
	 ***************/

	// block used for writing (only used on rank 0)
	// all ranks write their blocks to this write block on rank 0 (using one-sided communication)
	// This block is then written to get a single output file
	//printf("%d: Init write block with nxReq:%d, nyReq:%d, dxSim:%f, dySim:%f\n", myRank, nxRequested, nyRequested, dxSimulation, dySimulation);
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
	// Construct a netCDF writer
	outputFileName = outputBaseName;
	NetCdfWriter* writer;
	if(myRank == 0) {
		writer = new NetCdfWriter(
			outputFileName,
			writeBlock.getBathymetry(),
			boundarySize,
			writeBlock.getCellCountHorizontal(),
			writeBlock.getCellCountVertical(),
			dxSimulation,
			dySimulation,
			writeBlock.getOriginX(),
			writeBlock.getOriginY());
		//printf("%d: Init writer with CellCountHorizontal:%d, CellCountVertical:%d, OriginX:%d, getOriginX:%d\n", myRank, writeBlock.getCellCountHorizontal(), writeBlock.getCellCountVertical(), writeBlock.getOriginX(), writeBlock.getOriginY());
	}
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
	if(myRank == 0) {
		writer->writeTimeStep(
			writeBlock.getWaterHeight(),
			writeBlock.getMomentumHorizontal(),
			writeBlock.getMomentumVertical(),
			(float) 0.);
	}

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

			//TODO: exchange bathymetry

			//#pragma omp parallel for
			for(int x = xBounds[myRank%xRankCount]; x < xBounds[(myRank%xRankCount)+1]; x++) {
				for(int y = yBounds[myRank/xRankCount]; y < yBounds[(myRank/xRankCount)+1]; y++) {
					//printf("%d: x=%d, y=%d\n", myRank, x, y);
					// TODO: first send all MPI messages, then receive all
					// set values in ghost cells.
					// we need to sync here since block boundaries get exchanged over ranks
					blocks[x][y]->setGhostLayer();
				}
			}
			//#pragma omp parallel for
			for(int x = xBounds[myRank%xRankCount]; x < xBounds[(myRank%xRankCount)+1]; x++) {
				for(int y = yBounds[myRank/xRankCount]; y < yBounds[(myRank/xRankCount)+1]; y++) {
					blocks[x][y]->receiveGhostLayer();
				}
			}
			//chameleon_distributed_taskwait(0);

			// Accumulate comm time and start compute clock
			commClock = clock() - commClock;
			commTime += (float) commClock / CLOCKS_PER_SEC;
			computeClock = clock();

			#pragma omp parallel for
			for(int x = xBounds[myRank%xRankCount]; x < xBounds[(myRank%xRankCount)+1]; x++) {
				for(int y = yBounds[myRank/xRankCount]; y < yBounds[(myRank/xRankCount)+1]; y++) {
					//printf("%d: x=%d, y=%d\n", myRank, x, y);
					// compute numerical flux on each edge
					blocks[x][y]->computeNumericalFluxesHorizontal();
				}
			}
			chameleon_distributed_taskwait(0);

			for(int x = xBounds[myRank%xRankCount]; x < xBounds[(myRank%xRankCount)+1]; x++) {
				for(int y = yBounds[myRank/xRankCount]; y < yBounds[(myRank/xRankCount)+1]; y++) {
					//printf("%d: x=%d, y=%d\n", myRank, x, y);
					if(blocks[x][y]->getMaxTimestep() < timestep)
						timestep = blocks[x][y]->getMaxTimestep();
				}
			}

			// reduce over all ranks
			float maxTimestepGlobal;
			MPI_Allreduce(&timestep, &maxTimestepGlobal, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
			timestep = maxTimestepGlobal;

			#pragma omp parallel for
			for(int x = xBounds[myRank%xRankCount]; x < xBounds[(myRank%xRankCount)+1]; x++) {
				for(int y = yBounds[myRank/xRankCount]; y < yBounds[(myRank/xRankCount)+1]; y++) {
					//printf("%d: x=%d, y=%d\n", myRank, x, y);
					// compute numerical flux on each edge
					blocks[x][y]->maxTimestep = timestep;
					blocks[x][y]->computeNumericalFluxesVertical();
				}
			}
			chameleon_distributed_taskwait(0);
			
			for(int x = xBounds[myRank%xRankCount]; x < xBounds[(myRank%xRankCount)+1]; x++) {
				for(int y = yBounds[myRank/xRankCount]; y < yBounds[(myRank/xRankCount)+1]; y++) {
					// update the cell values
					blocks[x][y]->updateUnknowns(timestep);
				}
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
			//printf("%d: Step, current time:%f\n", myRank, t);
		}

		//printf("Write timestep to rank 0 (%fs)\n", t);

		// write output
		MPI_Win_fence(0, writeBlockWin_h);
		MPI_Win_fence(0, writeBlockWin_hu);
		MPI_Win_fence(0, writeBlockWin_hv);
		for(int x = xBounds[myRank%xRankCount]; x < xBounds[(myRank%xRankCount)+1]; x++) {
			for(int y = yBounds[myRank/xRankCount]; y < yBounds[(myRank/xRankCount)+1]; y++) {
				//printf("%d: x=%d, y=%d\n", myRank, x, y);
				// Send all data to rank 0, which will write it to a single file
				// send each column separately
				for(int j=1; j<blocks[x][y]->nx+1; j++) {
					int x_pos = x*x_blocksize;
					int y_pos = y*y_blocksize;
					MPI_Put(blocks[x][y]->h.getRawPointer()+1+(blocks[x][y]->ny+2)*j, blocks[x][y]->ny, MPI_FLOAT,
						0, 1+(nyRequested+2)*(1+j+x_pos)+y_pos, blocks[x][y]->ny, MPI_FLOAT, writeBlockWin_h);
					MPI_Put(blocks[x][y]->hu.getRawPointer()+1+(blocks[x][y]->ny+2)*j, blocks[x][y]->ny, MPI_FLOAT,
						0, 1+(nyRequested+2)*(1+j+x_pos)+y_pos, blocks[x][y]->ny, MPI_FLOAT, writeBlockWin_hu);
					MPI_Put(blocks[x][y]->hv.getRawPointer()+1+(blocks[x][y]->ny+2)*j, blocks[x][y]->ny, MPI_FLOAT,
						0, 1+(nyRequested+2)*(1+j+x_pos)+y_pos, blocks[x][y]->ny, MPI_FLOAT, writeBlockWin_hv);
				}
			}
		}

		MPI_Win_fence(0, writeBlockWin_h);
		MPI_Win_fence(0, writeBlockWin_hu);
		MPI_Win_fence(0, writeBlockWin_hv);

		if(myRank == 0) {
			//printf("Write timestep (%fs)\n", t);
			writer->writeTimeStep(
				writeBlock.getWaterHeight(),
				writeBlock.getMomentumHorizontal(),
				writeBlock.getMomentumVertical(),
				t);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}


	/************
	 * FINALIZE *
	 ************/


	//printf("SMP : Compute Time (CPU): %fs - (WALL): %fs | Total Time (Wall): %fs\n", blocks[xBounds[myRank%xRankCount]][myRank/xRankCount]->computeTime, blocks[xBounds[myRank%xRankCount]][myRank/xRankCount]->computeTimeWall, wallTime);
	printf("Chameleon: Computation ended\n");
	MPI_Barrier(MPI_COMM_WORLD);

    #pragma omp parallel
    {
        chameleon_thread_finalize();
    }
    chameleon_finalize();

	for(int x = xBounds[myRank%xRankCount]; x < xBounds[(myRank%xRankCount)+1]; x++) {
		for(int y = yBounds[myRank/xRankCount]; y < yBounds[(myRank/xRankCount)+1]; y++) {
			blocks[x][y]->freeMpiType();
		}
	}

	MPI_Finalize();

	return 0;
}
