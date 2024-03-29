/**
 * @file
 * This file is part of an SWE fork created for the Tsunami-Simulation Bachelor Lab Course.
 *
 * @author Jurek Olden (jurek.olden AT in.tum.de)
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
 * Implementation of SWE_DimensionalSplitting.hh
 *
 */
#include "SWE_DimensionalSplittingMpi.hh"

#include <cassert>
#include <algorithm>
#include <omp.h>

/*
 * Constructor of a SWE_DimensionalSplitting Block.
 * Computational domain is [1,...,nx]*[1,...,ny]
 * Ghost layer consists of two additional rows and columns
 *
 * State variables h, hu, hv and b are defined on the whole grid (including ghost layer)
 * Net updates coming from above/below/left/right are defined for each cell.
 *
 * Net updates are computed on all rows first, then on all columns, the total net updates are then composed
 * from the two 1D solutions.
 *
 * This strategy only works, if the timestep chosen w.r.t. to the maximum horizontal wave speeds
 * also satisfies the CFL-condition in y-direction.
 *
 * @param l_nx Size of the computational domain in x-direction
 * @param l_ny Size of the computational domain in y-direction
 * @param l_dx Cell width
 * @param l_dy Cell height
 */
SWE_DimensionalSplittingMpi::SWE_DimensionalSplittingMpi(int nx, int ny, float dx, float dy, float originX, float originY) :
		/*
		 * Important note concerning grid allocations:
		 * Since index shifts all over the place are bug-prone and maintenance unfriendly,
		 * an index of [x][y] is at the actual position x,y on the actual grid.
		 * This implies that the allocation size in any direction might be larger than the number of values needed.
		 * So if, for instance, array[x][y] needs to hold values in the domain [1,a][1,b],
		 * it will be allocated with size (a+1, b+1) instead of (a, b).
		 * array[0][0] is then unused.
		 */
		// Initialize grid metadata using the base class constructor
		SWE_Block(nx, ny, dx, dy, originX, originY),

		// intermediate state Q after x-sweep
		hStar(nx + 1, ny + 2),
		huStar(nx + 1, ny + 2),

		/*
		 * Temporary storage for the net updates per grid cell during a sweep.
		 * There are four update values per cell:
		 * Left-going wave from the right edge, analogue for the left edge.
		 * Down-going wave from the top edge, analogue for the bottom edge
		 */

		// For the x-sweep
		hNetUpdatesLeft(nx + 2, ny + 2),
		hNetUpdatesRight(nx + 2, ny + 2),

		huNetUpdatesLeft(nx + 2, ny + 2),
		huNetUpdatesRight(nx + 2, ny + 2),

		// For the y-sweep
		hNetUpdatesBelow(nx + 1, ny + 2),
		hNetUpdatesAbove(nx + 1, ny + 2),

		hvNetUpdatesBelow(nx + 1, ny + 2),
		hvNetUpdatesAbove(nx + 1, ny + 2) {

	MPI_Type_vector(nx, 1, ny + 2, MPI_FLOAT, &HORIZONTAL_BOUNDARY);
	MPI_Type_commit(&HORIZONTAL_BOUNDARY);

	computeTime = 0.;
	computeTimeWall = 0.;
}

void SWE_DimensionalSplittingMpi::freeMpiType() {
	MPI_Type_free(&HORIZONTAL_BOUNDARY);
}

void SWE_DimensionalSplittingMpi::connectNeighbours(int p_neighbourRankId[]) {
	for (int i = 0; i < 4; i++) {
		neighbourRankId[i] = p_neighbourRankId[i];
	}
}

void SWE_DimensionalSplittingMpi::exchangeBathymetry() {

	/*********
	 * SEND *
	 ********/

	// The requests generated by the Isends are immediately freed, since we will wait on the requests generated by the corresponding receives
	MPI_Request req;

	if (boundaryType[BND_LEFT] == CONNECT) {
		int startIndex = ny + 2 + 1;
		MPI_Isend(b.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_LEFT], MPI_TAG_OUT_B_LEFT, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);
	}
	if (boundaryType[BND_RIGHT] == CONNECT) {
		int startIndex = nx * (ny + 2) + 1;
		MPI_Isend(b.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_RIGHT], MPI_TAG_OUT_B_RIGHT, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);
	}
	if (boundaryType[BND_BOTTOM] == CONNECT) {
		for (int i = 1; i < nx + 1; i++) {
			MPI_Isend(&b[i][1], 1, MPI_FLOAT, neighbourRankId[BND_BOTTOM], MPI_TAG_OUT_B_BOTTOM, MPI_COMM_WORLD, &req); 
			MPI_Request_free(&req);
		}
	}
	if (boundaryType[BND_TOP] == CONNECT) {
		for (int i = 1; i < nx + 1; i++) {
			MPI_Isend(&b[i][ny], 1, MPI_FLOAT, neighbourRankId[BND_TOP], MPI_TAG_OUT_B_TOP, MPI_COMM_WORLD, &req); 
			MPI_Request_free(&req);
		}
	}

	/***********
	 * RECEIVE *
	 **********/

	MPI_Request recvReqs[4];
	MPI_Status stati[4];

	if (boundaryType[BND_LEFT] == CONNECT) {
		int startIndex = 1;
		MPI_Irecv(b.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_LEFT], MPI_TAG_OUT_B_RIGHT, MPI_COMM_WORLD, &recvReqs[BND_LEFT]);
	} else {
		recvReqs[BND_LEFT] = MPI_REQUEST_NULL;
	}

	if (boundaryType[BND_RIGHT] == CONNECT) {
		int startIndex = (nx + 1) * (ny + 2) + 1;
		MPI_Irecv(b.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_RIGHT], MPI_TAG_OUT_B_LEFT, MPI_COMM_WORLD, &recvReqs[BND_RIGHT]);
	} else {
		recvReqs[BND_RIGHT] = MPI_REQUEST_NULL;
	}

	if (boundaryType[BND_BOTTOM] == CONNECT) {
		for (int i = 1; i < nx + 1; i++) {
			MPI_Irecv(&b[i][0], 1, MPI_FLOAT, neighbourRankId[BND_BOTTOM], MPI_TAG_OUT_B_TOP, MPI_COMM_WORLD, &recvReqs[BND_BOTTOM]); 
		}
	} else {
		recvReqs[BND_BOTTOM] = MPI_REQUEST_NULL;
	}

	if (boundaryType[BND_TOP] == CONNECT) {
		for (int i = 1; i < nx + 1; i++) {
			MPI_Irecv(&b[i][ny + 1], 1, MPI_FLOAT, neighbourRankId[BND_TOP], MPI_TAG_OUT_B_BOTTOM, MPI_COMM_WORLD, &recvReqs[BND_TOP]); 
		}
	} else {
		recvReqs[BND_TOP] = MPI_REQUEST_NULL;
	}

	MPI_Waitall(4, recvReqs, stati);
}

void SWE_DimensionalSplittingMpi::setGhostLayer() {
	// Apply appropriate conditions for OUTFLOW/WALL boundaries
	SWE_Block::applyBoundaryConditions();

	MPI_Status status;

	assert(h.getRows() == ny + 2);
	assert(hu.getRows() == ny + 2);
	assert(hv.getRows() == ny + 2);
	assert(h.getCols() == nx + 2);
	assert(hu.getCols() == nx + 2);
	assert(hv.getCols() == nx + 2);

	/*********
	 * SEND *
	 ********/

	// The requests generated by the Isends are immediately freed, since we will wait on the requests generated by the corresponding receives
	MPI_Request req;

	if (boundaryType[BND_LEFT] == CONNECT) {
		int startIndex = ny + 2 + 1;

		MPI_Isend(h.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_LEFT], MPI_TAG_OUT_H_LEFT, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);

		MPI_Isend(hu.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_LEFT], MPI_TAG_OUT_HU_LEFT, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);

		MPI_Isend(hv.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_LEFT], MPI_TAG_OUT_HV_LEFT, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);
	}
	if (boundaryType[BND_RIGHT] == CONNECT) {
		int startIndex = nx * (ny + 2) + 1;

		MPI_Isend(h.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_RIGHT], MPI_TAG_OUT_H_RIGHT, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);

		MPI_Isend(hu.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_RIGHT], MPI_TAG_OUT_HU_RIGHT, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);

		MPI_Isend(hv.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_RIGHT], MPI_TAG_OUT_HV_RIGHT, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);
	}
	if (boundaryType[BND_BOTTOM] == CONNECT) {

		MPI_Isend(&h[1][1], 1, HORIZONTAL_BOUNDARY, neighbourRankId[BND_BOTTOM], MPI_TAG_OUT_H_BOTTOM, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);

		MPI_Isend(&hu[1][1], 1, HORIZONTAL_BOUNDARY, neighbourRankId[BND_BOTTOM], MPI_TAG_OUT_HU_BOTTOM, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);

		MPI_Isend(&hv[1][1], 1, HORIZONTAL_BOUNDARY, neighbourRankId[BND_BOTTOM], MPI_TAG_OUT_HV_BOTTOM, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);

	}
	if (boundaryType[BND_TOP] == CONNECT) {

		MPI_Isend(&h[1][ny], 1, HORIZONTAL_BOUNDARY, neighbourRankId[BND_TOP], MPI_TAG_OUT_H_TOP, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);

		MPI_Isend(&hu[1][ny], 1, HORIZONTAL_BOUNDARY, neighbourRankId[BND_TOP], MPI_TAG_OUT_HU_TOP, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);

		MPI_Isend(&hv[1][ny], 1, HORIZONTAL_BOUNDARY, neighbourRankId[BND_TOP], MPI_TAG_OUT_HV_TOP, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);

	}

	/***********
	 * RECEIVE *
	 **********/

	// 4 Boundaries times 3 arrays (h, hu, hv) means 12 requests
	// The requests corresponding to the h array will be at indices 0, 3, 6, 9
	// hu array requests will be at indices 1, 4, 5, 10 and so forth
	MPI_Request recvReqs[12];
	MPI_Status stati[12];

	if (boundaryType[BND_LEFT] == CONNECT) {
		int startIndex = 1;
		MPI_Irecv(h.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_LEFT], MPI_TAG_OUT_H_RIGHT, MPI_COMM_WORLD, &recvReqs[0]);
		MPI_Irecv(hu.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_LEFT], MPI_TAG_OUT_HU_RIGHT, MPI_COMM_WORLD, &recvReqs[1]);
		MPI_Irecv(hv.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_LEFT], MPI_TAG_OUT_HV_RIGHT, MPI_COMM_WORLD, &recvReqs[2]);
	} else {
		recvReqs[0] = MPI_REQUEST_NULL;
		recvReqs[1] = MPI_REQUEST_NULL;
		recvReqs[2] = MPI_REQUEST_NULL;
	}

	if (boundaryType[BND_RIGHT] == CONNECT) {
		int startIndex = (nx + 1) * (ny + 2) + 1;
		MPI_Irecv(h.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_RIGHT], MPI_TAG_OUT_H_LEFT, MPI_COMM_WORLD, &recvReqs[3]);
		MPI_Irecv(hu.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_RIGHT], MPI_TAG_OUT_HU_LEFT, MPI_COMM_WORLD, &recvReqs[4]);
		MPI_Irecv(hv.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourRankId[BND_RIGHT], MPI_TAG_OUT_HV_LEFT, MPI_COMM_WORLD, &recvReqs[5]);
	} else {
		recvReqs[3] = MPI_REQUEST_NULL;
		recvReqs[4] = MPI_REQUEST_NULL;
		recvReqs[5] = MPI_REQUEST_NULL;
	}

	if (boundaryType[BND_BOTTOM] == CONNECT) {

		MPI_Irecv(&h[1][0], 1, HORIZONTAL_BOUNDARY, neighbourRankId[BND_BOTTOM], MPI_TAG_OUT_H_TOP, MPI_COMM_WORLD, &recvReqs[6]); 
		MPI_Irecv(&hu[1][0], 1, HORIZONTAL_BOUNDARY, neighbourRankId[BND_BOTTOM], MPI_TAG_OUT_HU_TOP, MPI_COMM_WORLD, &recvReqs[7]); 
		MPI_Irecv(&hv[1][0], 1, HORIZONTAL_BOUNDARY, neighbourRankId[BND_BOTTOM], MPI_TAG_OUT_HV_TOP, MPI_COMM_WORLD, &recvReqs[8]); 

	} else {
		recvReqs[6] = MPI_REQUEST_NULL;
		recvReqs[7] = MPI_REQUEST_NULL;
		recvReqs[8] = MPI_REQUEST_NULL;
	}
	
	if (boundaryType[BND_TOP] == CONNECT) {

		MPI_Irecv(&h[1][ny + 1], 1, HORIZONTAL_BOUNDARY, neighbourRankId[BND_TOP], MPI_TAG_OUT_H_BOTTOM, MPI_COMM_WORLD, &recvReqs[9]); 
		MPI_Irecv(&hu[1][ny + 1], 1, HORIZONTAL_BOUNDARY, neighbourRankId[BND_TOP], MPI_TAG_OUT_HU_BOTTOM, MPI_COMM_WORLD, &recvReqs[10]); 
		MPI_Irecv(&hv[1][ny + 1], 1, HORIZONTAL_BOUNDARY, neighbourRankId[BND_TOP], MPI_TAG_OUT_HV_BOTTOM, MPI_COMM_WORLD, &recvReqs[11]); 

	} else {
		recvReqs[9] = MPI_REQUEST_NULL;
		recvReqs[10] = MPI_REQUEST_NULL;
		recvReqs[11] = MPI_REQUEST_NULL;
	}
	

	MPI_Waitall(12, recvReqs, stati);
}

/**
 * Compute net updates for the block.
 * The member variable #maxTimestep will be updated with the
 * maximum allowed time step size
 */
void SWE_DimensionalSplittingMpi::computeNumericalFluxes () {
	// Start compute clocks
	computeClock = clock();
	clock_gettime(CLOCK_MONOTONIC, &startTime);

	//maximum (linearized) wave speed within one iteration
	float maxHorizontalWaveSpeed = (float) 0.;
	float maxVerticalWaveSpeed = (float) 0.;

	#pragma omp parallel private(solver)
	{
		// x-sweep, compute the actual domain plus ghost rows above and below
		// iterate over cells on the x-axis, leave out the last column (two cells per computation)
		#pragma omp for reduction(max : maxHorizontalWaveSpeed)
		for (int x = 0; x < nx + 1; x++) {
			// iterate over all rows, including ghost layer
			for (int y = 0; y < ny + 2; y++) {
				solver.computeNetUpdates (
						h[x][y], h[x + 1][y],
						hu[x][y], hu[x + 1][y],
						b[x][y], b[x + 1][y],
						hNetUpdatesLeft[x][y], hNetUpdatesRight[x + 1][y],
						huNetUpdatesLeft[x][y], huNetUpdatesRight[x + 1][y],
						maxHorizontalWaveSpeed
						);
			}
		}
	}

	// Accumulate compute time -> exclude the reduction
	computeClock = clock() - computeClock;
	computeTime += (float) computeClock / CLOCKS_PER_SEC;

	clock_gettime(CLOCK_MONOTONIC, &endTime);
	computeTimeWall += (endTime.tv_sec - startTime.tv_sec);
	computeTimeWall += (float) (endTime.tv_nsec - startTime.tv_nsec) / 1E9;

	// compute max timestep according to cautious CFL-condition
	maxTimestep = (float) .4 * (dx / maxHorizontalWaveSpeed);
	MPI_Allreduce(&maxTimestep, &maxTimestepGlobal, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
	maxTimestep = maxTimestepGlobal;

	// restart compute clocks
	computeClock = clock();
	clock_gettime(CLOCK_MONOTONIC, &startTime);

	#pragma omp parallel private(solver)
	{
		// set intermediary Q* states
		#pragma omp for
		for (int x = 1; x < nx + 1; x++) {
			for (int y = 0; y < ny + 2; y++) {
				hStar[x][y] = h[x][y] - (maxTimestep / dx) * (hNetUpdatesLeft[x][y] + hNetUpdatesRight[x][y]);
				huStar[x][y] = hu[x][y] - (maxTimestep / dx) * (huNetUpdatesLeft[x][y] + huNetUpdatesRight[x][y]);
			}
		}

		// y-sweep
		#ifndef NDEBUG
		#pragma omp for
		#else
		#pragma omp for reduction(max : maxVerticalWaveSpeed)
		#endif
		for (int x = 1; x < nx + 1; x++) {
			for (int y = 0; y < ny + 1; y++) {
				solver.computeNetUpdates (
						h[x][y], h[x][y + 1],
						hv[x][y], hv[x][y + 1],
						b[x][y], b[x][y + 1],
						hNetUpdatesBelow[x][y], hNetUpdatesAbove[x][y + 1],
						hvNetUpdatesBelow[x][y], hvNetUpdatesAbove[x][y + 1],
						maxVerticalWaveSpeed
						);
			}
		}

		#ifndef NDEBUG
		#pragma omp single
		{
			// check if the cfl condition holds in the y-direction
			assert(maxTimestep < (float) .5 * (dy / maxVerticalWaveSpeed));
		}
		#endif // NDEBUG
	}

	// Accumulate compute time
	computeClock = clock() - computeClock;
	computeTime += (float) computeClock / CLOCKS_PER_SEC;
	
	clock_gettime(CLOCK_MONOTONIC, &endTime);
	computeTimeWall += (endTime.tv_sec - startTime.tv_sec);
	computeTimeWall += (float) (endTime.tv_nsec - startTime.tv_nsec) / 1E9;
}

/**
 * Updates the unknowns with the already computed net-updates.
 *
 * @param dt time step width used in the update. The timestep has to be equal to maxTimestep calculated by computeNumericalFluxes(),
 * since this is the step width used for the intermediary updates after the x-sweep.
 */
void SWE_DimensionalSplittingMpi::updateUnknowns (float dt) {
	// Start compute clocks
	computeClock = clock();
	clock_gettime(CLOCK_MONOTONIC, &startTime);

	// this assertion has to hold since the intermediary star states were calculated internally using a timestep width of maxTimestep
	assert(std::abs(dt - maxTimestep) < 0.00001);
	//update cell averages with the net-updates
	#pragma omp parallel for
	for (int x = 1; x < nx+1; x++) {
		for (int y = 1; y < ny + 1; y++) {
			h[x][y] = hStar[x][y] - (maxTimestep / dx) * (hNetUpdatesBelow[x][y] + hNetUpdatesAbove[x][y]);
			hu[x][y] = huStar[x][y];
			hv[x][y] = hv[x][y] - (maxTimestep / dx) * (hvNetUpdatesBelow[x][y] + hvNetUpdatesAbove[x][y]);
		}
	}

	// Accumulate compute time
	computeClock = clock() - computeClock;
	computeTime += (float) computeClock / CLOCKS_PER_SEC;
	
	clock_gettime(CLOCK_MONOTONIC, &endTime);
	computeTimeWall += (endTime.tv_sec - startTime.tv_sec);
	computeTimeWall += (float) (endTime.tv_nsec - startTime.tv_nsec) / 1E9;
}
