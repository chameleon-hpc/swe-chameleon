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
 * Implementation of SWE_DimensionalSplittingChameleon.hh
 *
 */
#include "SWE_DimensionalSplittingChameleon.hh"

#include <cassert>
#include <algorithm>
#include <omp.h>
#include <mpi.h>
#include "chameleon.h"

/*
 * Constructor of a SWE_DimensionalSplittingChameleon Block.
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
SWE_DimensionalSplittingChameleon::SWE_DimensionalSplittingChameleon (int nx, int ny, float dx, float dy, float originX, float originY) :
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
	hStar (nx + 1, ny + 2),
	huStar (nx + 1, ny + 2),

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

		computeTime = 0.;
		computeTimeWall = 0.;
	}

void SWE_DimensionalSplittingChameleon::setGhostLayer() {
	SWE_Block::applyBoundaryConditions();
}

void computeNumericalFluxesKernel(SWE_DimensionalSplittingChameleon* block, float* h_data, float* hu_data, float* hv_data, float* b_data,
								float* hNetUpdatesLeft_data, float* hNetUpdatesRight_data, float* huNetUpdatesLeft_data, float* huNetUpdatesRight_data,
								float* hNetUpdatesBelow_data, float* hNetUpdatesAbove_data, float* hvNetUpdatesBelow_data, float* hvNetUpdatesAbove_data,
								float* hStar_data, float* huStar_data) {
	// Set data pointers correctly
	block->h.rawData = h_data;
	block->hu.rawData = hu_data;
	block->hv.rawData = hv_data;
	block->b.rawData = b_data;
	block->hNetUpdatesLeft.rawData = hNetUpdatesLeft_data;
	block->hNetUpdatesRight.rawData = hNetUpdatesRight_data_data;
	block->huNetUpdatesLeft.rawData = huNetUpdatesLeft_data_data;
	block->huNetUpdatesRight.rawData = huNetUpdatesRight_data;
	block->hNetUpdatesBelow.rawData = hNetUpdatesBelow_data;
	block->hNetUpdatesAbove.rawData = hNetUpdatesAbove_data;
	block->hvNetUpdatesBelow.rawData = hvNetUpdatesBelow_data;
	block->hvNetUpdatesAbove.rawData = hvNetUpdatesAbove_data;
	block->hStar.rawData = hStar_data;
	block->huStar.rawData = huStar_data;
	
	// Start compute clocks
	computeClock = clock();
	clock_gettime(CLOCK_MONOTONIC, &startTime);

	//maximum (linearized) wave speed within one iteration
	float maxHorizontalWaveSpeed = (float) 0.;
	float maxVerticalWaveSpeed = (float) 0.;
	solver::Hybrid<float> localSolver = block->solver;

	// x-sweep, compute the actual domain plus ghost rows above and below
	// iterate over cells on the x-axis, leave out the last column (two cells per computation)
	#pragma omp for reduction(max : maxHorizontalWaveSpeed) collapse(2)
	for (int x = 0; x < block->nx + 1; x++) {
		// iterate over all rows, including ghost layer
		for (int y = 0; y < block->ny + 2; y++) {
			localSolver.computeNetUpdates (
					block->h[x][y], block->h[x + 1][y],
					block->hu[x][y], block->hu[x + 1][y],
					block->b[x][y], block->b[x + 1][y],
					block->hNetUpdatesLeft[x][y], block->hNetUpdatesRight[x + 1][y],
					block->huNetUpdatesLeft[x][y], block->huNetUpdatesRight[x + 1][y],
					maxHorizontalWaveSpeed
					);
		}
	}

	// compute max timestep according to cautious CFL-condition
	block->maxTimestep = (float) .4 * (block->dx / maxHorizontalWaveSpeed);

	// set intermediary Q* states
	#pragma omp for collapse(2)
	for (int x = 1; x < nx + 1; x++) {
		for (int y = 0; y < ny + 2; y++) {
			block->hStar[x][y] = block->h[x][y] - (block->maxTimestep / dx) * (block->hNetUpdatesLeft[x][y] + block->hNetUpdatesRight[x][y]);
			block->huStar[x][y] = block->hu[x][y] - (block->maxTimestep / dx) * (block->huNetUpdatesLeft[x][y] + block->huNetUpdatesRight[x][y]);
		}
	}

	// y-sweep
	#pragma omp for reduction(max : maxVerticalWaveSpeed) collapse(2)
	for (int x = 1; x < nx + 1; x++) {
		for (int y = 0; y < ny + 1; y++) {
			localSolver.computeNetUpdates (
					block->h[x][y], block->h[x][y + 1],
					block->hv[x][y], block->hv[x][y + 1],
					block->b[x][y], block->b[x][y + 1],
					block->hNetUpdatesBelow[x][y], block->hNetUpdatesAbove[x][y + 1],
					block->hvNetUpdatesBelow[x][y], block->hvNetUpdatesAbove[x][y + 1],
					maxVerticalWaveSpeed
					);
		}
	}
	
	#ifndef NDEBUG
	// check if the cfl condition holds in the y-direction
	assert(block->maxTimestep < (float) .5 * (block->dy / maxVerticalWaveSpeed));
	#endif // NDEBUG

	// Accumulate compute time
	computeClock = clock() - computeClock;
	computeTime += (float) computeClock / CLOCKS_PER_SEC;

	clock_gettime(CLOCK_MONOTONIC, &endTime);
	block->computeTimeWall += (endTime.tv_sec - startTime.tv_sec);
	block->computeTimeWall += (float) (endTime.tv_nsec - startTime.tv_nsec) / 1E9;
}

/**
 * Compute net updates for the block.
 * The member variable #maxTimestep will be updated with the
 * maximum allowed time step size
 */
void SWE_DimensionalSplittingChameleon::computeNumericalFluxes() {

	chameleon_map_data_entry_t* args = new chameleon_map_data_entry_t[15];
    args[0] = chameleon_map_data_entry_create(this, sizeof(SWE_DimensionalSplittingChameleon), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[1] = chameleon_map_data_entry_create(this->h.rawData, sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[2] = chameleon_map_data_entry_create(this->hu.rawData, sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[3] = chameleon_map_data_entry_create(this->hv.rawData, sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[4] = chameleon_map_data_entry_create(this->b.rawData, sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[5] = chameleon_map_data_entry_create(this->hNetUpdatesLeft.rawData, sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[6] = chameleon_map_data_entry_create(this->hNetUpdatesRight.rawData, sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[7] = chameleon_map_data_entry_create(this->huNetUpdatesLeft.rawData, sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[8] = chameleon_map_data_entry_create(this->huNetUpdatesRight.rawData, sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[9] = chameleon_map_data_entry_create(this->hNetUpdatesBelow.rawData, sizeof(float)*(nx + 1)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[10] = chameleon_map_data_entry_create(this->hNetUpdatesAbove.rawData, sizeof(float)*(nx + 1)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[11] = chameleon_map_data_entry_create(this->hvNetUpdatesBelow.rawData, sizeof(float)*(nx + 1)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[12] = chameleon_map_data_entry_create(this->hvNetUpdatesAbove.rawData, sizeof(float)*(nx + 1)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[13] = chameleon_map_data_entry_create(this->hStar.rawData, sizeof(float)*(nx + 1)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[14] = chameleon_map_data_entry_create(this->huStar.rawData, sizeof(float)*(nx + 1)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
	
	int32_t res = chameleon_add_task_manual(
        (void *)&computeNumericalFluxesKernel,
        15, // number of args
        args);
}

/**
 * Updates the unknowns with the already computed net-updates.
 *
 * @param dt time step width used in the update. The timestep has to be equal to maxTimestep calculated by computeNumericalFluxes(),
 * since this is the step width used for the intermediary updates after the x-sweep.
 */
void SWE_DimensionalSplittingChameleon::updateUnknowns (float dt) {
	// Start compute clocks
	computeClock = clock();
	clock_gettime(CLOCK_MONOTONIC, &startTime);

	// this assertion has to hold since the intermediary star states were calculated internally using a timestep width of maxTimestep
	assert(std::abs(dt - maxTimestep) < 0.00001);
	//update cell averages with the net-updates
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
