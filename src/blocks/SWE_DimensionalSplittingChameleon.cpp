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
	hvNetUpdatesAbove(nx + 1, ny + 2),
	left(NULL),
	right(NULL) {
		computeTime = 0.;
		computeTimeWall = 0.;

		MPI_Type_vector(nx, 1, ny + 2, MPI_FLOAT, &HORIZONTAL_BOUNDARY);
		MPI_Type_commit(&HORIZONTAL_BOUNDARY);
	}

void SWE_DimensionalSplittingChameleon::setLeft(SWE_DimensionalSplittingChameleon* argLeft) {
	left = argLeft;
}

void SWE_DimensionalSplittingChameleon::setRight(SWE_DimensionalSplittingChameleon* argRight) {
	right = argRight;
}

void SWE_DimensionalSplittingChameleon::setRank(int rank) {
	myRank = rank;
}

void SWE_DimensionalSplittingChameleon::setGhostLayer() {
	if(right != NULL) {
		for(int i = 1; i < ny+1; i++) {
			h[nx+1][i] = left->getWaterHeight()[nx+1][i];
			hu[nx+1][i] = left->getMomentumHorizontal()[nx+1][i];
			hv[nx+1][i] = left->getMomentumVertical()[nx+1][i];
		}
	}
	if(left != NULL) {
		for(int i = 1; i < ny+1; i++) {
			h[0][i] = left->getWaterHeight()[0][i];
			hu[0][i] = left->getMomentumHorizontal()[0][i];
			hv[0][i] = left->getMomentumVertical()[0][i];
		}
	}
	// TODO: Use MPI Communication to set top and bottom GhostLayers
}

void computeNumericalFluxesKernel(SWE_DimensionalSplittingChameleon* block, float* h_data, float* hu_data, float* hv_data, float* b_data,
								float* hNetUpdatesLeft_data, float* hNetUpdatesRight_data, float* huNetUpdatesLeft_data, float* huNetUpdatesRight_data,
								float* hNetUpdatesBelow_data, float* hNetUpdatesAbove_data, float* hvNetUpdatesBelow_data, float* hvNetUpdatesAbove_data,
								float* hStar_data, float* huStar_data) {
	// Set data pointers correctly
	block->getModifiableWaterHeight().setRawPointer(h_data);
	block->getModifiableMomentumHorizontal().setRawPointer(hu_data);
	block->getModifiableMomentumVertical().setRawPointer(hv_data);
	block->getModifiableBathymetry().setRawPointer(b_data);
	block->hNetUpdatesLeft.setRawPointer(hNetUpdatesLeft_data);
	block->hNetUpdatesRight.setRawPointer(hNetUpdatesRight_data);
	block->huNetUpdatesLeft.setRawPointer(huNetUpdatesLeft_data);
	block->huNetUpdatesRight.setRawPointer(huNetUpdatesRight_data);
	block->hNetUpdatesBelow.setRawPointer(hNetUpdatesBelow_data);
	block->hNetUpdatesAbove.setRawPointer(hNetUpdatesAbove_data);
	block->hvNetUpdatesBelow.setRawPointer(hvNetUpdatesBelow_data);
	block->hvNetUpdatesAbove.setRawPointer(hvNetUpdatesAbove_data);
	block->hStar.setRawPointer(hStar_data);
	block->huStar.setRawPointer(huStar_data);
	
	// Start compute clocks
	block->computeClock = clock();
	clock_gettime(CLOCK_MONOTONIC, &(block->startTime));

	//maximum (linearized) wave speed within one iteration
	float maxHorizontalWaveSpeed = (float) 0.;
	float maxVerticalWaveSpeed = (float) 0.;
	solver::Hybrid<float> localSolver = block->solver;

	// x-sweep, compute the actual domain plus ghost rows above and below
	// iterate over cells on the x-axis, leave out the last column (two cells per computation)
	//#pragma omp for reduction(max : maxHorizontalWaveSpeed) collapse(2)
	for (int x = 0; x < block->nx + 1; x++) {
		// iterate over all rows, including ghost layer
		for (int y = 0; y < block->ny + 2; y++) {
			localSolver.computeNetUpdates (
					block->getWaterHeight()[x][y], block->getWaterHeight()[x + 1][y],
					block->getMomentumHorizontal()[x][y], block->getMomentumHorizontal()[x + 1][y],
					block->getBathymetry()[x][y], block->getBathymetry()[x + 1][y],
					block->hNetUpdatesLeft[x][y], block->hNetUpdatesRight[x + 1][y],
					block->huNetUpdatesLeft[x][y], block->huNetUpdatesRight[x + 1][y],
					maxHorizontalWaveSpeed
					);
		}
	}

	// compute max timestep according to cautious CFL-condition
	block->maxTimestep = (float) .4 * (block->dx / maxHorizontalWaveSpeed);

	// set intermediary Q* states
	//#pragma omp for collapse(2)
	for (int x = 1; x < block->nx + 1; x++) {
		for (int y = 0; y < block->ny + 2; y++) {
			block->hStar[x][y] = block->getWaterHeight()[x][y] - (block->maxTimestep / block->dx) * (block->hNetUpdatesLeft[x][y] + block->hNetUpdatesRight[x][y]);
			block->huStar[x][y] = block->getMomentumHorizontal()[x][y] - (block->maxTimestep / block->dx) * (block->huNetUpdatesLeft[x][y] + block->huNetUpdatesRight[x][y]);
		}
	}

	// y-sweep
	//#pragma omp for reduction(max : maxVerticalWaveSpeed) collapse(2)
	for (int x = 1; x < block->nx + 1; x++) {
		for (int y = 0; y < block->ny + 1; y++) {
			localSolver.computeNetUpdates (
					block->getWaterHeight()[x][y], block->getWaterHeight()[x][y + 1],
					block->getMomentumVertical()[x][y], block->getMomentumVertical()[x][y + 1],
					block->getBathymetry()[x][y], block->getBathymetry()[x][y + 1],
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
	block->computeClock = clock() - block->computeClock;
	block->computeTime += (float) block->computeClock / CLOCKS_PER_SEC;

	clock_gettime(CLOCK_MONOTONIC, &(block->endTime));
	block->computeTimeWall += (block->endTime.tv_sec - block->startTime.tv_sec);
	block->computeTimeWall += (float) (block->endTime.tv_nsec - block->startTime.tv_nsec) / 1E9;
}

/**
 * Compute net updates for the block.
 * The member variable #maxTimestep will be updated with the
 * maximum allowed time step size
 */
void SWE_DimensionalSplittingChameleon::computeNumericalFluxes() {

	chameleon_map_data_entry_t* args = new chameleon_map_data_entry_t[15];
    args[0] = chameleon_map_data_entry_create(this, sizeof(SWE_DimensionalSplittingChameleon), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[1] = chameleon_map_data_entry_create(this->getWaterHeight().getRawPointer(), sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[2] = chameleon_map_data_entry_create(this->hu.getRawPointer(), sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[3] = chameleon_map_data_entry_create(this->hv.getRawPointer(), sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[4] = chameleon_map_data_entry_create(this->b.getRawPointer(), sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[5] = chameleon_map_data_entry_create(this->hNetUpdatesLeft.getRawPointer(), sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[6] = chameleon_map_data_entry_create(this->hNetUpdatesRight.getRawPointer(), sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[7] = chameleon_map_data_entry_create(this->huNetUpdatesLeft.getRawPointer(), sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[8] = chameleon_map_data_entry_create(this->huNetUpdatesRight.getRawPointer(), sizeof(float)*(nx + 2)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[9] = chameleon_map_data_entry_create(this->hNetUpdatesBelow.getRawPointer(), sizeof(float)*(nx + 1)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[10] = chameleon_map_data_entry_create(this->hNetUpdatesAbove.getRawPointer(), sizeof(float)*(nx + 1)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[11] = chameleon_map_data_entry_create(this->hvNetUpdatesBelow.getRawPointer(), sizeof(float)*(nx + 1)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[12] = chameleon_map_data_entry_create(this->hvNetUpdatesAbove.getRawPointer(), sizeof(float)*(nx + 1)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[13] = chameleon_map_data_entry_create(this->hStar.getRawPointer(), sizeof(float)*(nx + 1)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
    args[14] = chameleon_map_data_entry_create(this->huStar.getRawPointer(), sizeof(float)*(nx + 1)*(ny + 2), CHAM_OMP_TGT_MAPTYPE_TO | CHAM_OMP_TGT_MAPTYPE_FROM);
	
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
