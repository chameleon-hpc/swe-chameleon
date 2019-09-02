#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include <thread>
#include <chrono>
#include <random>

volatile bool run = true;
thread_local std::mt19937 gen(1024);
thread_local std::uniform_int_distribution<> dis(1, 6);

int power(int base, int exp, int limit) {
	int res = 1;
	for(int i = 0; i < exp; i++) {
		res *= base;
		res = res % limit;
		//printf("res=%d\n", res);
	}
	return res;
}

int pseudo_random(int limit, int seed) {
	std::uniform_int_distribution<> dis(0, limit-1);
	return dis(gen);
	/*
	int low_seed = seed%8 + 2;
	int high_seed = seed%29 + 3;
	int first = power(low_seed, high_seed, 1000000000);
	printf("first = %d\n", first);
	//printf("%d\n", power(power((25+seed)%30, 2+(seed%8))%30 + (seed%13), (5+seed)%8)%limit);
	printf("res = %d\n", first%limit);
	return first%limit;
	*/
}

bool is_active(int round, int num_total_interference_threads, int rank, int numRanks, int thread, int numThreads) {
	for(int i = 0; i < num_total_interference_threads; i++) {
		int random_rank = pseudo_random(numRanks, round*i);
		int random_thread = pseudo_random(numThreads, round*i+1);
		if(rank == random_rank && thread == random_thread)
			return true;
	}
	return false;
}

void kernel() {
	int vec_len = 1 << 16;
	double vector[vec_len];
	for(int i = 0; i < vec_len; i++)
		vector[i] = 1.0;
	while(run) {
		for(int i = 0; i < vec_len; i++)
			vector[i] = vector[i]+1.000001 + 1;
	}
	//printf("Exit kernel\n");
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

	int num_total_interference_threads = atoi(argv[1]);
	printf("num_total_interference_threads=%d\n", num_total_interference_threads);

    int myRank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    int round = 0;
    while(true) {
    	run = true;
    	#pragma omp parallel
    	{
	    	//if(is_active(round, num_total_interference_threads, myRank, numRanks, omp_get_thread_num(), omp_get_num_threads()) && omp_get_thread_num() != 0) {
	    	if(is_active(round, num_total_interference_threads, myRank, numRanks, omp_get_thread_num(), omp_get_num_threads())) {
	    		//printf("active on rank %d, thread%d\n", myRank, omp_get_thread_num());
	    		kernel();
	    	}
	    	else
		    	std::this_thread::sleep_for(std::chrono::milliseconds(50));

			//printf("%d, %d:Set run to false\n", myRank, omp_get_thread_num());
		    run = false;
		    #pragma omp barrier
    	}

    	round++;
    }

    MPI_Finalize();
    return 0;
}
 
