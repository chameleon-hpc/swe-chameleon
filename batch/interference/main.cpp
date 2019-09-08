#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <thread>
#include <chrono>

double getTime() {
    struct timespec time;
    clock_gettime(CLOCK_MONOTONIC, &time);
    return (double) time.tv_sec + ((double)time.tv_nsec)/1E9;
}

int main(int argc, char **argv)
{
	float time_fraction = atof(argv[1]);
	bool all_cores = true;
	if(time_fraction > 1.0) {
		all_cores = false;
		time_fraction -= 1.0;
	}
	if(all_cores)
		printf("all cores, time_fraction=%f\n, ", time_fraction);
	else
		printf("one core, time_fraction=%f\n", time_fraction);

    while(true) {
		printf("step\n");
    	#pragma omp parallel
    	{
    		double start_time = getTime();
    		//if(omp_get_thread_num() == 0)
			//	printf("start time = %lf\n", start_time);
    		if(all_cores || omp_get_thread_num() == 0) {
	    		while(getTime() - start_time < 0.1*time_fraction) {
		    			int vec_len = 1 << 10;
		    			int iteration_count = 1 << 10;
						double vector[vec_len];
						for(int i = 0; i < vec_len; i++)
							vector[i] = 1.0;
						for(int j = 0; j < iteration_count; j++)
							for(int i = 0; i < vec_len; i++)
								vector[i] = vector[i]+1.000001 + 1;
		    	}
	    	}
	    	#pragma omp barrier
    		//if(omp_get_thread_num() == 0)
			//	printf("wait start time = %lf\n", getTime()-start_time);
	    	int num_milliseconds = (int)100*(1-time_fraction);
	    	std::this_thread::sleep_for(std::chrono::milliseconds(num_milliseconds));
    		//if(omp_get_thread_num() == 0)
			//	printf("wait end time = %lf\n", getTime()-start_time);
    	}
    }
    return 0;
}
