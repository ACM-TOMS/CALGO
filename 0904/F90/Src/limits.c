#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
//#include <mpi.h>

#define LIMITS 10

typedef struct {
    int resource;
    char *name;
} rlimits;

/* defined in /usr/include/asm-generic/resource.h and/or <bits/resource.h> */
static rlimits limits[] = {
//    { RLIMIT_CPU, "cpu time"},
//    { RLIMIT_FSIZE, "file size"},
//    { RLIMIT_DATA, "data seg size"},
//    { RLIMIT_STACK, "stack size"},
    { RLIMIT_CORE, "core file size"},
//    { RLIMIT_RSS, "resident set size"},
//    { RLIMIT_NPROC, "max user processes"},
//    { RLIMIT_NOFILE, "open files"},
//    { RLIMIT_MEMLOCK, "max locked memory" },
//    { RLIMIT_AS, "max memory size"},
//    { RLIMIT_LOCKS, "File locks"},
//    { RLIMIT_SIGPENDING, "pending signals"},
//    { RLIMIT_MSGQUEUE, "POSIX message queues"},
//    { RLIMIT_NICE, "nice priority"},
//    { RLIMIT_RTPRIO, "real-time priority"},
//    { RLIMIT_RTTIME, "timeout for RT tasks in us"},
//    { RLIMIT_RLIM_NLIMITS, "nlimits"},
};

void print_limits() {
    struct rlimit rlim;
    rlimits *limit;
    for(limit = limits; limit->name; limit++) {
	if (!getrlimit(limit->resource,&rlim)) {
	    printf("%.20s(soft): %20ld \n",limit->name, rlim.rlim_cur);	
	    printf("%.20s(hard): %20ld \n", limit->name, rlim.rlim_max);	
	}
    }
}

/*int main (int argc, char *argv[]) 
{
    struct rlimit rlim;
    rlimits *limit;
    int myid=0, numprocs=0;
    FILE *outfile;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    printf("unlimited=%ld\n",RLIM_INFINITY);	
    printf("First time:\n");
    print_limits();
    // Set core to unlimited
    rlim.rlim_cur = RLIM_INFINITY;
    rlim.rlim_max = RLIM_INFINITY;
    setrlimit(RLIMIT_CORE, &rlim);
    // Print the limits again
    printf("Second time:\n");
    print_limits();
    MPI_Finalize();
    return 0;
    } */
