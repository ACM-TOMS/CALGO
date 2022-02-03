#include <sys/resource.h>
#include <error.h>
#include <stdlib.h>
#include <stdio.h>

void core_unlimited__(void)
{
    struct rlimit rlim;
    rlim.rlim_cur = RLIM_INFINITY;
    rlim.rlim_max = RLIM_INFINITY;
    if (setrlimit(RLIMIT_CORE, &rlim) == -1 ) {
        perror("setrlimit error");
        exit(1);
    }
}

void print_core_limit__(void)
{
    struct rlimit rlim;
    if (!getrlimit(RLIMIT_CORE, &rlim)) {
        printf("Core limit(soft): %20ld \n", rlim.rlim_cur); 
        printf("Core limit(hard): %20ld \n", rlim.rlim_max);        
    }
}
