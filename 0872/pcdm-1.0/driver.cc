/*********************************************************************/
/*                                                                   */
/* Andrey Chernikov, January 2006                                    */
/*                                                                   */
/* driver.cc                                                         */
/*                                                                   */
/* This file is an example on how to use PCDM as a library           */
/*                                                                   */
/*********************************************************************/

#include "pcdm.h"
#include "pcdmio.h"
#ifdef MEM_CHECK
#include <mcheck.h>
#endif

int main (int argc, char* argv[])
{
    PcdmIO      pcdmio;

#ifdef MEM_CHECK
    mtrace();
#endif

    // ** load the data **
    
    pcdmio.oneFlag = 1;
    pcdmio.read("./models/pipe4000_madd.dat", "./models/pipe4000_madd.cdt");
    pcdmio.areaBound = 0.5;
    pcdmio.polyDir = "./poly";
    pcdmio.statDir = "./stat";

    // ** call the PCDM library **
    
    pcdmInit(argc, argv);
    
    pcdmRefine(&pcdmio);
    
    pcdmFinalize();

    return 0;
}
