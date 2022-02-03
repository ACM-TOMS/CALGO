/*********************************************************************/
/*                                                                   */
/* Copyright (C) 2006  Andrey Nikolayevich Chernikov                 */
/*                                                                   */
/* pcdm.h                                                            */
/*                                                                   */
/* The interface for calling PCDM as a library                       */
/*                                                                   */
/*********************************************************************/

#ifndef __PCDM_H__
#define __PCDM_H__

#include "defs.h"
#include "pcdmio.h"

void pcdmInit (int argc, char* argv[]);

void pcdmRefine (struct PcdmIO * pcdmio);

void pcdmFinalize ();

#endif
