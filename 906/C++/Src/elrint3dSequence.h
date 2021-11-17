#ifndef  _ELRINT3D_SEQUENCE_H_
#define  _ELRINT3D_SEQUENCE_H_

// Application include
#include    <embeddedSequence.h>

extern const EmbeddedSequence<double>* singleton_elrint3dSequence();
// Return a global singleton instance of the lattice sequence used in 'elrint3d',
// to be shared by all cubature rules based on this sequence. 

#endif  // _ELRINT3D_SEQUENCE_H_ 


