#ifndef _PRECOMPUTED_COORDS_H_
#define	_PRECOMPUTED_COORDS_H_

extern void get_coords_by_index(double * const coords, int index);
// Retrieves the coordinates and associated weight of the cubature point with the
// given 'index' from arrays of precomputed values. The implementation of this
// function is hard-coded based on those precomputed values.  It is also
// specifically for dimension 3.

// Four elements will be returned in 'coords': the first three elements are the
// coordinates of the cubature point and the fourth is the weight associated
// with that point.  The weight is the product of the weights associated with
// each coordinate value after application of the periodization transformation.

#endif  // _PRECOMPUTED_COORDS_H_

