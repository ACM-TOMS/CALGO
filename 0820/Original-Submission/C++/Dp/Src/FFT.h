//***************************************************************
// FFT.h
//****************************************************************

#ifndef FFT_LIBRARY_H
#define FFT_LIBRARY_H

integer rev(integer In, integer J);
void    bitReverse(Interval &In);
void Weights(Interval& Rew, Interval& Imw, integer flag);
void CoreFFT(Interval& ReA, Interval& ImA, Interval& C, Interval& S);

void FFT(Interval &ReIn, Interval &ImIn);

#endif
