#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* fft/signals.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

int  FUNCTION(fft_signal,complex_pulse) (const size_t k, 
                                        const size_t n,
                                        const size_t stride,
                                        const BASE z_real, 
                                        const BASE z_imag,
                                        MpIeee data[],
                                        MpIeee fft[]);

int  FUNCTION(fft_signal,complex_constant) (const size_t n,
                                           const size_t stride,
                                           const BASE z_real,
                                           const BASE z_imag,
                                           MpIeee data[],
                                           MpIeee fft[]);

int  FUNCTION(fft_signal,complex_exp) (const int k,
                                      const size_t n,
                                      const size_t stride,
                                      const BASE z_real,
                                      const BASE z_imag,
                                      MpIeee data[],
                                      MpIeee fft[]);


int  FUNCTION(fft_signal,complex_exppair) (const int k1,
                                          const int k2,
                                          const size_t n,
                                          const size_t stride,
                                          const BASE z1_real,
                                          const BASE z1_imag,
                                          const BASE z2_real,
                                          const BASE z2_imag,
                                          MpIeee data[],
                                          MpIeee fft[]);

int  FUNCTION(fft_signal,complex_noise) (const size_t n,
                                        const size_t stride,
                                        MpIeee data[],
                                        MpIeee fft[]);

int  FUNCTION(fft_signal,real_noise) (const size_t n,
                                     const size_t stride,
                                     MpIeee data[],
                                     MpIeee fft[]);

