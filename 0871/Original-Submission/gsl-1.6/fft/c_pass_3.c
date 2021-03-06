/* fft/c_pass_3.c
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

static int
FUNCTION(fft_complex,pass_3) (const BASE in[],
                              const size_t istride,
                              BASE out[],
                              const size_t ostride,
                              const gsl_fft_direction sign,
                              const size_t product,
                              const size_t n,
                              const TYPE(gsl_complex) * twiddle1,
                              const TYPE(gsl_complex) * twiddle2)
{
  size_t i = 0, j = 0;
  size_t k, k1;

  const size_t factor = 3;
  const size_t m = n / factor;
  const size_t q = n / product;
  const size_t product_1 = product / factor;
  const size_t jump = (factor - 1) * product_1;

  const ATOMIC tau = sqrt (3.0) / 2.0;

  for (k = 0; k < q; k++)
    {
      ATOMIC w1_real, w1_imag, w2_real, w2_imag;

      if (k == 0)
        {
          w1_real = 1.0;
          w1_imag = 0.0;
          w2_real = 1.0;
          w2_imag = 0.0;
        }
      else
        {
          if (sign == forward)
            {
              /* forward tranform */
              w1_real = GSL_REAL(twiddle1[k - 1]);
              w1_imag = GSL_IMAG(twiddle1[k - 1]);
              w2_real = GSL_REAL(twiddle2[k - 1]);
              w2_imag = GSL_IMAG(twiddle2[k - 1]);
            }
          else
            {
              /* backward tranform: w -> conjugate(w) */
              w1_real = GSL_REAL(twiddle1[k - 1]);
              w1_imag = -GSL_IMAG(twiddle1[k - 1]);
              w2_real = GSL_REAL(twiddle2[k - 1]);
              w2_imag = -GSL_IMAG(twiddle2[k - 1]);
            }
        }

      for (k1 = 0; k1 < product_1; k1++)
        {
          const ATOMIC z0_real = REAL(in,istride,i);
          const ATOMIC z0_imag = IMAG(in,istride,i);
          const ATOMIC z1_real = REAL(in,istride,i+m);
          const ATOMIC z1_imag = IMAG(in,istride,i+m);
          const ATOMIC z2_real = REAL(in,istride,i+2*m);
          const ATOMIC z2_imag = IMAG(in,istride,i+2*m);

          /* compute x = W(3) z */

          /* t1 = z1 + z2 */
          const ATOMIC t1_real = z1_real + z2_real;
          const ATOMIC t1_imag = z1_imag + z2_imag;
          
          /* t2 = z0 - t1/2 */
          const ATOMIC t2_real = z0_real - t1_real / 2.0;
          const ATOMIC t2_imag = z0_imag - t1_imag / 2.0;
          
          /* t3 = (+/-) sin(pi/3)*(z1 - z2) */
          const ATOMIC t3_real = ((int) sign) * tau * (z1_real - z2_real);
          const ATOMIC t3_imag = ((int) sign) * tau * (z1_imag - z2_imag);
          
          /* x0 = z0 + t1 */
          const ATOMIC x0_real = z0_real + t1_real;
          const ATOMIC x0_imag = z0_imag + t1_imag;
          
          /* x1 = t2 + i t3 */
          const ATOMIC x1_real = t2_real - t3_imag;
          const ATOMIC x1_imag = t2_imag + t3_real;
          
          /* x2 = t2 - i t3 */
          const ATOMIC x2_real = t2_real + t3_imag;
          const ATOMIC x2_imag = t2_imag - t3_real;

          /* apply twiddle factors */

          /* to0 = 1 * x0 */
          REAL(out,ostride,j) = x0_real;
          IMAG(out,ostride,j) = x0_imag;
          
          /* to1 = w1 * x1 */
          REAL(out,ostride,j+product_1) = w1_real * x1_real - w1_imag * x1_imag;
          IMAG(out,ostride,j+product_1) = w1_real * x1_imag + w1_imag * x1_real;
          
          /* to2 = w2 * x2 */
          REAL(out,ostride,j+2*product_1) = w2_real * x2_real - w2_imag * x2_imag;
          IMAG(out,ostride,j+2*product_1) = w2_real * x2_imag + w2_imag * x2_real;

          i++; j++;
        }
      j += jump;
    }
  return 0;
}
