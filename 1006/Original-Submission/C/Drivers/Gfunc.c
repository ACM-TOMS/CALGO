/**********************************************************************************

  DELTAGAMMAINC Fast and Accurate Evaluation of a Generalized Incomplete Gamma
  Function. Copyright (C) 2016 Remy Abergel (remy.abergel AT gmail.com), Lionel
  Moisan (Lionel.Moisan AT parisdescartes.fr).

  This file is a part of the DELTAGAMMAINC software, dedicated to the
  computation of a generalized incomplete gammafunction. See the Companion paper
  for a complete description of the algorithm.

  ``Fast and accurate evaluation of a generalized incomplete gamma function''
  (Rémy Abergel, Lionel Moisan), preprint MAP5 nº2016-14, revision 1.

  This program is free software: you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************/

/*
 * Compilation (with gcc): gcc -O3 kernel.c Gfunc.c -lm -o Gfunc
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern void G_func(double*,double,double); // see source file 'kernel.c'

static void display_usage()
{
  printf("Usage: Gfunc [--help] x p\n\n");
  printf("   --help : display help\n");
  printf("   x      : (double) a real number, possibly infinite (-inf < x <= inf)\n");
  printf("   p      : (double) a positive real number (p > 0)\n\n");
  printf("Description: compute G(p,x) such as\n\n");
  printf("  if x <= p: G(p,x) = exp(x-p*log(|x|)) * integral over [0,|x|] of s^{p-1} * exp(-sign(x)*s) ds\n");
  printf("  otherwise: G(p,x) =  exp(x-p*log(x))  * integral over [x,inf] of s^{p-1} * exp(-s) ds\n\n");
}

int main(int argc, char **argv)
{
  char *x_value = NULL,*p_value = NULL;
  double G,x,p;
  int err,i;

  /* parser */
  for(i=1;i<argc;i++) {
    if (strcmp(argv[i],"--help") == 0) {
      display_usage();
      return EXIT_SUCCESS;
    }
    else if(NULL == x_value) {
      x_value = argv[i];
      err = sscanf(x_value,"%lf",&x);
      if(err != 1) {
	printf("\nError: could not retrieve properly input argument 'x'.\n");
	display_usage();
	return EXIT_FAILURE;
      }
    }
    else if(NULL == p_value) {
      p_value = argv[i];
      err = sscanf(p_value,"%lf",&p);
      if(err != 1) {
	printf("\nError: could not retrieve properly input argument 'p'.\n");
	display_usage();
	return EXIT_FAILURE;
      }
    }
  }
  if ((NULL == x_value) || (NULL == p_value)) {
    printf("\nError: invalid number of arguments.\n");
    display_usage();
    return EXIT_FAILURE;
  }

  /* consistency checks */
  if(p <= 0) { printf("\nError: input 'p' must be positive (p > 0)!\nType Gfunc --help to display help.\n\n"); return EXIT_FAILURE; }

  /*** evaluate G(p,x) ***/
  G_func(&G,p,x);

  /*** printf results ***/
  if (p >= x) printf("\nComputing G(p,x) = exp(x-p*log(|x|)) * integral over [0,|x|] of s^{p-1} * exp(-sign(x)*s) ds,\n");
  else printf("\nComputing G(p,x) = exp(x-p*log(x)) * integral over [x,inf] of s^{p-1} * exp(-s) ds,\n");

  printf("where (hexadecimal): x=%a, p=%a,\n",x,p);
  printf("(in decimal approx): x=%.34g, p=%.34g.\n\n",x,p);

  printf("Computed G(p,x) = %.17e\n\n",G);

  return EXIT_SUCCESS;
}
