#include <iostream>
#include <iomanip>
using namespace std;

#include "BigInt.hh"
#include "MpIeee.hh"
#include "ArithmosIO.hh"
#include <vector>


BigInt matrix_print_off(BigInt nr, BigInt nc, MpIeee **A);
BigInt vector_print_off(BigInt nr, MpIeee *x);
void gauss(MpIeee **a, MpIeee *b, MpIeee *x, BigInt n);

BigInt matrix_print_off(BigInt nr, BigInt nc, MpIeee **A);
BigInt vector_print_off(BigInt nr, MpIeee x);
void gauss(MpIeee **a, MpIeee *b, MpIeee *x, BigInt n);

