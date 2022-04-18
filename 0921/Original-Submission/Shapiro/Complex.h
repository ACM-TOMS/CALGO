/*
    Jonathan Hauenstein and Frank Sottile
    Copyright 2010
*/
#include <iostream>
#include <gmp.h>

class complex_rational {
public:
  mpq_t re, im;

  // constructor (initialize to 0)
  complex_rational() { mpq_init(re); mpq_set_ui(re, 0, 1); mpq_init(im); mpq_set_ui(im, 0, 1); }
  // constructor (initialize to a)
  complex_rational(const complex_rational& a) { mpq_init(re); mpq_set(re, a.re); mpq_init(im); mpq_set(im, a.im); }
  // constructor (initialize to a)
  complex_rational(const mpq_t& a) { mpq_init(re); mpq_set(re, a); mpq_init(im); mpq_set_ui(im, 0, 1); }
  // constructor (initialize to a)
  complex_rational(const int& a) { mpq_init(re); mpq_set_si(re, a, 1); mpq_init(im); mpq_set_ui(im, 0, 1); }

  // initialize & clear
  void init() { mpq_init(this->re); mpq_set_ui(this->re, 0, 1); mpq_init(this->im); mpq_set_ui(this->im, 0, 1); }
  void init(const complex_rational& a) { mpq_init(this->re); mpq_set(this->re, a.re); mpq_init(this->im); mpq_set(this->im, a.im); }
  void init(const mpq_t& a) { mpq_init(this->re); mpq_set(this->re, a); mpq_init(this->im); mpq_set_ui(this->im, 0, 1); }
  void init(const int& a) { mpq_init(this->re); mpq_set_ui(this->re, a, 1); mpq_init(this->im); mpq_set_ui(this->im, 0, 1); }
  void clear() { mpq_clear(this->re); mpq_clear(this->im); }

  // set zero, one, I
  complex_rational& set_zero() { mpq_set_ui(this->re, 0, 1); mpq_set_ui(this->im, 0, 1); return *this; }
  complex_rational& set_one() { mpq_set_ui(this->re, 1, 1); mpq_set_ui(this->im, 0, 1); return *this; }
  complex_rational& set_I() { mpq_set_ui(this->re, 0, 1); mpq_set_ui(this->im, 1, 1); return *this; }
  complex_rational& set_mpq(mpq_t re, mpq_t im) { mpq_set(this->re, re); mpq_set(this->im, im); return *this; }

  // conjugate
  complex_rational& conjugate() { mpq_neg(this->im, this->im); return *this; }

  // print string
  void print_str(FILE*& ptr) { mpq_out_str(ptr,10,this->re); fprintf(ptr, " "); mpq_out_str(ptr,10,this->im); fprintf(ptr, "\n"); }
  void print_str(FILE*& ptr, const int& base) { mpq_out_str(ptr,base,this->re); fprintf(ptr, " "); mpq_out_str(ptr,base,this->im); fprintf(ptr, "\n"); }

  // =
  complex_rational& operator=(const complex_rational &rhs) { mpq_set(this->re, rhs.re); mpq_set(this->im, rhs.im); return *this; }

  // +=
  complex_rational& operator+=(const complex_rational &rhs) { mpq_add(this->re, this->re, rhs.re); mpq_add(this->im, this->im, rhs.im); return *this; }

  // +
  const complex_rational operator+(const complex_rational &other) const { complex_rational result(*this); result += other; return result; }

  // -=
  complex_rational& operator-=(const complex_rational &rhs) { mpq_sub(this->re, this->re, rhs.re); mpq_sub(this->im, this->im, rhs.im); return *this; }

  // -
  const complex_rational operator-(const complex_rational &other) const { complex_rational result(*this); result -= other; return result; }

  // *=
  complex_rational& operator*=(const complex_rational &rhs) { complex_rational temp, result(*this); mpq_mul(temp.re, result.re, rhs.re); mpq_mul(temp.im, result.im, rhs.im); mpq_sub(temp.re, temp.re, temp.im); mpq_mul(temp.im, result.re, rhs.im); mpq_mul(this->im, result.im, rhs.re); mpq_add(this->im, this->im, temp.im); mpq_set(this->re, temp.re); return *this; }

  // *
  const complex_rational operator*(const complex_rational &other) const { complex_rational result(*this); result *= other; return result; }

  // = (mpq_t)
  complex_rational& operator=(const mpq_t rhs) { mpq_set(this->re, rhs); mpq_set_ui(this->im, 0, 1); return *this; }

  // += (mpq_t)
  complex_rational& operator+=(const mpq_t rhs) { mpq_add(this->re, this->re, rhs); return *this; }

  // + (mpq_t)
  const complex_rational operator+(const mpq_t other) const { complex_rational result(*this); result += other; return result; }

  // -= (mpq_t)
  complex_rational& operator-=(const mpq_t rhs) { mpq_sub(this->re, this->re, rhs); return *this; }

  // - (mpq_t)
  const complex_rational operator-(const mpq_t other) const { complex_rational result(*this); result -= other; return result; }

  // *= (mpq_t)
  complex_rational& operator*=(const mpq_t rhs) { mpq_mul(this->re, this->re, rhs); mpq_mul(this->im, this->im, rhs); return *this; }

  // * (mpq_t)
  const complex_rational operator*(const mpq_t other) const { complex_rational result(*this); result *= other; return result; }

  // = (int)
  complex_rational& operator=(const int rhs) { mpq_set_si(this->re, rhs, 1); mpq_set_ui(this->im, 0, 1); return *this; }

  // += (int)
  complex_rational& operator+=(const int rhs) { complex_rational result(rhs); mpq_add(this->re, this->re, result.re); return *this; }

  // + (int)
  const complex_rational operator+(const int other) const { complex_rational result(*this); result += other; return result; }

  // -= (int)
  complex_rational& operator-=(const int rhs) { complex_rational result(rhs); mpq_sub(this->re, this->re, result.re); return *this; }

  // - (int)
  const complex_rational operator-(const int other) const { complex_rational result(*this); result -= other; return result; }

  // *= (int)
  complex_rational& operator*=(const int rhs) { complex_rational result(rhs); mpq_mul(this->re, this->re, result.re); mpq_mul(this->im, this->im, result.re); return *this; }

  // * (int)
  const complex_rational operator*(const int other) const { complex_rational result(*this); result *= other; return result; }

  // destructor (clear)
  ~complex_rational() { mpq_clear(re); mpq_clear(im); }
};

// - (negate)
complex_rational operator-(const complex_rational a) { complex_rational result; mpq_neg(result.re, a.re); mpq_neg(result.im, a.im); return result; }

// other operations
complex_rational operator+(mpq_t rat, complex_rational cmp) { return cmp + rat; }
complex_rational operator-(mpq_t rat, complex_rational cmp) { return cmp - rat; }
complex_rational operator*(mpq_t rat, complex_rational cmp) { return cmp * rat; }
complex_rational conjugate(complex_rational cmp) { return cmp.conjugate(); }
complex_rational operator+(int j, complex_rational cmp) { return cmp + j; }
complex_rational operator-(int j, complex_rational cmp) { return cmp - j; }
complex_rational operator*(int j, complex_rational cmp) { return cmp * j; }

void print_str(FILE* ptr, const int digits, complex_rational cmp) { mpq_out_str(ptr, digits, cmp.re); fprintf(ptr, " "); mpq_out_str(ptr, digits,cmp.im); fprintf(ptr, "\n"); }




