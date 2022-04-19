//-*- C++ -*-
#include <stdlib.h>

#include "iosetting.h"
#include "mat2.h"

void Vec2::print(const char *msg) const {
  if (*msg) cout << msg << ":" << endl;
  cout << WID6 << PREC3 << x << WID6 << PREC3 << y << endl;
}

ostream& operator<<(ostream& s, const Vec2& v) {
  return s << '(' << v.x << ',' << v.y << ')';
}

istream& operator>>(istream& s, Vec2& v) {
/*
   Input Format
   FLOAT
   ( FLOAT )
   ( FLOAT, FLOAT )
*/
  double x=0.0, y=0.0;
  char c=0;
  s >> c;
  if ( c=='(' ) {
    s >> x >> c;
    if ( c==',' ) {
      s >> y >> c;
      if ( c!=')' ) s.clear(ios::badbit);
    }
    if ( c!=')' ) s.clear(ios::badbit);
  } else {
    s.putback(c);
    s >> x;
  }
  if (s) v = Vec2(x,y);
  return s;
}

void Mat2::print(const char *msg) const {
  if (*msg) cout << msg << ":" << endl;
  cout << WID6 << PREC3 << xx << WID6 << PREC3 << xy << endl;
  cout << WID6 << PREC3 << yx << WID6 << PREC3 << yy << endl;
}
