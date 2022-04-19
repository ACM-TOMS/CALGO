//-*- C++ -*-
#include <stdlib.h>

#include "iosetting.h"
#include "mat3.h"

void Vec3::print(const char *msg) const {
  if (*msg) cout << msg << ":" << endl;
  cout << WID6 << PREC3 << x << WID6 << PREC3 << y << WID6 << PREC3 << z << endl;
}

ostream& operator<<(ostream& s, const Vec3& v) {
  return s << '(' << v.x << ',' << v.y << ',' << v.z << ')';
}

istream& operator>>(istream& s, Vec3& v) {
/*
   Input Format
   FLOAT
   ( FLOAT )
   ( FLOAT, FLOAT )
   ( FLOAT, FLOAT, FLOAT)
*/
  double x=0.0, y=0.0, z=0.0;
  char c=0;
  s >> c;
  if ( c=='(' ) {
    s >> x >> c;
    if ( c==',' ) {
      s >> y >> c;
      if ( c==',' ) s >> z >> c;
      if ( c!=')' ) s.clear(ios::badbit);
    }
    if ( c!=')' ) s.clear(ios::badbit);
  } else {
    s.putback(c);
    s >> x;
  }
  if (s) v = Vec3(x,y,z);
  return s;
}

void Mat3::print(const char *msg) const {
  if (*msg) cout << msg << ":" << endl;
  cout << WID6 << PREC3 << xx << WID6 << PREC3 << xy << WID6 << PREC3 << xz << endl;
  cout << WID6 << PREC3 << yx << WID6 << PREC3 << yy << WID6 << PREC3 << yz << endl;
  cout << WID6 << PREC3 << zx << WID6 << PREC3 << zy << WID6 << PREC3 << zz << endl;
}

void DiagMat3::print(const char *msg) const {
  if (*msg) cout << msg << ":" << endl;
  cout << WID6 << PREC3 << xx << endl;
  cout << WID6 << PREC3 << yy << endl;
  cout << WID6 << PREC3 << zz << endl;
}

