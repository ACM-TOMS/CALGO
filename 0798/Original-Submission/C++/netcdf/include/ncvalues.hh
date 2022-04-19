/*********************************************************************
 *   Copyright 1992, University Corporation for Atmospheric Research
 *   See netcdf/README file for copying and redistribution conditions.
 *
 *   Purpose:	interface for classes of typed arrays for netCDF
 *
 *   $Header: /upc/share/CVS/netcdf-3/cxx/ncvalues.hh,v 1.1 1997/04/21 15:44:53 russ Exp $
 *********************************************************************/

#ifndef Ncvalues_def
#define Ncvalues_def

#include <generic.h>
#include <iostream.h>
#ifdef STRSTREAM_H_SPEC
#   include STRSTREAM_H_SPEC
#else
#   include <strstream.h>
#endif
#include <limits.h>
#include <string.h>
#include "netcdf.h"

typedef unsigned char ncbyte;

#define NC_UNSPECIFIED ((nc_type)0)

enum NcType 
{
  ncNoType = NC_UNSPECIFIED, 
  ncByte = NC_BYTE, 
  ncChar = NC_CHAR, 
  ncShort = NC_SHORT, 
  ncLong = NC_LONG, 
  ncFloat = NC_FLOAT, 
  ncDouble = NC_DOUBLE
};

#define ncBad_ncbyte ncBad_byte
static const ncbyte ncBad_byte = FILL_BYTE;
static const char ncBad_char = FILL_CHAR;
static const short ncBad_short = FILL_SHORT;
static const nclong ncBad_nclong = FILL_LONG;
static const long ncBad_long = FILL_LONG;
static const float ncBad_float = FILL_FLOAT;
static const double ncBad_double = FILL_DOUBLE;

// This is the same as the name2 macro from generic.h, but we need to define
// our own version since rescanning something generated with the name2 macro
// won't necessarily cause name2 to be expanded again.
#define makename2(z, y)		makename2_x(z, y)
#define makename2_x(z, y)		z##y

#define NcVal(TYPE) makename2(NcValues_,TYPE)

#define NcValuesdeclare(TYPE)						      \
class NcVal(TYPE) : public NcValues					      \
{									      \
  public:								      \
    NcVal(TYPE)( void );						      \
    NcVal(TYPE)(long num);						      \
    NcVal(TYPE)(long num, const TYPE* vals);				      \
    NcVal(TYPE)(const NcVal(TYPE)&);					      \
    virtual NcVal(TYPE)& operator=(const NcVal(TYPE)&);			      \
    virtual ~NcVal(TYPE)( void );					      \
    virtual void* base( void ) const;					      \
    virtual int bytes_for_one( void ) const;				      \
    virtual ncbyte as_ncbyte( long n ) const;				      \
    virtual char as_char( long n ) const;				      \
    virtual short as_short( long n ) const;				      \
    virtual nclong as_nclong( long n ) const;				      \
    virtual long as_long( long n ) const;				      \
    virtual float as_float( long n ) const;				      \
    virtual double as_double( long n ) const;				      \
    virtual char* as_string( long n ) const;				      \
    virtual int invalid( void ) const;					      \
  private:								      \
    TYPE* the_values;							      \
    ostream& print(ostream&) const;					      \
};

#define NcTypeEnum(TYPE) makename2(_nc__,TYPE)
#define _nc__ncbyte ncByte
#define _nc__char ncChar
#define _nc__short ncShort
#define _nc__nclong ncLong
#define _nc__float ncFloat
#define _nc__double ncDouble
#define NcValuesimplement(TYPE)						      \
NcVal(TYPE)::NcVal(TYPE)( void )					      \
	: NcValues(NcTypeEnum(TYPE), 0), the_values(0)			      \
{}									      \
									      \
NcVal(TYPE)::NcVal(TYPE)(long num, const TYPE* vals)			      \
	: NcValues(NcTypeEnum(TYPE), num)				      \
{									      \
    the_values = new TYPE[num];						      \
    for(int i = 0; i < num; i++)					      \
      the_values[i] = vals[i];						      \
}									      \
									      \
NcVal(TYPE)::NcVal(TYPE)(long num)					      \
	: NcValues(NcTypeEnum(TYPE), num), the_values(new TYPE[num])	      \
{}									      \
									      \
NcVal(TYPE)::NcVal(TYPE)(const NcVal(TYPE)& v)				      \
{									      \
    delete[] the_values;						      \
    the_values = new TYPE[v.the_number];				      \
    for(int i = 0; i < v.the_number; i++)				      \
      the_values[i] = v.the_values[i];					      \
}									      \
									      \
NcVal(TYPE)& NcVal(TYPE)::operator=(const NcVal(TYPE)& v)		      \
{									      \
    delete[] the_values;						      \
    the_values = new TYPE[v.the_number];				      \
    for(int i = 0; i < v.the_number; i++)				      \
      the_values[i] = v.the_values[i];					      \
    return *this;							      \
}									      \
									      \
void* NcVal(TYPE)::base( void ) const					      \
{									      \
    return the_values;							      \
}									      \
									      \
NcVal(TYPE)::~NcVal(TYPE)( void )					      \
{									      \
    delete[] the_values;						      \
}									      \
									      \
int NcVal(TYPE)::invalid( void ) const					      \
{									      \
    for(int i=0;i<the_number;i++)					      \
	if (the_values[i] == makename2(ncBad_,TYPE)) return 1; 	       	      \
    return 0;                                                                 \
}                                                                             \


#define Ncbytes_for_one_implement(TYPE)					      \
int NcVal(TYPE)::bytes_for_one( void ) const				      \
{									      \
    return nctypelen((nc_type) NcTypeEnum(TYPE));			      \
}

#define as_ncbyte_implement(TYPE)					      \
ncbyte NcVal(TYPE)::as_ncbyte( long n ) const				      \
{									      \
    if (the_values[n] < 0 || the_values[n] > UCHAR_MAX)		              \
      return ncBad_byte;						      \
    return (ncbyte) the_values[n];			                      \
}

#define as_char_implement(TYPE)						      \
char NcVal(TYPE)::as_char( long n ) const				      \
{									      \
    if (the_values[n] < CHAR_MIN || the_values[n] > CHAR_MAX)		      \
      return ncBad_char;						      \
    return (char) the_values[n];					      \
}

#define as_short_implement(TYPE)					      \
short NcVal(TYPE)::as_short( long n ) const				      \
{									      \
    if (the_values[n] < SHRT_MIN || the_values[n] > SHRT_MAX)		      \
      return ncBad_short;						      \
    return (short) the_values[n];				              \
}

#define NCLONG_MIN INT_MIN
#define NCLONG_MAX INT_MAX
#define as_nclong_implement(TYPE)					      \
nclong NcVal(TYPE)::as_nclong( long n ) const				      \
{									      \
    if (the_values[n] < NCLONG_MIN || the_values[n] > NCLONG_MAX)	      \
      return ncBad_nclong;						      \
    return (nclong) the_values[n];				              \
}

#define as_long_implement(TYPE)					              \
long NcVal(TYPE)::as_long( long n ) const				      \
{									      \
    if (the_values[n] < LONG_MIN || the_values[n] > LONG_MAX)	              \
      return ncBad_long;						      \
    return (long) the_values[n];		       	                      \
}

#define as_float_implement(TYPE)					      \
inline float NcVal(TYPE)::as_float( long n ) const			      \
{									      \
    return (float) the_values[n];				              \
}

#define as_double_implement(TYPE)					      \
inline double NcVal(TYPE)::as_double( long n ) const			      \
{									      \
    return (double) the_values[n];				              \
}

#define as_string_implement(TYPE)					      \
char* NcVal(TYPE)::as_string( long n ) const				      \
{									      \
    char* s = new char[32];						      \
    ostrstream(s, sizeof(s)) << the_values[n] << ends ;			      \
    return s;								      \
}

class NcValues			// ABC for value blocks
{
  public:
    NcValues( void );
    NcValues(NcType, long);
    virtual ~NcValues( void );
    virtual long num( void );
    virtual ostream& print(ostream&) const = 0;
    virtual void* base( void ) const = 0;
    virtual int bytes_for_one( void ) const = 0;

    // The following member functions provide conversions from the value
    // type to a desired basic type.  If the value is out of range, the
    // default "fill-value" for the appropriate type is returned.
    virtual ncbyte as_ncbyte( long n ) const = 0; // nth value as a byte
    virtual char as_char( long n ) const = 0;     // nth value as char
    virtual short as_short( long n ) const = 0;   // nth value as short
    virtual nclong as_nclong( long n ) const = 0; // nth value as nclong
    virtual long as_long( long n ) const = 0;     // nth value as long
    virtual float as_float( long n ) const = 0;   // nth value as floating-point
    virtual double as_double( long n ) const = 0; // nth value as double
    virtual char* as_string( long n ) const = 0;  // value as string
    
  protected:
    NcType the_type;
    long the_number;
    friend ostream& operator<< (ostream&, const NcValues&);
};

declare(NcValues,ncbyte)
declare(NcValues,char)
declare(NcValues,short)
declare(NcValues,nclong)
declare(NcValues,float)
declare(NcValues,double)

#endif
