//GENIAL - GENeric Image & Array Library
//Copyright (C) 2005  IENT - RWTH Aachen
//
//This program is free software; you can redistribute it and/or
//modify it under the terms of the GNU General Public License
//as published by the Free Software Foundation; either version 2
//of the License, or (at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#ifndef SSE3_H
#define SSE3_H

#ifndef NO_SIMD

#ifndef SSE3
#define SSE3
#endif

#include <pmmintrin.h>
#include "sse2.h"

inline void loadu(m128b &x, const unsigned char *p) { x.vec = _mm_lddqu_si128((__m128i *)p); }

extern __m128 _mm_hadd_ps(__m128 a, __m128 b);
extern __m128 _mm_hsub_ps(__m128 a, __m128 b);

inline m128f hadd (const m128f &x, const m128f &y) { return m128f(_mm_hadd_ps(x.vec,y.vec)); }
inline m128f hsub (const m128f &x, const m128f &y) { return m128f(_mm_hsub_ps(x.vec,y.vec)); }
inline m128d hadd (const m128d &x, const m128d &y) { return m128d(_mm_hadd_pd(x.vec,y.vec)); }
inline m128d hsub (const m128d &x, const m128d &y) { return m128d(_mm_hsub_pd(x.vec,y.vec)); }

inline m128f sum  (const m128f &x0) { m128f a=hadd(x0,x0); return hadd(a,a);  }
inline m128f sumlo(const m128f &x0) { return hadd(x0,x0); }
inline m128f sumhi(const m128f &x0) { m128f a=movehl(x0,x0); return hadd(a,a); }
inline m128f sum  (const m128f &x0, const m128f &x1) { m128f a=hadd(x0,x1); return hadd(a,a); }
inline m128f sum  (const m128f &x0, const m128f &x1, const m128f &x2, const m128f &x3) { return hadd(hadd(x0,x1),hadd(x2,x3)); }

inline m128d sum(const m128d &x0) { return hadd(x0,x0); }
inline m128d sum(const m128d &x0, const m128d &x1) { return hadd(x0,x1); }

inline m128cf moveldup(const m128cf &x) { return m128cf(_mm_moveldup_ps(x.vec)); }
inline m128cf movehdup(const m128cf &x) { return m128cf(_mm_movehdup_ps(x.vec)); }

inline m128cd moveldup(const m128cd &x) { return m128cd(_mm_movedup_pd(x.vec)); }
inline m128cd movehdup(const m128cd &x) { return m128cd(_mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(1,1))); }


inline m64f  norm(const m128cf &x) { __m128  a=_mm_mul_ps(x.vec,x.vec); return m64f(_mm_hadd_ps(a,a)); }
inline m64d  norm(const m128cd &x) { __m128d a=_mm_mul_pd(x.vec,x.vec); return m64d(_mm_hadd_pd(a,a)); }
inline m128f norm(const m128cf &x,const m128cf &y) { __m128  a0=_mm_mul_ps(x.vec,x.vec); __m128  a1=_mm_mul_ps(y.vec,y.vec); return m128f(_mm_hadd_ps(a0,a1)); }
inline m128d norm(const m128cd &x,const m128cd &y) { __m128d a0=_mm_mul_pd(x.vec,x.vec); __m128d a1=_mm_mul_pd(y.vec,y.vec); return m128d(_mm_hadd_pd(a0,a1)); }

inline m128cf operator*(const m128cf &x, const m128cf &y) 
{
	__m128 a=_mm_mul_ps(moveldup(x).vec,y.vec);
	__m128 b=_mm_mul_ps(movehdup(x).vec, flip_ri(y).vec );
	return m128cf( _mm_addsub_ps(a,b) );
}
inline m128cd operator*(const m128cd &x, const m128cd &y)
{
	__m128d a=_mm_mul_pd(moveldup(x).vec,y.vec);
	__m128d b=_mm_mul_pd(movehdup(x).vec,flip_ri(y).vec);
	return m128cd( _mm_addsub_pd(a,b) );
}

inline m128cf cmul(const m128cf &x, const m128cf &y) 
{
	__m128 a=_mm_mul_ps(moveldup(x).vec,y.vec);
	__m128 b=_mm_mul_ps(movehdup(x).vec, flip_ri(y).vec );
	return m128cf( _mm_addsub_ps(a,_mm_sub_ps(_mm_setzero_ps(),b)) );
}
inline m128cd cmul(const m128cd &x, const m128cd &y)
{
	__m128d a=_mm_mul_pd(moveldup(x).vec,y.vec);
	__m128d b=_mm_mul_pd(movehdup(x).vec,flip_ri(y).vec);
	return m128cd( _mm_addsub_pd(a,_mm_sub_pd(_mm_setzero_pd(),b)) );
}

inline m128cf icmul(const m128cf &x, const m128cf &y) 
{
	__m128 a=_mm_mul_ps(movehdup(x).vec,y.vec);
	__m128 b=_mm_mul_ps(moveldup(x).vec, flip_ri(y).vec );
	return m128cf( _mm_addsub_ps(a,b) );
}
inline m128cd icmul(const m128cd &x, const m128cd &y)
{
	__m128d a=_mm_mul_pd(movehdup(x).vec,y.vec);
	__m128d b=_mm_mul_pd(moveldup(x).vec,flip_ri(y).vec);
	return m128cd( _mm_addsub_pd(a,b) );
}

inline m128cf mimul(const m128cf &x, const m128cf &y) 
{
	__m128 a=_mm_mul_ps(movehdup(x).vec,y.vec);
	__m128 b=_mm_mul_ps(moveldup(x).vec, flip_ri(y).vec );
	return m128cf( _mm_addsub_ps(a,_mm_sub_ps(_mm_setzero_ps(),b)) );
}
inline m128cd mimul(const m128cd &x, const m128cd &y)
{
	__m128d a=_mm_mul_pd(movehdup(x).vec,y.vec);
	__m128d b=_mm_mul_pd(moveldup(x).vec,flip_ri(y).vec);
	return m128cd( _mm_addsub_pd(a,_mm_sub_pd(_mm_setzero_pd(),b)) );
}

inline m128cf imul (const m128cf &x, const m128cf &y) 
{
	__m128 a=_mm_mul_ps(movehdup(x).vec,y.vec);
	__m128 b=_mm_mul_ps(moveldup(x).vec, flip_ri(y).vec );
	return m128cf( _mm_addsub_ps(_mm_sub_ps(_mm_setzero_ps(),a),b) );
}
inline m128cd imul (const m128cd &x, const m128cd &y)
{
	__m128d a=_mm_mul_pd(movehdup(x).vec,y.vec);
	__m128d b=_mm_mul_pd(moveldup(x).vec,flip_ri(y).vec);
	return m128cd( _mm_addsub_pd(_mm_sub_pd(_mm_setzero_pd(),a),b) );
}

inline m128cf subadd(const m128cf &x, const m128cf &y) { return m128cf(_mm_addsub_ps(x.vec,y.vec)); }
inline m128cd subadd(const m128cd &x, const m128cd &y) { return m128cd(_mm_addsub_pd(x.vec,y.vec)); }

inline m128cf addsub(const m128cf &x, const m128cf &y) { return m128cf(subadd(x,-y)); }
inline m128cd addsub(const m128cd &x, const m128cd &y) { return m128cd(subadd(x,-y)); }

inline m128cf mimul(const m128cf &x) { return flip_ri(x); }
inline m128cd mimul(const m128cd &x) { return flip_ri(x); }

inline m128cf imul(const m128cf &x) { return subadd(m128cf(_mm_setzero_ps()), flip_ri(x) ); }
inline m128cd imul(const m128cd &x) { return subadd(m128cd(_mm_setzero_pd()), flip_ri(x) ); }

inline m128cf add_imul (const m128cf &a, const m128cf &x) { return subadd(a,flip_ri(x)); }
inline m128cd add_imul (const m128cd &a, const m128cd &x) { return subadd(a,flip_ri(x)); }

inline m128cf sub_imul (const m128cf &a, const m128cf &x) { return a-imul(x); }
inline m128cd sub_imul (const m128cd &a, const m128cd &x) { return a-imul(x); }


//inline m128cf operator/(const m128cf &x, const m128cf & y)
//{
//	__m128 a=_mm_mul_ps( _mm_movehdup_ps(x.vec), y.vec);
//	__m128 b=_mm_mul_ps( _mm_moveldup_ps(x.vec), flip_ri(y));
//	__m128 c=_mm_add_ps( _mm_mul_ps(y.vec,y.vec), _mm_mul_ps(flip_ri(y),flip_ri(y)) );
//	return m128cf( _mm_div_ps(flip_ri(_mm_addsub_ps(a,b)),c) );
//}
//
//inline m128cd operator/(const m128cd &x, const m128cd &y)
//{
//	__m128d a=_mm_mul_pd( _mm_shuffle_pd(x.vec,x.vec,_MM_SHUFFLE2(1,1)),y.vec);
//	__m128d b=_mm_mul_pd( flip_ri(y).vec, _mm_movedup_pd(x.vec));
//	__m128d c=flip_ri(_mm_addsub_pd(a,b));
//  __m128d d=_mm_add_pd(_mm_mul_pd(flip_ri(y).vec,flip_ri(y).vec),_mm_mul_pd(y.vec,y.vec));
//	return m128cd( _mm_div_pd(c,d) );
//}



#endif // NO_SIMD

#endif

