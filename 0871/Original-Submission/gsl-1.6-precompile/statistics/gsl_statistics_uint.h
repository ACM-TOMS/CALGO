#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/gsl_statistics_uint.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Jim Davies, Brian Gough
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

#ifndef __GSL_STATISTICS_UINT_H__
#define __GSL_STATISTICS_UINT_H__

#include <stddef.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

MpIeee gsl_stats_uint_mean (const unsigned int  data[], const size_t stride, const size_t n);
MpIeee gsl_stats_uint_variance(const unsigned int  data[], const size_t stride, const size_t n);
MpIeee gsl_stats_uint_sd(const unsigned int  data[], const size_t stride, const size_t n);
MpIeee gsl_stats_uint_variance_with_fixed_mean(const unsigned int  data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_uint_sd_with_fixed_mean(const unsigned int  data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_uint_absdev(const unsigned int  data[], const size_t stride, const size_t n);
MpIeee gsl_stats_uint_skew(const unsigned int  data[], const size_t stride, const size_t n);
MpIeee gsl_stats_uint_kurtosis(const unsigned int  data[], const size_t stride, const size_t n);
MpIeee gsl_stats_uint_lag1_autocorrelation(const unsigned int  data[], const size_t stride, const size_t n);

MpIeee gsl_stats_uint_covariance(const unsigned int  data1[], const size_t stride1,const unsigned int  data2[], const size_t stride2, const size_t n);

MpIeee gsl_stats_uint_variance_m(const unsigned int  data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_uint_sd_m(const unsigned int  data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_uint_absdev_m(const unsigned int  data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_uint_skew_m_sd(const unsigned int  data[], const size_t stride, const size_t n, const MpIeee mean, const MpIeee sd);
MpIeee gsl_stats_uint_kurtosis_m_sd(const unsigned int  data[], const size_t stride, const size_t n, const MpIeee mean, const MpIeee sd);
MpIeee gsl_stats_uint_lag1_autocorrelation_m(const unsigned int  data[], const size_t stride, const size_t n, const MpIeee mean);

MpIeee gsl_stats_uint_covariance_m(const unsigned int  data1[], const size_t stride1,const unsigned int  data2[], const size_t stride2, const size_t n, const MpIeee mean1, const MpIeee mean2);


MpIeee gsl_stats_uint_pvariance(const unsigned int  data1[], const size_t stride1, const size_t n1, const unsigned int  data2[], const size_t stride2, const size_t n2);
MpIeee gsl_stats_uint_ttest(const unsigned int  data1[], const size_t stride1, const size_t n1, const unsigned int  data2[], const size_t stride2, const size_t n2);

unsigned int  gsl_stats_uint_max(const unsigned int  data[], const size_t stride, const size_t n);
unsigned int  gsl_stats_uint_min(const unsigned int  data[], const size_t stride, const size_t n);
void gsl_stats_uint_minmax (unsigned int  * min, unsigned int  * max, const unsigned int  data[], const size_t stride, const size_t n);

size_t gsl_stats_uint_max_index (const unsigned int  data[], const size_t stride, const size_t n);
size_t gsl_stats_uint_min_index (const unsigned int  data[], const size_t stride, const size_t n);
void gsl_stats_uint_minmax_index (size_t * min_index, size_t * max_index, const unsigned int  data[], const size_t stride, const size_t n);

MpIeee gsl_stats_uint_median_from_sorted_data(const unsigned int  sorted_data[], const size_t stride, const size_t n) ;
MpIeee gsl_stats_uint_quantile_from_sorted_data(const unsigned int  sorted_data[], const size_t stride, const size_t n, const MpIeee f) ;

__END_DECLS

#endif /* __GSL_STATISTICS_UINT_H__ */
