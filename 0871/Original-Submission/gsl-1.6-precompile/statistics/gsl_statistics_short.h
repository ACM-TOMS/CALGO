#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/gsl_statistics_short.h
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

#ifndef __GSL_STATISTICS_SHORT_H__
#define __GSL_STATISTICS_SHORT_H__

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

MpIeee gsl_stats_short_mean (const short data[], const size_t stride, const size_t n);
MpIeee gsl_stats_short_variance(const short data[], const size_t stride, const size_t n);
MpIeee gsl_stats_short_sd(const short data[], const size_t stride, const size_t n);
MpIeee gsl_stats_short_variance_with_fixed_mean(const short data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_short_sd_with_fixed_mean(const short data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_short_absdev(const short data[], const size_t stride, const size_t n);
MpIeee gsl_stats_short_skew(const short data[], const size_t stride, const size_t n);
MpIeee gsl_stats_short_kurtosis(const short data[], const size_t stride, const size_t n);
MpIeee gsl_stats_short_lag1_autocorrelation(const short data[], const size_t stride, const size_t n);

MpIeee gsl_stats_short_covariance(const short data1[], const size_t stride1,const short data2[], const size_t stride2, const size_t n);

MpIeee gsl_stats_short_variance_m(const short data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_short_sd_m(const short data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_short_absdev_m(const short data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_short_skew_m_sd(const short data[], const size_t stride, const size_t n, const MpIeee mean, const MpIeee sd);
MpIeee gsl_stats_short_kurtosis_m_sd(const short data[], const size_t stride, const size_t n, const MpIeee mean, const MpIeee sd);
MpIeee gsl_stats_short_lag1_autocorrelation_m(const short data[], const size_t stride, const size_t n, const MpIeee mean);

MpIeee gsl_stats_short_covariance_m(const short data1[], const size_t stride1,const short data2[], const size_t stride2, const size_t n, const MpIeee mean1, const MpIeee mean2);


MpIeee gsl_stats_short_pvariance(const short data1[], const size_t stride1, const size_t n1, const short data2[], const size_t stride2, const size_t n2);
MpIeee gsl_stats_short_ttest(const short data1[], const size_t stride1, const size_t n1, const short data2[], const size_t stride2, const size_t n2);

short gsl_stats_short_max (const short data[], const size_t stride, const size_t n);
short gsl_stats_short_min (const short data[], const size_t stride, const size_t n);
void gsl_stats_short_minmax (short * min, short * max, const short data[], const size_t stride, const size_t n);

size_t gsl_stats_short_max_index (const short data[], const size_t stride, const size_t n);
size_t gsl_stats_short_min_index (const short data[], const size_t stride, const size_t n);
void gsl_stats_short_minmax_index (size_t * min_index, size_t * max_index, const short data[], const size_t stride, const size_t n);

MpIeee gsl_stats_short_median_from_sorted_data(const short sorted_data[], const size_t stride, const size_t n) ;
MpIeee gsl_stats_short_quantile_from_sorted_data(const short sorted_data[], const size_t stride, const size_t n, const MpIeee f) ;

__END_DECLS

#endif /* __GSL_STATISTICS_SHORT_H__ */
