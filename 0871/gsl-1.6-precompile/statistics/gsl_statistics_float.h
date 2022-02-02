#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/gsl_statistics_float.h
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

#ifndef __GSL_STATISTICS_FLOAT_H__
#define __GSL_STATISTICS_FLOAT_H__

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

MpIeee gsl_stats_float_mean (const MpIeee data[], const size_t stride, const size_t n);
MpIeee gsl_stats_float_variance(const MpIeee data[], const size_t stride, const size_t n);
MpIeee gsl_stats_float_sd(const MpIeee data[], const size_t stride, const size_t n);
MpIeee gsl_stats_float_variance_with_fixed_mean(const MpIeee data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_float_sd_with_fixed_mean(const MpIeee data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_float_absdev(const MpIeee data[], const size_t stride, const size_t n);
MpIeee gsl_stats_float_skew(const MpIeee data[], const size_t stride, const size_t n);
MpIeee gsl_stats_float_kurtosis(const MpIeee data[], const size_t stride, const size_t n);
MpIeee gsl_stats_float_lag1_autocorrelation(const MpIeee data[], const size_t stride, const size_t n);

MpIeee gsl_stats_float_covariance(const MpIeee data1[], const size_t stride1,const MpIeee data2[], const size_t stride2, const size_t n);

MpIeee gsl_stats_float_variance_m(const MpIeee data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_float_sd_m(const MpIeee data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_float_absdev_m(const MpIeee data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_float_skew_m_sd(const MpIeee data[], const size_t stride, const size_t n, const MpIeee mean, const MpIeee sd);
MpIeee gsl_stats_float_kurtosis_m_sd(const MpIeee data[], const size_t stride, const size_t n, const MpIeee mean, const MpIeee sd);
MpIeee gsl_stats_float_lag1_autocorrelation_m(const MpIeee data[], const size_t stride, const size_t n, const MpIeee mean);

MpIeee gsl_stats_float_covariance_m(const MpIeee data1[], const size_t stride1,const MpIeee data2[], const size_t stride2, const size_t n, const MpIeee mean1, const MpIeee mean2);

/* DEFINED FOR FLOATING POINT TYPES ONLY */

MpIeee gsl_stats_float_wmean(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n);
MpIeee gsl_stats_float_wvariance(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n);
MpIeee gsl_stats_float_wsd(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n);
MpIeee gsl_stats_float_wvariance_with_fixed_mean(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_float_wsd_with_fixed_mean(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n, const MpIeee mean);
MpIeee gsl_stats_float_wabsdev(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n);
MpIeee gsl_stats_float_wskew(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n);
MpIeee gsl_stats_float_wkurtosis(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n);

MpIeee gsl_stats_float_wvariance_m(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n, const MpIeee wmean);
MpIeee gsl_stats_float_wsd_m(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n, const MpIeee wmean);
MpIeee gsl_stats_float_wabsdev_m(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n, const MpIeee wmean);
MpIeee gsl_stats_float_wskew_m_sd(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n, const MpIeee wmean, const MpIeee wsd);
MpIeee gsl_stats_float_wkurtosis_m_sd(const MpIeee w[], const size_t wstride, const MpIeee data[], const size_t stride, const size_t n, const MpIeee wmean, const MpIeee wsd);

/* END OF FLOATING POINT TYPES */

MpIeee gsl_stats_float_pvariance(const MpIeee data1[], const size_t stride1, const size_t n1, const MpIeee data2[], const size_t stride2, const size_t n2);
MpIeee gsl_stats_float_ttest(const MpIeee data1[], const size_t stride1, const size_t n1, const MpIeee data2[], const size_t stride2, const size_t n2);

MpIeee gsl_stats_float_max(const MpIeee data[], const size_t stride, const size_t n);
MpIeee gsl_stats_float_min(const MpIeee data[], const size_t stride, const size_t n);
void gsl_stats_float_minmax (MpIeee * min, MpIeee * max, const MpIeee data[], const size_t stride, const size_t n);

size_t gsl_stats_float_max_index (const MpIeee data[], const size_t stride, const size_t n);
size_t gsl_stats_float_min_index (const MpIeee data[], const size_t stride, const size_t n);
void gsl_stats_float_minmax_index (size_t * min_index, size_t * max_index, const MpIeee data[], const size_t stride, const size_t n);

MpIeee gsl_stats_float_median_from_sorted_data(const MpIeee sorted_data[], const size_t stride, const size_t n) ;
MpIeee gsl_stats_float_quantile_from_sorted_data(const MpIeee sorted_data[], const size_t stride, const size_t n, const MpIeee f) ;

__END_DECLS

#endif /* __GSL_STATISTICS_FLOAT_H__ */
