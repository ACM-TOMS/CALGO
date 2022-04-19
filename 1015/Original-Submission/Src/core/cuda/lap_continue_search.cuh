#pragma once
#include "lap_search.cuh"

namespace lap
{
	namespace cuda
	{
		template <class SC, class TC>
		__device__ __forceinline__ void continueSearchJMinMinRead(SC &t_min, int &t_jmin, int &t_colsol, char *colactive, int *colsol, int *pred, SC *d, int j, int i, TC *tt, SC *v, SC min, int jmin, SC tt_jmin, SC v_jmin, SC max, int size, int dim2)
		{
			t_min = max;
			t_jmin = dim2;
			t_colsol = 0;

			if (j < size)
			{
				if (j == jmin) colactive[jmin] = 0;
				else if (colactive[j] != 0)
				{
					SC h = d[j];
					SC v2 = ((SC)tt[j] - tt_jmin) - (v[j] - v_jmin) + min;

					bool is_smaller = (v2 < h);
					if (is_smaller)
					{
						pred[j] = i;
						d[j] = h = v2;
					}

					t_min = h;
					t_jmin = j;
					t_colsol = colsol[j];
				}
			}
		}

		template <class SC>
		__device__ __forceinline__ void continueSearchJMinMinRead(SC &t_min, int &t_jmin, int &t_colsol, char *colactive, int *colsol, int *pred, SC *d, int j, int i, SC *v, SC min, int jmin, SC v_jmin, SC max, int size, int dim, int dim2)
		{
			t_min = max;
			t_jmin = dim2;
			t_colsol = 0;

			if (j < size)
			{
				// ignore any columns assigned to virtual rows
				if ((j == jmin) || ((colsol[j] >= dim) && (d[j] <= min))) colactive[j] = 0;
				else if (colactive[j] != 0)
				{
					SC h = d[j];
					SC v2 = -(v[j] - v_jmin) + min;

					bool is_smaller = (v2 < h);
					if (is_smaller)
					{
						pred[j] = i;
						d[j] = h = v2;
					}

					t_min = h;
					t_jmin = j;
					t_colsol = colsol[j];
				}
			}
		}

		template <class SC, class TC>
		__device__ __forceinline__ void continueSearchMinRead(SC &t_min, int &t_jmin, int &t_colsol, char *colactive, int *colsol, int *pred, SC *d, int j, int i, TC *tt, SC *v, SC min, SC tt_jmin, SC v_jmin, int jmin, SC max, int size, int dim2)
		{
			t_min = max;
			t_jmin = dim2;
			t_colsol = 0;

			if (j < size)
			{
				if (j == jmin) colactive[j] = 0;
				else if (colactive[j] != 0)
				{
					SC h = d[j];
					SC v2 = ((SC)tt[j] - tt_jmin) - (v[j] - v_jmin) + min;

					bool is_smaller = (v2 < h);
					if (is_smaller)
					{
						pred[j] = i;
						d[j] = h = v2;
					}

					t_min = h;
					t_jmin = j;
					t_colsol = colsol[j];
				}
			}
		}

		template <class SC>
		__device__ __forceinline__ void continueSearchMinRead(SC &t_min, int &t_jmin, int &t_colsol, char *colactive, int *colsol, int *pred, SC *d, int j, int i, SC *v, SC min, SC v_jmin, int jmin, SC max, int size, int dim, int dim2)
		{
			t_min = max;
			t_jmin = dim2;
			t_colsol = 0;

			if (j < size)
			{
				// ignore any columns assigned to virtual rows
				if ((j == jmin) || ((colsol[j] >= dim) && (d[j] <= min))) colactive[j] = 0;
				else if (colactive[j] != 0)
				{
					SC h = d[j];
					SC v2 = -(v[j] - v_jmin) + min;

					bool is_smaller = (v2 < h);
					if (is_smaller)
					{
						pred[j] = i;
						d[j] = h = v2;
					}

					t_min = h;
					t_jmin = j;
					t_colsol = colsol[j];
				}
			}
		}

		template <class MS, class SC, class TC>
		__global__ void continueSearchJMinMinSmall_kernel(MS *s, unsigned int *semaphore, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int i, int jmin, SC min, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_tt_jmin, b_v_jmin;
			if (threadIdx.x == 0)
			{
				b_tt_jmin = (SC)tt[jmin];
				b_v_jmin = v[jmin];
			}
			__syncthreads();
			SC tt_jmin = b_tt_jmin;
			SC v_jmin = b_v_jmin;

			SC t_min;
			int t_jmin, t_colsol;

			continueSearchJMinMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, tt, v, min, jmin, tt_jmin, v_jmin, max, size, dim2);
			searchSmall(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void continueSearchJMinMinMedium_kernel(MS *s, unsigned int *semaphore, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int i, int jmin, SC min, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_tt_jmin, b_v_jmin;
			if (threadIdx.x == 0)
			{
				b_tt_jmin = (SC)tt[jmin];
				b_v_jmin = v[jmin];
			}
			__syncthreads();
			SC tt_jmin = b_tt_jmin;
			SC v_jmin = b_v_jmin;

			SC t_min;
			int t_jmin, t_colsol;

			continueSearchJMinMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, tt, v, min, jmin, tt_jmin, v_jmin, max, size, dim2);
			searchMedium(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void continueSearchJMinMinLarge_kernel(MS *s, unsigned int *semaphore, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int i, int jmin, SC min, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_tt_jmin, b_v_jmin;
			if (threadIdx.x == 0)
			{
				b_tt_jmin = (SC)tt[jmin];
				b_v_jmin = v[jmin];
			}
			__syncthreads();
			SC tt_jmin = b_tt_jmin;
			SC v_jmin = b_v_jmin;

			SC t_min;
			int t_jmin, t_colsol;

			continueSearchJMinMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, tt, v, min, jmin, tt_jmin, v_jmin, max, size, dim2);
			searchLarge(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void continueSearchJMinMinSmallVirtual_kernel(MS *s, unsigned int *semaphore, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int i, int jmin, SC min, SC max, int size, int dim, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_v_jmin;
			if (threadIdx.x == 0)
			{
				b_v_jmin = v[jmin];
			}
			__syncthreads();
			SC v_jmin = b_v_jmin;

			SC t_min;
			int t_jmin, t_colsol;

			continueSearchJMinMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, v, min, jmin, v_jmin, max, size, dim, dim2);
			searchSmall(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void continueSearchJMinMinMediumVirtual_kernel(MS *s, unsigned int *semaphore, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int i, int jmin, SC min, SC max, int size, int dim, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_v_jmin;
			if (threadIdx.x == 0)
			{
				b_v_jmin = v[jmin];
			}
			__syncthreads();
			SC v_jmin = b_v_jmin;

			SC t_min;
			int t_jmin, t_colsol;

			continueSearchJMinMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, v, min, jmin, v_jmin, max, size, dim, dim2);
			searchMedium(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void continueSearchJMinMinLargeVirtual_kernel(MS *s, unsigned int *semaphore, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int i, int jmin, SC min, SC max, int size, int dim, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_v_jmin;
			if (threadIdx.x == 0)
			{
				b_v_jmin = v[jmin];
			}
			__syncthreads();
			SC v_jmin = b_v_jmin;

			SC t_min;
			int t_jmin, t_colsol;

			continueSearchJMinMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, v, min, jmin, v_jmin, max, size, dim, dim2);
			searchLarge(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__device__ __forceinline__ void continueSearchMinPeerExchange(volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *v, TC *tt, SC *v2, TC *tt2, SC &tt_jmin, SC &v_jmin, int jmin, int size)
		{
			__shared__ SC b_tt_jmin, b_v_jmin;

			if ((jmin >= 0) && (jmin < size))
			{
				if (threadIdx.x == 0)
				{
					b_tt_jmin = tt_jmin = (SC)tt[jmin];
					b_v_jmin = v_jmin = v[jmin];
					int sem = atomicInc(semaphore, gridDim.x - 1);
					if (sem == 0) data_valid[0] = 1;
				}
			}
			else
			{
				if (threadIdx.x == 0)
				{
					int sem = atomicInc(semaphore, gridDim.x - 1);
					if (sem == 0)
					{
						while (data_valid[0] == 0) {}
						__threadfence();
						b_tt_jmin = s2->min = (SC)tt2[0];
						b_v_jmin = s2->max = v2[0];
						__threadfence();
						s2->data_valid = 1;
					}
					else
					{
						while (s2->data_valid == 0) {}
						__threadfence();
						b_tt_jmin = s2->min;
						b_v_jmin = s2->max;
					}
				}
			}
			__syncthreads();
			tt_jmin = b_tt_jmin;
			v_jmin = b_v_jmin;
		}

		template <class MS, class SC, class TC>
		__device__ __forceinline__ void continueSearchMinPeerExchangeSmall(volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *v, TC *tt, SC *v2, TC *tt2, SC &tt_jmin, SC &v_jmin, int jmin, int size)
		{
			if ((jmin >= 0) && (jmin < size))
			{
				if (threadIdx.x == 0)
				{
					tt_jmin = (SC)tt[jmin];
					v_jmin = v[jmin];
					int sem = atomicInc(semaphore, gridDim.x - 1);
					if (sem == 0) data_valid[0] = 1;
				}
			}
			else
			{
				if (threadIdx.x == 0)
				{
					int sem = atomicInc(semaphore, gridDim.x - 1);
					if (sem == 0)
					{
						while (data_valid[0] == 0) {}
						__threadfence_system();
						s2->min = tt_jmin = (SC)tt2[0];
						s2->max = v_jmin = v2[0];
						__threadfence();
						s2->data_valid = 1;
					}
					else
					{
						while (s2->data_valid == 0) {}
						__threadfence();
						tt_jmin = s2->min;
						v_jmin = s2->max;
					}
				}
			}
			tt_jmin = __shfl_sync(0xffffffff, tt_jmin, 0, 32);
			v_jmin = __shfl_sync(0xffffffff, v_jmin, 0, 32);
		}

		template <class SC>
		__device__ __forceinline__ void continueSearchMinPeerExchange(SC *v2, SC &v_jmin)
		{
			__shared__ SC b_v_jmin;

			if (threadIdx.x == 0) b_v_jmin = v2[0];
			__syncthreads();
			v_jmin = b_v_jmin;
		}

		template <class SC>
		__device__ __forceinline__ void continueSearchMinPeerExchangeSmall(SC *v2, SC &v_jmin)
		{
			if (threadIdx.x == 0) v_jmin = v2[0];
			v_jmin = __shfl_sync(0xffffffff, v_jmin, 0, 32);
		}

		template <class MS, class SC, class TC>
		__device__ __forceinline__ void continueSearchMinExchange(volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *v, TC *tt, volatile SC *v2, volatile SC *tt2, SC &tt_jmin, SC &v_jmin, int jmin, int size)
		{
			__shared__ SC b_tt_jmin, b_v_jmin;

			if ((jmin >= 0) && (jmin < size))
			{
				if (threadIdx.x == 0)
				{
					b_tt_jmin = tt_jmin = (SC)tt[jmin];
					b_v_jmin = v_jmin = v[jmin];
					int sem = atomicInc(semaphore, gridDim.x - 1);
					if (sem == 0)
					{
						tt2[0] = tt_jmin;
						v2[0] = v_jmin;
						__threadfence_system();
						data_valid[0] = 1;
					}
				}
			}
			else
			{
				if (threadIdx.x == 0)
				{
					int sem = atomicInc(semaphore, gridDim.x - 1);
					if (sem == 0)
					{
						while (data_valid[0] == 0) {}
						__threadfence_system();
						b_tt_jmin = s2->min = (SC)tt2[0];
						b_v_jmin = s2->max = v2[0];
						__threadfence();
						s2->data_valid = 1;
					}
					else
					{
						while (s2->data_valid == 0) {}
						__threadfence();
						b_tt_jmin = s2->min;
						b_v_jmin = s2->max;
					}
				}
			}
			__syncthreads();
			tt_jmin = b_tt_jmin;
			v_jmin = b_v_jmin;
		}

		template <class MS, class SC, class TC>
		__device__ __forceinline__ void continueSearchMinExchangeSmall(volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *v, TC *tt, volatile SC *v2, volatile SC *tt2, SC &tt_jmin, SC &v_jmin, int jmin, int size)
		{
			if ((jmin >= 0) && (jmin < size))
			{
				if (threadIdx.x == 0)
				{
					tt_jmin = (SC)tt[jmin];
					v_jmin = v[jmin];
					int sem = atomicInc(semaphore, gridDim.x - 1);
					if (sem == 0)
					{
						tt2[0] = tt_jmin;
						v2[0] = v_jmin;
						__threadfence_system();
						data_valid[0] = 1;
					}
				}
			}
			else
			{
				if (threadIdx.x == 0)
				{
					int sem = atomicInc(semaphore, gridDim.x - 1);
					if (sem == 0)
					{
						while (data_valid[0] == 0) {}
						__threadfence_system();
						s2->min = tt_jmin = (SC)tt2[0];
						s2->max = v_jmin = v2[0];
						__threadfence();
						s2->data_valid = 1;
					}
					else
					{
						while (s2->data_valid == 0) {}
						__threadfence();
						tt_jmin = s2->min;
						v_jmin = s2->max;
					}
				}
			}
			tt_jmin = __shfl_sync(0xffffffff, tt_jmin, 0, 32);
			v_jmin = __shfl_sync(0xffffffff, v_jmin, 0, 32);
		}

		template <class MS, class SC>
		__device__ __forceinline__ void continueSearchMinExchange(volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *v, volatile SC *v2, SC &v_jmin, int jmin, int size)
		{
			__shared__ SC b_v_jmin;

			if ((jmin >= 0) && (jmin < size))
			{
				if (threadIdx.x == 0)
				{
					b_v_jmin = v_jmin = v[jmin];
					int sem = atomicInc(semaphore, gridDim.x - 1);
					if (sem == 0)
					{
						v2[0] = v_jmin;
						__threadfence_system();
						data_valid[0] = 1;
					}
				}
			}
			else
			{
				if (threadIdx.x == 0)
				{
					int sem = atomicInc(semaphore, gridDim.x - 1);
					if (sem == 0)
					{
						while (data_valid[0] == 0) {}
						__threadfence_system();
						b_v_jmin = s2->max = v2[0];
						__threadfence();
						s2->data_valid = 1;
					}
					else
					{
						while (s2->data_valid == 0) {}
						__threadfence();
						b_v_jmin = s2->max;
					}
				}
			}
			__syncthreads();
			v_jmin = b_v_jmin;
		}

		template <class MS, class SC>
		__device__ __forceinline__ void continueSearchMinExchangeSmall(volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *v, volatile SC *v2, SC &v_jmin, int jmin, int size)
		{
			if ((jmin >= 0) && (jmin < size))
			{
				if (threadIdx.x == 0)
				{
					v_jmin = v[jmin];
					int sem = atomicInc(semaphore, gridDim.x - 1);
					if (sem == 0)
					{
						v2[0] = v_jmin;
						__threadfence_system();
						data_valid[0] = 1;
					}
				}
			}
			else
			{
				if (threadIdx.x == 0)
				{
					int sem = atomicInc(semaphore, gridDim.x - 1);
					if (sem == 0)
					{
						while (data_valid[0] == 0) {}
						__threadfence_system();
						s2->max = v_jmin = v2[0];
						__threadfence();
						s2->data_valid = 1;
					}
					else
					{
						while (s2->data_valid == 0) {}
						__threadfence();
						v_jmin = s2->max;
					}
				}
			}
			v_jmin = __shfl_sync(0xffffffff, v_jmin, 0, 32);
		}

		template <class MS, class SC>
		__global__ void continueSearchMinPeerSmall_kernel(MS *s, unsigned int *semaphore, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int i, SC *v2, int jmin, SC min, SC max, int size, int dim, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;
			SC v_jmin;

			continueSearchMinPeerExchangeSmall(v2, v_jmin);
			continueSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, v, min, v_jmin, jmin, max, size, dim, dim2);
			searchSmall(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void continueSearchMinPeerMedium_kernel(MS *s, unsigned int *semaphore, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int i, SC *v2, int jmin, SC min, SC max, int size, int dim, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;
			SC v_jmin;

			continueSearchMinPeerExchange(v2, v_jmin);
			continueSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, v, min, v_jmin, jmin, max, size, dim, dim2);
			searchMedium(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void continueSearchMinPeerLarge_kernel(MS *s, unsigned int *semaphore, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int i, SC *v2, int jmin, SC min, SC max, int size, int dim, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;
			SC v_jmin;

			continueSearchMinPeerExchange(v2, v_jmin);
			continueSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, v, min, v_jmin, jmin, max, size, dim, dim2);
			searchLarge(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC, class TC2>
		__global__ void continueSearchMinPeerSmall_kernel(MS *s, volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int i, TC2 *tt2, SC *v2, int jmin, SC min, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;
			SC tt_jmin, v_jmin;

			continueSearchMinPeerExchangeSmall(s2, semaphore + 1, data_valid, v, tt, v2, tt2, tt_jmin, v_jmin, jmin, size);
			continueSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, tt, v, min, tt_jmin, v_jmin, jmin, max, size, dim2);
			searchSmall(s, s2, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC, class TC2>
		__global__ void continueSearchMinPeerMedium_kernel(MS *s, volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int i, TC2 *tt2, SC *v2, int jmin, SC min, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;
			SC tt_jmin, v_jmin;

			continueSearchMinPeerExchange(s2, semaphore + 1, data_valid, v, tt, v2, tt2, tt_jmin, v_jmin, jmin, size);
			continueSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, tt, v, min, tt_jmin, v_jmin, jmin, max, size, dim2);
			searchMedium(s, s2, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC, class TC2>
		__global__ void continueSearchMinPeerLarge_kernel(MS *s, volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int i, TC2 *tt2, SC *v2, int jmin, SC min, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;
			SC tt_jmin, v_jmin;

			continueSearchMinPeerExchange(s2, semaphore + 1, data_valid, v, tt, v2, tt2, tt_jmin, v_jmin, jmin, size);
			continueSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, tt, v, min, tt_jmin, v_jmin, jmin, max, size, dim2);
			searchLarge(s, s2, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void continueSearchMinSmall_kernel(MS *s, volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int i, SC *v2, int jmin, SC min, SC max, int size, int dim, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;
			SC v_jmin;

			continueSearchMinExchangeSmall(s2, semaphore + 1, data_valid, v, v2, v_jmin, jmin, size);
			continueSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, v, min, v_jmin, jmin, max, size, dim, dim2);
			searchSmall(s, s2, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void continueSearchMinMedium_kernel(MS *s, volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int i, SC *v2, int jmin, SC min, SC max, int size, int dim, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;
			SC v_jmin;

			continueSearchMinExchange(s2, semaphore + 1, data_valid, v, v2, v_jmin, jmin, size);
			continueSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, v, min, v_jmin, jmin, max, size, dim, dim2);
			searchMedium(s, s2, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void continueSearchMinLarge_kernel(MS *s, volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int i, SC *v2, int jmin, SC min, SC max, int size, int dim, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;
			SC v_jmin;

			continueSearchMinExchange(s2, semaphore + 1, data_valid, v, v2, v_jmin, jmin, size);
			continueSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, v, min, v_jmin, jmin, max, size, dim, dim2);
			searchLarge(s, s2, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void continueSearchMinSmall_kernel(MS *s, volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int i, volatile SC *tt2, volatile SC *v2, int jmin, SC min, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;
			SC tt_jmin, v_jmin;

			continueSearchMinExchangeSmall(s2, semaphore + 1, data_valid, v, tt, v2, tt2, tt_jmin, v_jmin, jmin, size);
			continueSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, tt, v, min, tt_jmin, v_jmin, jmin, max, size, dim2);
			searchSmall(s, s2, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void continueSearchMinMedium_kernel(MS *s, volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int i, volatile SC *tt2, volatile SC *v2, int jmin, SC min, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;
			SC tt_jmin, v_jmin;

			continueSearchMinExchange(s2, semaphore + 1, data_valid, v, tt, v2, tt2, tt_jmin, v_jmin, jmin, size);
			continueSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, tt, v, min, tt_jmin, v_jmin, jmin, max, size, dim2);
			searchMedium(s, s2, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void continueSearchMinLarge_kernel(MS *s, volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, SC *o_min, int *o_jmin, int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int i, volatile SC *tt2, volatile SC *v2, int jmin, SC min, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;
			SC tt_jmin, v_jmin;

			continueSearchMinExchange(s2, semaphore + 1, data_valid, v, tt, v2, tt2, tt_jmin, v_jmin, jmin, size);
			continueSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, i, tt, v, min, tt_jmin, v_jmin, jmin, max, size, dim2);
			searchLarge(s, s2, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}
	}
}
