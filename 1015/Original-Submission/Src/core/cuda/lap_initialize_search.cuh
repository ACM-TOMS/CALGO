#pragma once
#include "lap_search.cuh"

namespace lap
{
	namespace cuda
	{
		template <class SC, class TC>
		__device__ __forceinline__ void initializeSearchMinRead(SC &t_min, int &t_jmin, int &t_colsol, char *colactive, int *colsol, int *pred, SC *d, int j, int f, TC *tt, SC *v, SC max, int size, int dim2)
		{
			t_min = max;
			t_jmin = dim2;
			t_colsol = 0;

			if (j < size)
			{
				colactive[j] = 1;
				pred[j] = f;

				d[j] = t_min = (SC)tt[j] - v[j];
				t_jmin = j;
				t_colsol = colsol[j];
			}
		}

		template <class SC, class TC>
		__device__ __forceinline__ void initializeSearchMinRead(SC &t_min, int &t_jmin, int &t_colsol, char *colactive, int *colsol, int *colsol_in, int *pred, SC *d, int j, int f, TC *tt, SC *v, SC max, int size, int dim2)
		{
			t_min = max;
			t_jmin = dim2;
			t_colsol = 0;

			if (j < size)
			{
				colactive[j] = 1;
				pred[j] = f;

				d[j] = t_min = (SC)tt[j] - v[j];
				t_jmin = j;
				t_colsol = colsol[j] = colsol_in[j];
			}
		}

		template <class SC>
		__device__ __forceinline__ void initializeSearchMinRead(SC &t_min, int &t_jmin, int &t_colsol, char *colactive, int *colsol, int *pred, SC *d, int j, int f, SC *v, SC max, int size, int dim2)
		{
			t_min = max;
			t_jmin = dim2;
			t_colsol = 0;

			if (j < size)
			{
				colactive[j] = 1;
				pred[j] = f;

				d[j] = t_min = -v[j];
				t_jmin = j;
				t_colsol = colsol[j];
			}
		}

		template <class SC>
		__device__ __forceinline__ void initializeSearchMinRead(SC &t_min, int &t_jmin, int &t_colsol, char *colactive, int *colsol, int *colsol_in, int *pred, SC *d, int j, int f, SC *v, SC max, int size, int dim2)
		{
			t_min = max;
			t_jmin = dim2;
			t_colsol = 0;

			if (j < size)
			{
				colactive[j] = 1;
				pred[j] = f;

				d[j] = t_min = -v[j];
				t_jmin = j;
				t_colsol = colsol[j] = colsol_in[j];
			}
		}

		// normal
		template <class MS, class SC, class TC>
		__global__ void initializeSearchMinSmall_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, f, tt, v, max, size, dim2);
			searchSmall(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void initializeSearchMinMedium_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, f, tt, v, max, size, dim2);
			searchMedium(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void initializeSearchMinLarge_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, f, tt, v, max, size, dim2);
			searchLarge(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		// copy colsol
		template <class MS, class SC, class TC>
		__global__ void initializeSearchMinSmallCopy_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *colsol_in, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, colsol_in, pred, d, j, f, tt, v, max, size, dim2);
			searchSmall(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void initializeSearchMinMediumCopy_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *colsol_in, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, colsol_in, pred, d, j, f, tt, v, max, size, dim2);
			searchMedium(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void initializeSearchMinLargeCopy_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *colsol_in, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, colsol_in, pred, d, j, f, tt, v, max, size, dim2);
			searchLarge(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		// virtual row
		template <class MS, class SC>
		__global__ void initializeSearchMinSmallVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, f, v, max, size, dim2);
			searchSmall(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void initializeSearchMinMediumVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, f, v, max, size, dim2);
			searchMedium(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void initializeSearchMinLargeVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, f, v, max, size, dim2);
			searchLarge(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		// copy colsol, virtual row
		template <class MS, class SC>
		__global__ void initializeSearchMinSmallVirtualCopy_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *colsol_in, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, colsol_in, pred, d, j, f, v, max, size, dim2);
			searchSmall(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void initializeSearchMinMediumVirtualCopy_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *colsol_in, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, colsol_in, pred, d, j, f, v, max, size, dim2);
			searchMedium(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void initializeSearchMinLargeVirtualCopy_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *colsol_in, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, colsol_in, pred, d, j, f, v, max, size, dim2);
			searchLarge(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		// normal + initialize second struct
		template <class MS, class SC, class TC>
		__global__ void initializeSearchMinSmallRemote_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j == 0) s2->data_valid = 0;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, f, tt, v, max, size, dim2);
			searchSmall(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void initializeSearchMinMediumRemote_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j == 0) s2->data_valid = 0;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, f, tt, v, max, size, dim2);
			searchMedium(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void initializeSearchMinLargeRemote_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j == 0) s2->data_valid = 0;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, f, tt, v, max, size, dim2);
			searchLarge(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		// copy colsol + initialize second struct
		template <class MS, class SC, class TC>
		__global__ void initializeSearchMinSmallCopyRemote_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *colsol_in, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j == 0) s2->data_valid = 0;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, colsol_in, pred, d, j, f, tt, v, max, size, dim2);
			searchSmall(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void initializeSearchMinMediumCopyRemote_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *colsol_in, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j == 0) s2->data_valid = 0;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, colsol_in, pred, d, j, f, tt, v, max, size, dim2);
			searchMedium(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC, class TC>
		__global__ void initializeSearchMinLargeCopyRemote_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, TC *tt, char *colactive, int *colsol, int *colsol_in, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j == 0) s2->data_valid = 0;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, colsol_in, pred, d, j, f, tt, v, max, size, dim2);
			searchLarge(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		// virtual row + initialize second struct
		template <class MS, class SC>
		__global__ void initializeSearchMinSmallVirtualRemote_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j == 0) s2->data_valid = 0;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, f, v, max, size, dim2);
			searchSmall(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void initializeSearchMinMediumVirtualRemote_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j == 0) s2->data_valid = 0;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, f, v, max, size, dim2);
			searchMedium(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void initializeSearchMinLargeVirtualRemote_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j == 0) s2->data_valid = 0;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, pred, d, j, f, v, max, size, dim2);
			searchLarge(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		// copy colsol, virtual row + initialize second struct
		template <class MS, class SC>
		__global__ void initializeSearchMinSmallVirtualCopyRemote_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *colsol_in, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j == 0) s2->data_valid = 0;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, colsol_in, pred, d, j, f, v, max, size, dim2);
			searchSmall(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void initializeSearchMinMediumVirtualCopyRemote_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *colsol_in, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j == 0) s2->data_valid = 0;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, colsol_in, pred, d, j, f, v, max, size, dim2);
			searchMedium(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}

		template <class MS, class SC>
		__global__ void initializeSearchMinLargeVirtualCopyRemote_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC *v, SC *d, char *colactive, int *colsol, int *colsol_in, int *pred, int f, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j == 0) s2->data_valid = 0;

			SC t_min;
			int t_jmin, t_colsol;

			initializeSearchMinRead(t_min, t_jmin, t_colsol, colactive, colsol, colsol_in, pred, d, j, f, v, max, size, dim2);
			searchLarge(s, semaphore, o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol, max, size, dim2);
		}
	}
}
