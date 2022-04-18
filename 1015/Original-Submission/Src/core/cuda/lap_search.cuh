#pragma once

namespace lap
{
	namespace cuda
	{
		template <class SC>
		__device__ __forceinline__ void searchWriteShared(SC *b_min, int *b_jmin, int *b_colsol, SC t_min, int t_jmin, int t_colsol)
		{
			int bidx = threadIdx.x >> 5;
			b_min[bidx] = t_min;
			b_jmin[bidx] = t_jmin;
			b_colsol[bidx] = t_colsol;
		}

		template <class SC>
		__device__ __forceinline__ void searchReadShared(SC &t_min, int &t_jmin, int &t_colsol, SC *b_min, int *b_jmin, int *b_colsol)
		{
			t_min = b_min[threadIdx.x];
			t_jmin = b_jmin[threadIdx.x];
			t_colsol = b_colsol[threadIdx.x];
		}

		template <class SC>
		__device__ __forceinline__ void searchCombineMedium(SC &t_min, int &t_jmin, int &t_colsol, SC *b_min, int *b_jmin, int *b_colsol)
		{
			minWarpIndex(t_min, t_jmin, t_colsol);
			if ((threadIdx.x & 0x1f) == 0)
			{
				searchWriteShared(b_min, b_jmin, b_colsol, t_min, t_jmin, t_colsol);
			}
			__syncthreads();
			if (threadIdx.x < 8)
			{
				searchReadShared(t_min, t_jmin, t_colsol, b_min, b_jmin, b_colsol);
				minWarpIndex8(t_min, t_jmin, t_colsol);
			}
		}

		template <class SC>
		__device__ __forceinline__ void searchCombineLarge(SC &t_min, int &t_jmin, int &t_colsol, SC *b_min, int *b_jmin, int *b_colsol)
		{
			minWarpIndex(t_min, t_jmin, t_colsol);
			if ((threadIdx.x & 0x1f) == 0)
			{
				searchWriteShared(b_min, b_jmin, b_colsol, t_min, t_jmin, t_colsol);
			}
			__syncthreads();
			if (threadIdx.x < 32)
			{
				searchReadShared(t_min, t_jmin, t_colsol, b_min, b_jmin, b_colsol);
				minWarpIndex(t_min, t_jmin, t_colsol);
			}
		}

		template <class SC>
		__device__ __forceinline__ void searchWriteTemp(volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC t_min, int t_jmin, int t_colsol)
		{
			if (threadIdx.x == 0)
			{
				o_min[blockIdx.x] = t_min;
				o_jmin[blockIdx.x] = t_jmin;
				o_colsol[blockIdx.x] = t_colsol;
			}
		}

		template <class SC>
		__device__ __forceinline__ void searchReadTemp(SC &t_min, int &t_jmin, int &t_colsol, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC max, int dim2)
		{
			if (threadIdx.x < gridDim.x)
			{
				t_min = o_min[threadIdx.x];
				t_jmin = o_jmin[threadIdx.x];
				t_colsol = o_colsol[threadIdx.x];
			}
			else
			{
				t_min = max;
				t_jmin = dim2;
				t_colsol = 0;
			}
		}

		template <class SC>
		__device__ __forceinline__ void searchReadTempLarge(SC &t_min, int &t_jmin, int &t_colsol, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC max, int dim2)
		{
			if (threadIdx.x < gridDim.x)
			{
				t_min = o_min[threadIdx.x];
				t_jmin = o_jmin[threadIdx.x];
				t_colsol = o_colsol[threadIdx.x];
				// read additional values
				for (int i = threadIdx.x + blockDim.x; i < gridDim.x; i += blockDim.x)
				{
					SC c_min = o_min[threadIdx.x];
					int c_jmin = o_jmin[threadIdx.x];
					int c_colsol = o_colsol[threadIdx.x];
					if ((c_min < t_min) || ((c_min == t_min) && ((c_colsol < t_colsol) || ((c_colsol == t_colsol) && (c_jmin < t_jmin)))))
					{
						t_min = c_min;
						t_jmin = c_jmin;
						t_colsol = c_colsol;
					}
				}
			}
			else
			{
				t_min = max;
				t_jmin = dim2;
				t_colsol = 0;
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void searchWrite(MS *s, SC t_min, int t_jmin, int t_colsol)
		{
			if (threadIdx.x == 0)
			{
				s->min = t_min;
				s->jmin = t_jmin;
				s->colsol = t_colsol;
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void searchSmall(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC &t_min, int &t_jmin, int &t_colsol, SC max, int size, int dim2)
		{
			minWarpIndex(t_min, t_jmin, t_colsol);
			searchWriteTemp(o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol);

			if (semaphoreWarp(semaphore))
			{
				searchReadTemp(t_min, t_jmin, t_colsol, o_min, o_jmin, o_colsol, max, dim2);
				minWarpIndex(t_min, t_jmin, t_colsol);
				searchWrite(s, t_min, t_jmin, t_colsol);
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void searchMedium(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC &t_min, int &t_jmin, int &t_colsol, SC max, int size, int dim2)
		{
			__shared__ SC b_min[8];
			__shared__ int b_jmin[8], b_colsol[8];

			searchCombineMedium(t_min, t_jmin, t_colsol, b_min, b_jmin, b_colsol);
			searchWriteTemp(o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol);

			if (semaphoreBlock(semaphore))
			{
				searchReadTemp(t_min, t_jmin, t_colsol, o_min, o_jmin, o_colsol, max, dim2);
				searchCombineMedium(t_min, t_jmin, t_colsol, b_min, b_jmin, b_colsol);
				searchWrite(s, t_min, t_jmin, t_colsol);
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void searchLarge(MS *s, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC &t_min, int &t_jmin, int &t_colsol, SC max, int size, int dim2)
		{
			__shared__ SC b_min[32];
			__shared__ int b_jmin[32], b_colsol[32];

			searchCombineLarge(t_min, t_jmin, t_colsol, b_min, b_jmin, b_colsol);
			searchWriteTemp(o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol);

			if (semaphoreBlock(semaphore))
			{
				searchReadTempLarge(t_min, t_jmin, t_colsol, o_min, o_jmin, o_colsol, max, dim2);
				searchCombineLarge(t_min, t_jmin, t_colsol, b_min, b_jmin, b_colsol);
				searchWrite(s, t_min, t_jmin, t_colsol);
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void searchSmall(MS *s, volatile MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC &t_min, int &t_jmin, int &t_colsol, SC max, int size, int dim2)
		{
			minWarpIndex(t_min, t_jmin, t_colsol);
			searchWriteTemp(o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol);

			if (semaphoreWarp(semaphore))
			{
				searchReadTemp(t_min, t_jmin, t_colsol, o_min, o_jmin, o_colsol, max, dim2);
				minWarpIndex(t_min, t_jmin, t_colsol);
				searchWrite(s, t_min, t_jmin, t_colsol);
				if (threadIdx.x == 0) s2->data_valid = 0;
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void searchMedium(MS *s, volatile MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC &t_min, int &t_jmin, int &t_colsol, SC max, int size, int dim2)
		{
			__shared__ SC b_min[8];
			__shared__ int b_jmin[8], b_colsol[8];

			searchCombineMedium(t_min, t_jmin, t_colsol, b_min, b_jmin, b_colsol);
			searchWriteTemp(o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol);

			if (semaphoreBlock(semaphore))
			{
				searchReadTemp(t_min, t_jmin, t_colsol, o_min, o_jmin, o_colsol, max, dim2);
				searchCombineMedium(t_min, t_jmin, t_colsol, b_min, b_jmin, b_colsol);
				searchWrite(s, t_min, t_jmin, t_colsol);
				if (threadIdx.x == 0) s2->data_valid = 0;
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void searchLarge(MS *s, volatile MS *s2, unsigned int *semaphore, volatile SC *o_min, volatile int *o_jmin, volatile int *o_colsol, SC &t_min, int &t_jmin, int &t_colsol, SC max, int size, int dim2)
		{
			__shared__ SC b_min[32];
			__shared__ int b_jmin[32], b_colsol[32];

			searchCombineLarge(t_min, t_jmin, t_colsol, b_min, b_jmin, b_colsol);
			searchWriteTemp(o_min, o_jmin, o_colsol, t_min, t_jmin, t_colsol);

			if (semaphoreBlock(semaphore))
			{
				searchReadTempLarge(t_min, t_jmin, t_colsol, o_min, o_jmin, o_colsol, max, dim2);
				searchCombineLarge(t_min, t_jmin, t_colsol, b_min, b_jmin, b_colsol);
				searchWrite(s, t_min, t_jmin, t_colsol);
				if (threadIdx.x == 0) s2->data_valid = 0;
			}
		}
	}
}
