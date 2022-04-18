#pragma once

namespace lap
{
	namespace cuda
	{
		template <class SC, class TC>
		__device__ __forceinline__ void getMinSecondBestRead(SC &t_min_cost, SC &t_second_cost, SC &t_picked_cost, int &t_jmin, int i, int j, TC *tt, SC *v, int *picked, SC max, int size, int dim2)
		{
			t_min_cost = max;
			t_second_cost = max;
			t_picked_cost = max;
			t_jmin = dim2;

			if (j < size)
			{
				SC t_cost = (SC)tt[j] - v[j];
				t_min_cost = t_cost;
				if (picked[j] == 0)
				{
					t_jmin = j;
					t_picked_cost = t_cost;
				}
			}
		}

		template <class SC, class TC>
		__device__ __forceinline__ void getMinSecondBestRead(SC &t_min_cost, SC &t_second_cost, SC &t_picked_cost, int &t_jmin, int i, int j, TC *tt, SC *v, int *picked, int last_picked, SC max, int size, int dim2)
		{
			t_min_cost = max;
			t_second_cost = max;
			t_picked_cost = max;
			t_jmin = dim2;

			if (j < size)
			{
				SC t_cost = (SC)tt[j] - v[j];
				t_min_cost = t_cost;
				if (j == last_picked)
				{
					picked[j] = 1;
				}
				else if (picked[j] == 0)
				{
					t_jmin = j;
					t_picked_cost = t_cost;
				}
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinSecondBestRead(SC &t_min_cost, SC &t_second_cost, SC &t_picked_cost, int &t_jmin, int i, int j, SC *v, int *picked, SC max, int size, int dim2)
		{
			t_min_cost = max;
			t_second_cost = max;
			t_picked_cost = max;
			t_jmin = dim2;

			if (j < size)
			{
				SC t_cost = -v[j];
				t_min_cost = t_cost;
				if (picked[j] == 0)
				{
					t_jmin = j;
					t_picked_cost = t_cost;
				}
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinSecondBestRead(SC &t_min_cost, SC &t_second_cost, SC &t_picked_cost, int &t_jmin, int i, int j, SC *v, int *picked, int last_picked, SC max, int size, int dim2)
		{
			t_min_cost = max;
			t_second_cost = max;
			t_picked_cost = max;
			t_jmin = dim2;

			if (j < size)
			{
				SC t_cost = -v[j];
				t_min_cost = t_cost;
				if (j == last_picked)
				{
					picked[j] = 1;
				}
				else if (picked[j] == 0)
				{
					t_jmin = j;
					t_picked_cost = t_cost;
				}
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinSecondBestCombineSmall(SC &t_min_cost, SC &t_second_cost, SC &t_picked_cost, int &t_jmin)
		{
			minWarpIndex(t_picked_cost, t_jmin);
			SC old_min_cost = t_min_cost;
			minWarp(t_min_cost);
			bool is_min = (t_min_cost == old_min_cost);
			int mask = __ballot_sync(0xffffffff, is_min);
			is_min &= (__popc(mask) == 1);
			if (!is_min) t_second_cost = old_min_cost;
			minWarp(t_second_cost);
		}

		template <class SC>
		__device__ __forceinline__ void getMinSecondBestCombineTiny(SC &t_min_cost, SC &t_second_cost, SC &t_picked_cost, int &t_jmin)
		{
			minWarpIndex8(t_picked_cost, t_jmin);
			SC old_min_cost = t_min_cost;
			minWarp8(t_min_cost);
			bool is_min = (t_min_cost == old_min_cost);
			int mask = __ballot_sync(0xff, is_min);
			is_min &= (__popc(mask) == 1);
			if (!is_min) t_second_cost = old_min_cost;
			minWarp8(t_second_cost);
		}

		template <class SC>
		__device__ __forceinline__ void getMinSecondBestWriteShared(SC *b_min_cost, SC *b_max_cost, SC *b_picked_cost, int *b_jmin, SC t_min_cost, SC t_second_cost, SC t_picked_cost, int t_jmin)
		{
			int bidx = threadIdx.x >> 5;
			b_min_cost[bidx] = t_min_cost;
			b_max_cost[bidx] = t_second_cost;
			b_picked_cost[bidx] = t_picked_cost;
			b_jmin[bidx] = t_jmin;
		}

		template <class SC>
		__device__ __forceinline__ void getMinSecondBestReadShared(SC &t_min_cost, SC &t_second_cost, SC &t_picked_cost, int &t_jmin, SC *b_min_cost, SC *b_max_cost, SC *b_picked_cost, int *b_jmin)
		{
			t_min_cost = b_min_cost[threadIdx.x];
			t_second_cost = b_max_cost[threadIdx.x];
			t_picked_cost = b_picked_cost[threadIdx.x];
			t_jmin = b_jmin[threadIdx.x];
		}

		template <class SC>
		__device__ __forceinline__ void getMinSecondBestCombineMedium(SC &t_min_cost, SC &t_second_cost, SC &t_picked_cost, int &t_jmin, SC *b_min_cost, SC *b_max_cost, SC *b_picked_cost, int *b_jmin)
		{
			getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			if ((threadIdx.x & 0x1f) == 0)
			{
				getMinSecondBestWriteShared(b_min_cost, b_max_cost, b_picked_cost, b_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			}
			__syncthreads();
			if (threadIdx.x < 8)
			{
				getMinSecondBestReadShared(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestCombineTiny(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinSecondBestCombineLarge(SC &t_min_cost, SC &t_second_cost, SC &t_picked_cost, int &t_jmin, SC *b_min_cost, SC *b_max_cost, SC *b_picked_cost, int *b_jmin)
		{
			getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			if ((threadIdx.x & 0x1f) == 0)
			{
				getMinSecondBestWriteShared(b_min_cost, b_max_cost, b_picked_cost, b_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			}
			__syncthreads();
			if (threadIdx.x < 32)
			{
				getMinSecondBestReadShared(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinSecondBestWriteTemp(volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC t_min_cost, SC t_second_cost, SC t_picked_cost, int t_jmin)
		{
			if (threadIdx.x == 0)
			{
				o_min_cost[blockIdx.x] = t_min_cost;
				o_max_cost[blockIdx.x] = t_second_cost;
				o_picked_cost[blockIdx.x] = t_picked_cost;
				o_jmin[blockIdx.x] = t_jmin;
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinSecondBestReadTemp(SC &t_min_cost, SC &t_second_cost, SC &t_picked_cost, int &t_jmin, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC max, int dim2)
		{
			if (threadIdx.x < gridDim.x)
			{
				t_min_cost = o_min_cost[threadIdx.x];
				t_second_cost = o_max_cost[threadIdx.x];
				t_picked_cost = o_picked_cost[threadIdx.x];
				t_jmin = o_jmin[threadIdx.x];
			}
			else
			{
				t_min_cost = max;
				t_second_cost = max;
				t_picked_cost = max;
				t_jmin = dim2;
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinSecondBestReadTempLarge(SC &t_min_cost, SC &t_second_cost, SC &t_picked_cost, int &t_jmin, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC max, int dim2)
		{
			if (threadIdx.x < gridDim.x)
			{
				t_min_cost = o_min_cost[threadIdx.x];
				t_second_cost = o_max_cost[threadIdx.x];
				t_picked_cost = o_picked_cost[threadIdx.x];
				t_jmin = o_jmin[threadIdx.x];
				// read additional values
				for (int i = threadIdx.x + blockDim.x; i < gridDim.x; i += blockDim.x)
				{
					SC c_min_cost = o_min_cost[i];
					SC c_second_cost = o_max_cost[i];
					SC c_picked_cost = o_picked_cost[i];
					int c_jmin = o_jmin[i];
					if (c_min_cost < t_min_cost)
					{
						if (t_min_cost < c_second_cost) t_second_cost = t_min_cost;
						else t_second_cost = c_second_cost;
						t_min_cost = c_min_cost;
					}
					else if (c_min_cost < t_second_cost) t_second_cost = c_min_cost;
					if ((c_picked_cost < t_picked_cost) || ((c_picked_cost == t_picked_cost) && (c_jmin < t_jmin)))
					{
						t_jmin = c_jmin;
						t_picked_cost = c_picked_cost;
					}
				}
			}
			else
			{
				t_min_cost = max;
				t_second_cost = max;
				t_picked_cost = max;
				t_jmin = dim2;
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void getMinSecondBestWrite(MS *s, SC t_min_cost, SC t_second_cost, SC t_picked_cost, int t_jmin, int start, SC *v, SC max, int dim2)
		{
			if (threadIdx.x == 0)
			{
				s->min = t_min_cost;
				s->max = t_second_cost;
				s->picked = t_picked_cost;
				if (t_jmin < dim2)
				{
					s->v_jmin = v[t_jmin];
				}
				else
				{
					s->v_jmin = max;
				}
				__threadfence_system();
				s->jmin = t_jmin + start;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinSecondBestSmall_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, SC *v, int *picked, SC max, int i, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, tt, v, picked, max, size, dim2);
			getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreWarp(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinSecondBestMedium_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, SC *v, int *picked, SC max, int i, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_max_cost[8], b_picked_cost[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, tt, v, picked, max, size, dim2);
			getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinSecondBestLarge_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, SC *v, int *picked, SC max, int i, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_max_cost[32], b_picked_cost[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, tt, v, picked, max, size, dim2);
			getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTempLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinSecondBestSmallSmall_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, SC *v, int *picked, volatile MS *s_old, SC max, int i, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			int last_picked = getLastPickedSmall(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, tt, v, picked, last_picked, max, size, dim2);
			getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreWarp(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0) s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinSecondBestSmallMedium_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, SC *v, int *picked, volatile MS *s_old, SC max, int i, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_max_cost[8], b_picked_cost[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			int last_picked = getLastPicked(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, tt, v, picked, last_picked, max, size, dim2);
			getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinSecondBestSmallLarge_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, SC *v, int *picked, volatile MS *s_old, SC max, int i, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_max_cost[32], b_picked_cost[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			int last_picked = getLastPicked(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, tt, v, picked, last_picked, max, size, dim2);
			getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTempLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinSecondBestLargeSmall_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, SC *v, int *picked, volatile MS *s_old, SC max, int i, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			int last_picked = getLastPickedLarge(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, tt, v, picked, last_picked, max, size, dim2);
			getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreWarp(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinSecondBestLargeMedium_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, SC *v, int *picked, volatile MS *s_old, SC max, int i, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_max_cost[8], b_picked_cost[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			int last_picked = getLastPickedLarge(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, tt, v, picked, last_picked, max, size, dim2);
			getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinSecondBestLargeLarge_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, SC *v, int *picked, volatile MS *s_old, SC max, int i, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_max_cost[32], b_picked_cost[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			int last_picked = getLastPickedLarge(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, tt, v, picked, last_picked, max, size, dim2);
			getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTempLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinSecondBestSingleSmall_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, SC *v, int *picked, SC max, int i, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, tt, v, picked, max, size, dim2);
			getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreWarp(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, 0, v, max, dim2);
				if (threadIdx.x == 0) picked[t_jmin] = 1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinSecondBestSingleMedium_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, SC *v, int *picked, SC max, int i, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_max_cost[8], b_picked_cost[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, tt, v, picked, max, size, dim2);
			getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, 0, v, max, dim2);
				if (threadIdx.x == 0) picked[t_jmin] = 1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinSecondBestSingleLarge_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, SC *v, int *picked, SC max, int i, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_max_cost[32], b_picked_cost[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, tt, v, picked, max, size, dim2);
			getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTempLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, 0, v, max, dim2);
				if (threadIdx.x == 0) picked[t_jmin] = 1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinSecondBestSmallVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC *v, int *picked, SC max, int i, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, v, picked, max, size, dim2);
			getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreWarp(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinSecondBestMediumVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC *v, int *picked, SC max, int i, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_max_cost[8], b_picked_cost[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, v, picked, max, size, dim2);
			getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinSecondBestLargeVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC *v, int *picked, SC max, int i, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_max_cost[32], b_picked_cost[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, v, picked, max, size, dim2);
			getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTempLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinSecondBestSmallSmallVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC *v, int *picked, volatile MS *s_old, SC max, int i, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			int last_picked = getLastPickedSmall(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, v, picked, last_picked, max, size, dim2);
			getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreWarp(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0) s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinSecondBestSmallMediumVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC *v, int *picked, volatile MS *s_old, SC max, int i, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_max_cost[8], b_picked_cost[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			int last_picked = getLastPicked(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, v, picked, last_picked, max, size, dim2);
			getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinSecondBestSmallLargeVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC *v, int *picked, volatile MS *s_old, SC max, int i, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_max_cost[32], b_picked_cost[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			int last_picked = getLastPicked(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, v, picked, last_picked, max, size, dim2);
			getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTempLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinSecondBestLargeSmallVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC *v, int *picked, volatile MS *s_old, SC max, int i, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			int last_picked = getLastPickedLarge(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, v, picked, last_picked, max, size, dim2);
			getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreWarp(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinSecondBestLargeMediumVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC *v, int *picked, volatile MS *s_old, SC max, int i, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_max_cost[8], b_picked_cost[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			int last_picked = getLastPickedLarge(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, v, picked, last_picked, max, size, dim2);
			getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinSecondBestLargeLargeVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC *v, int *picked, volatile MS *s_old, SC max, int i, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_max_cost[32], b_picked_cost[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			int last_picked = getLastPickedLarge(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, v, picked, last_picked, max, size, dim2);
			getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTempLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, start, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinSecondBestSingleSmallVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC *v, int *picked, SC max, int i, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, v, picked, max, size, dim2);
			getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreWarp(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineSmall(t_min_cost, t_second_cost, t_picked_cost, t_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, 0, v, max, dim2);
				if (threadIdx.x == 0) picked[t_jmin] = 1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinSecondBestSingleMediumVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC *v, int *picked, SC max, int i, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_max_cost[8], b_picked_cost[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, v, picked, max, size, dim2);
			getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTemp(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineMedium(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, 0, v, max, dim2);
				if (threadIdx.x == 0) picked[t_jmin] = 1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinSecondBestSingleLargeVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC *v, int *picked, SC max, int i, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_max_cost[32], b_picked_cost[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_second_cost, t_picked_cost;
			int t_jmin;

			getMinSecondBestRead(t_min_cost, t_second_cost, t_picked_cost, t_jmin, i, j, v, picked, max, size, dim2);
			getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinSecondBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_second_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinSecondBestReadTempLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, max, dim2);
				getMinSecondBestCombineLarge(t_min_cost, t_second_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinSecondBestWrite(s, t_min_cost, t_second_cost, t_picked_cost, t_jmin, 0, v, max, dim2);
				if (threadIdx.x == 0) picked[t_jmin] = 1;
			}
		}
	}
}
