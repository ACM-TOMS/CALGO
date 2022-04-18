#pragma once

namespace lap
{
	namespace cuda
	{
		template <class SC, class TC>
		__device__ __forceinline__ void getMinMaxBestRead(SC &t_min_cost, SC &t_max_cost, SC &t_picked_cost, int &t_jmin, int i, int j, TC *tt, int *picked, SC min, SC max, int size, int dim2)
		{
			t_min_cost = max;
			t_max_cost = min;
			t_picked_cost = max;
			t_jmin = dim2;

			if (j < size)
			{
				SC t_cost = (SC)tt[j];
				t_min_cost = t_cost;
				if (i == j) t_max_cost = t_cost;
				if (picked[j] == 0)
				{
					t_jmin = j;
					t_picked_cost = t_cost;
				}
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinMaxBestRead(SC &t_min_cost, SC &t_max_cost, SC &t_picked_cost, int &t_jmin, int i, int j, int *picked, SC min, SC max, int size, int dim2)
		{
			t_min_cost = max;
			t_max_cost = min;
			t_picked_cost = max;
			t_jmin = dim2;

			if (j < size)
			{
				SC t_cost = SC(0);
				t_min_cost = t_cost;
				if (i == j) t_max_cost = t_cost;
				if (picked[j] == 0)
				{
					t_jmin = j;
					t_picked_cost = t_cost;
				}
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinMaxBestCombineSmall(SC &t_min_cost, SC &t_max_cost, SC &t_picked_cost, int &t_jmin)
		{
			minWarpIndex(t_picked_cost, t_jmin);
			minWarp(t_min_cost);
			maxWarp(t_max_cost);
		}

		template <class SC>
		__device__ __forceinline__ void getMinMaxBestCombineTiny(SC &t_min_cost, SC &t_max_cost, SC &t_picked_cost, int &t_jmin)
		{
			minWarpIndex8(t_picked_cost, t_jmin);
			minWarp8(t_min_cost);
			maxWarp8(t_max_cost);
		}

		template <class SC>
		__device__ __forceinline__ void getMinMaxBestWriteShared(SC *b_min_cost, SC *b_max_cost, SC *b_picked_cost, int *b_jmin, SC t_min_cost, SC t_max_cost, SC t_picked_cost, int t_jmin)
		{
			int bidx = threadIdx.x >> 5;
			b_min_cost[bidx] = t_min_cost;
			b_max_cost[bidx] = t_max_cost;
			b_picked_cost[bidx] = t_picked_cost;
			b_jmin[bidx] = t_jmin;
		}

		template <class SC>
		__device__ __forceinline__ void getMinMaxBestReadShared(SC &t_min_cost, SC &t_max_cost, SC &t_picked_cost, int &t_jmin, SC *b_min_cost, SC *b_max_cost, SC *b_picked_cost, int *b_jmin)
		{
			t_min_cost = b_min_cost[threadIdx.x];
			t_max_cost = b_max_cost[threadIdx.x];
			t_picked_cost = b_picked_cost[threadIdx.x];
			t_jmin = b_jmin[threadIdx.x];
		}

		template <class SC>
		__device__ __forceinline__ void getMinMaxBestCombineMedium(SC &t_min_cost, SC &t_max_cost, SC &t_picked_cost, int &t_jmin, SC *b_min_cost, SC *b_max_cost, SC *b_picked_cost, int *b_jmin)
		{
			getMinMaxBestCombineSmall(t_min_cost, t_max_cost, t_picked_cost, t_jmin);
			if ((threadIdx.x & 0x1f) == 0)
			{
				getMinMaxBestWriteShared(b_min_cost, b_max_cost, b_picked_cost, b_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);
			}
			__syncthreads();
			if (threadIdx.x < 8)
			{
				getMinMaxBestReadShared(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinMaxBestCombineTiny(t_min_cost, t_max_cost, t_picked_cost, t_jmin);
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinMaxBestCombineLarge(SC &t_min_cost, SC &t_max_cost, SC &t_picked_cost, int &t_jmin, SC *b_min_cost, SC *b_max_cost, SC *b_picked_cost, int *b_jmin)
		{
			getMinMaxBestCombineSmall(t_min_cost, t_max_cost, t_picked_cost, t_jmin);
			if ((threadIdx.x & 0x1f) == 0)
			{
				getMinMaxBestWriteShared(b_min_cost, b_max_cost, b_picked_cost, b_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);
			}
			__syncthreads();
			if (threadIdx.x < 32)
			{
				getMinMaxBestReadShared(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinMaxBestCombineSmall(t_min_cost, t_max_cost, t_picked_cost, t_jmin);
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinMaxBestWriteTemp(volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC t_min_cost, SC t_max_cost, SC t_picked_cost, int t_jmin)
		{
			if (threadIdx.x == 0)
			{
				o_min_cost[blockIdx.x] = t_min_cost;
				o_max_cost[blockIdx.x] = t_max_cost;
				o_picked_cost[blockIdx.x] = t_picked_cost;
				o_jmin[blockIdx.x] = t_jmin;
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinMaxBestReadTemp(SC &t_min_cost, SC &t_max_cost, SC &t_picked_cost, int &t_jmin, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC min, SC max, int dim2)
		{
			if (threadIdx.x < gridDim.x)
			{
				t_min_cost = o_min_cost[threadIdx.x];
				t_max_cost = o_max_cost[threadIdx.x];
				t_picked_cost = o_picked_cost[threadIdx.x];
				t_jmin = o_jmin[threadIdx.x];
			}
			else
			{
				t_min_cost = max;
				t_max_cost = min;
				t_picked_cost = max;
				t_jmin = dim2;
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinMaxBestReadTempLarge(SC &t_min_cost, SC &t_max_cost, SC &t_picked_cost, int &t_jmin, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, SC min, SC max, int dim2)
		{
			if (threadIdx.x < gridDim.x)
			{
				t_min_cost = o_min_cost[threadIdx.x];
				t_max_cost = o_max_cost[threadIdx.x];
				t_picked_cost = o_picked_cost[threadIdx.x];
				t_jmin = o_jmin[threadIdx.x];
				// read additional values
				for (int i = threadIdx.x + blockDim.x; i < gridDim.x; i += blockDim.x)
				{
					SC c_min_cost = o_min_cost[i];
					SC c_max_cost = o_max_cost[i];
					SC c_picked_cost = o_picked_cost[i];
					int c_jmin = o_jmin[i];
					if (c_min_cost < t_min_cost) t_min_cost = c_min_cost;
					if (c_max_cost > t_max_cost) t_max_cost = c_max_cost;
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
				t_max_cost = min;
				t_picked_cost = max;
				t_jmin = dim2;
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void getMinMaxBestWrite(MS *s, SC t_min_cost, SC t_max_cost, SC t_picked_cost, int t_jmin, int *data_valid)
		{
			if (threadIdx.x == 0)
			{
				s->min = t_min_cost;
				s->max = t_max_cost;
				s->picked = t_picked_cost;
				s->jmin = t_jmin;

				__threadfence_system();
				data_valid[0] = 1;
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void getMinMaxBestSingleWrite(MS *s, volatile SC *o_min_cost, volatile int *o_jmin, SC t_min_cost, SC t_max_cost, SC t_picked_cost, int t_jmin, int i)
		{
			if (threadIdx.x == 0)
			{
				s->min = t_min_cost;
				s->max = t_max_cost;
				s->picked = t_picked_cost;
				s->jmin = t_jmin;

				o_min_cost[0] = t_min_cost;
				o_jmin[0] = t_jmin;
			}
		}

		// 32 threads per block, up to 32 blocks, requires no shared memory and no thread synchronization
		template <class MS, class SC, class TC>
		__global__ void getMinMaxBestSmall_kernel(MS *s, MS *s2, unsigned int *semaphore, int * data_valid, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, int *picked, SC min, SC max, int i, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_max_cost, t_picked_cost;
			int t_jmin;

			getMinMaxBestRead(t_min_cost, t_max_cost, t_picked_cost, t_jmin, i - start, j, tt, picked, min, max, size, dim2);
			getMinMaxBestCombineSmall(t_min_cost, t_max_cost, t_picked_cost, t_jmin);
			getMinMaxBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);

			if (semaphoreWarp(semaphore))
			{
				getMinMaxBestReadTemp(t_min_cost, t_max_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, min, max, dim2);
				getMinMaxBestCombineSmall(t_min_cost, t_max_cost, t_picked_cost, t_jmin);
				getMinMaxBestWrite(s, t_min_cost, t_max_cost, t_picked_cost, t_jmin + start, data_valid);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		// 256 threads per block, up to 256 blocks
		template <class MS, class SC, class TC>
		__global__ void getMinMaxBestMedium_kernel(MS *s, MS *s2, unsigned int *semaphore, int * data_valid, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, int *picked, SC min, SC max, int i, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_max_cost[8], b_picked_cost[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_max_cost, t_picked_cost;
			int t_jmin;

			getMinMaxBestRead(t_min_cost, t_max_cost, t_picked_cost, t_jmin, i - start, j, tt, picked, min, max, size, dim2);
			getMinMaxBestCombineMedium(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinMaxBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinMaxBestReadTemp(t_min_cost, t_max_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, min, max, dim2);
				getMinMaxBestCombineMedium(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinMaxBestWrite(s, t_min_cost, t_max_cost, t_picked_cost, t_jmin + start, data_valid);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		// 1024 threads per block, can be more than 1024 blocks
		template <class MS, class SC, class TC>
		__global__ void getMinMaxBestLarge_kernel(MS *s, MS *s2, unsigned int *semaphore, int * data_valid, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, int *picked, SC min, SC max, int i, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_max_cost[32], b_picked_cost[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_max_cost, t_picked_cost;
			int t_jmin;

			getMinMaxBestRead(t_min_cost, t_max_cost, t_picked_cost, t_jmin, i - start, j, tt, picked, min, max, size, dim2);
			getMinMaxBestCombineLarge(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinMaxBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinMaxBestReadTempLarge(t_min_cost, t_max_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, min, max, dim2);
				getMinMaxBestCombineLarge(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinMaxBestWrite(s, t_min_cost, t_max_cost, t_picked_cost, t_jmin + start, data_valid);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		// 32 threads per block, up to 32 blocks, requires no shared memory and no thread synchronization
		template <class MS, class SC, class TC>
		__global__ void getMinMaxBestSingleSmall_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, int *picked, SC min, SC max, int i, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_max_cost, t_picked_cost;
			int t_jmin;

			getMinMaxBestRead(t_min_cost, t_max_cost, t_picked_cost, t_jmin, i, j, tt, picked, min, max, size, dim2);
			getMinMaxBestCombineSmall(t_min_cost, t_max_cost, t_picked_cost, t_jmin);
			getMinMaxBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);

			if (semaphoreWarp(semaphore))
			{
				getMinMaxBestReadTemp(t_min_cost, t_max_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, min, max, dim2);
				getMinMaxBestCombineSmall(t_min_cost, t_max_cost, t_picked_cost, t_jmin);
				getMinMaxBestSingleWrite(s, o_min_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin, i);
			}
		}

		// 256 threads per block, up to 256 blocks
		template <class MS, class SC, class TC>
		__global__ void getMinMaxBestSingleMedium_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, int *picked, SC min, SC max, int i, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_max_cost[8], b_picked_cost[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_max_cost, t_picked_cost;
			int t_jmin;

			getMinMaxBestRead(t_min_cost, t_max_cost, t_picked_cost, t_jmin, i, j, tt, picked, min, max, size, dim2);
			getMinMaxBestCombineMedium(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinMaxBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinMaxBestReadTemp(t_min_cost, t_max_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, min, max, dim2);
				getMinMaxBestCombineMedium(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinMaxBestSingleWrite(s, o_min_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin, i);
			}
		}

		// 1024 threads per block, can be more than 1024 blocks
		template <class MS, class SC, class TC>
		__global__ void getMinMaxBestSingleLarge_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, TC *tt, int *picked, SC min, SC max, int i, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_max_cost[32], b_picked_cost[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_max_cost, t_picked_cost;
			int t_jmin;

			getMinMaxBestRead(t_min_cost, t_max_cost, t_picked_cost, t_jmin, i, j, tt, picked, min, max, size, dim2);
			getMinMaxBestCombineLarge(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinMaxBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinMaxBestReadTempLarge(t_min_cost, t_max_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, min, max, dim2);
				getMinMaxBestCombineLarge(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinMaxBestSingleWrite(s, o_min_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin, i);
			}
		}

		// 32 threads per block, up to 32 blocks, requires no shared memory and no thread synchronization
		template <class MS, class SC>
		__global__ void getMinMaxBestSmallVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, int * data_valid, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, int *picked, SC min, SC max, int i, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_max_cost, t_picked_cost;
			int t_jmin;

			getMinMaxBestRead(t_min_cost, t_max_cost, t_picked_cost, t_jmin, i - start, j, picked, min, max, size, dim2);
			getMinMaxBestCombineSmall(t_min_cost, t_max_cost, t_picked_cost, t_jmin);
			getMinMaxBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);

			if (semaphoreWarp(semaphore))
			{
				getMinMaxBestReadTemp(t_min_cost, t_max_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, min, max, dim2);
				getMinMaxBestCombineSmall(t_min_cost, t_max_cost, t_picked_cost, t_jmin);
				getMinMaxBestWrite(s, t_min_cost, t_max_cost, t_picked_cost, t_jmin + start, data_valid);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		// 256 threads per block, up to 256 blocks
		template <class MS, class SC>
		__global__ void getMinMaxBestMediumVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, int * data_valid, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, int *picked, SC min, SC max, int i, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_max_cost[8], b_picked_cost[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_max_cost, t_picked_cost;
			int t_jmin;

			getMinMaxBestRead(t_min_cost, t_max_cost, t_picked_cost, t_jmin, i - start, j, picked, min, max, size, dim2);
			getMinMaxBestCombineMedium(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinMaxBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinMaxBestReadTemp(t_min_cost, t_max_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, min, max, dim2);
				getMinMaxBestCombineMedium(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinMaxBestWrite(s, t_min_cost, t_max_cost, t_picked_cost, t_jmin + start, data_valid);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		// 1024 threads per block, can be more than 1024 blocks
		template <class MS, class SC>
		__global__ void getMinMaxBestLargeVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, int * data_valid, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, int *picked, SC min, SC max, int i, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_max_cost[32], b_picked_cost[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_max_cost, t_picked_cost;
			int t_jmin;

			getMinMaxBestRead(t_min_cost, t_max_cost, t_picked_cost, t_jmin, i - start, j, picked, min, max, size, dim2);
			getMinMaxBestCombineLarge(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinMaxBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinMaxBestReadTempLarge(t_min_cost, t_max_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, min, max, dim2);
				getMinMaxBestCombineLarge(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinMaxBestWrite(s, t_min_cost, t_max_cost, t_picked_cost, t_jmin + start, data_valid);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		// 32 threads per block, up to 32 blocks, requires no shared memory and no thread synchronization
		template <class MS, class SC>
		__global__ void getMinMaxBestSingleSmallVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, int *picked, SC min, SC max, int i, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_max_cost, t_picked_cost;
			int t_jmin;

			getMinMaxBestRead(t_min_cost, t_max_cost, t_picked_cost, t_jmin, i, j, picked, min, max, size, dim2);
			getMinMaxBestCombineSmall(t_min_cost, t_max_cost, t_picked_cost, t_jmin);
			getMinMaxBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);

			if (semaphoreWarp(semaphore))
			{
				getMinMaxBestReadTemp(t_min_cost, t_max_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, min, max, dim2);
				getMinMaxBestCombineSmall(t_min_cost, t_max_cost, t_picked_cost, t_jmin);
				getMinMaxBestSingleWrite(s, o_min_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin, i);
			}
		}

		// 256 threads per block, up to 256 blocks
		template <class MS, class SC>
		__global__ void getMinMaxBestSingleMediumVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, int *picked, SC min, SC max, int i, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_max_cost[8], b_picked_cost[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_max_cost, t_picked_cost;
			int t_jmin;

			getMinMaxBestRead(t_min_cost, t_max_cost, t_picked_cost, t_jmin, i, j, picked, min, max, size, dim2);
			getMinMaxBestCombineMedium(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinMaxBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinMaxBestReadTemp(t_min_cost, t_max_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, min, max, dim2);
				getMinMaxBestCombineMedium(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinMaxBestSingleWrite(s, o_min_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin, i);
			}
		}

		// 1024 threads per block, can be more than 1024 blocks
		template <class MS, class SC>
		__global__ void getMinMaxBestSingleLargeVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_max_cost, volatile SC *o_picked_cost, volatile int *o_jmin, int *picked, SC min, SC max, int i, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_max_cost[32], b_picked_cost[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_max_cost, t_picked_cost;
			int t_jmin;

			getMinMaxBestRead(t_min_cost, t_max_cost, t_picked_cost, t_jmin, i, j, picked, min, max, size, dim2);
			getMinMaxBestCombineLarge(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
			getMinMaxBestWriteTemp(o_min_cost, o_max_cost, o_picked_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin);

			if (semaphoreBlock(semaphore))
			{
				getMinMaxBestReadTempLarge(t_min_cost, t_max_cost, t_picked_cost, t_jmin, o_min_cost, o_max_cost, o_picked_cost, o_jmin, min, max, dim2);
				getMinMaxBestCombineLarge(t_min_cost, t_max_cost, t_picked_cost, t_jmin, b_min_cost, b_max_cost, b_picked_cost, b_jmin);
				getMinMaxBestSingleWrite(s, o_min_cost, o_jmin, t_min_cost, t_max_cost, t_picked_cost, t_jmin, i);
			}
		}
	}
}
