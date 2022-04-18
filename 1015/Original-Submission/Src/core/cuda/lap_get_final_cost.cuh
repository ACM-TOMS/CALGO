#pragma once

namespace lap
{
	namespace cuda
	{
		template <class SC, class TC>
		__device__ __forceinline__ void getFinalCostRead(SC &t_min_cost, SC &t_picked_cost, SC &t_picked_v, int j, TC *tt, SC *v, int j_picked, SC max, int size)
		{
			t_min_cost = max;
			t_picked_cost = max;
			t_picked_v = max;

			if (j < size)
			{
				SC t_cost = (SC)tt[j] - v[j];
				t_min_cost = t_cost;
				if (j == j_picked)
				{
					t_picked_cost = (SC)tt[j];
					t_picked_v = v[j];
				}
			}
		}

		template <class SC>
		__device__ __forceinline__ void getFinalCostRead(SC &t_min_cost, SC &t_picked_cost, SC &t_picked_v, int j, SC *v, int j_picked, SC max, int size)
		{
			t_min_cost = max;
			t_picked_cost = max;
			t_picked_v = max;

			if (j < size)
			{
				SC t_cost = -v[j];
				t_min_cost = t_cost;
				if (j == j_picked)
				{
					t_picked_cost = SC(0);
					t_picked_v = v[j];
				}
			}
		}

		template <class SC>
		__device__ __forceinline__ void getFinalCostCombineSmall(SC &t_min_cost, SC &t_picked_cost, SC &t_picked_v)
		{
			minWarp(t_min_cost);
			minWarp(t_picked_cost);
			minWarp(t_picked_v);
		}

		template <class SC>
		__device__ __forceinline__ void getFinalCostCombineTiny(SC &t_min_cost, SC &t_picked_cost, SC &t_picked_v)
		{
			minWarp8(t_min_cost);
			minWarp8(t_picked_cost);
			minWarp8(t_picked_v);
		}

		template <class SC>
		__device__ __forceinline__ void getFinalCostWriteShared(SC *b_min_cost, SC *b_picked_cost, SC *b_picked_v, SC t_min_cost, SC t_picked_cost, SC t_picked_v)
		{
			int bidx = threadIdx.x >> 5;
			b_min_cost[bidx] = t_min_cost;
			b_picked_cost[bidx] = t_picked_cost;
			b_picked_v[bidx] = t_picked_v;
		}

		template <class SC>
		__device__ __forceinline__ void getFinalCostReadShared(SC &t_min_cost, SC &t_picked_cost, SC &t_picked_v, SC *b_min_cost, SC *b_picked_cost, SC *b_picked_v)
		{
			t_min_cost = b_min_cost[threadIdx.x];
			t_picked_cost = b_picked_cost[threadIdx.x];
			t_picked_v = b_picked_v[threadIdx.x];
		}

		template <class SC>
		__device__ __forceinline__ void getFinalCostCombineMedium(SC &t_min_cost, SC &t_picked_cost, SC &t_picked_v, SC *b_min_cost, SC *b_picked_cost, SC *b_picked_v)
		{
			getFinalCostCombineSmall(t_min_cost, t_picked_cost, t_picked_v);
			if ((threadIdx.x & 0x1f) == 0)
			{
				getFinalCostWriteShared(b_min_cost, b_picked_cost, b_picked_v, t_min_cost, t_picked_cost, t_picked_v);
			}
			__syncthreads();
			if (threadIdx.x < 8)
			{
				getFinalCostReadShared(t_min_cost, t_picked_cost, t_picked_v, b_min_cost, b_picked_cost, b_picked_v);
				getFinalCostCombineTiny(t_min_cost, t_picked_cost, t_picked_v);
			}
		}

		template <class SC>
		__device__ __forceinline__ void getFinalCostCombineLarge(SC &t_min_cost, SC &t_picked_cost, SC &t_picked_v, SC *b_min_cost, SC *b_picked_cost, SC *b_picked_v)
		{
			getFinalCostCombineSmall(t_min_cost, t_picked_cost, t_picked_v);
			if ((threadIdx.x & 0x1f) == 0)
			{
				getFinalCostWriteShared(b_min_cost, b_picked_cost, b_picked_v, t_min_cost, t_picked_cost, t_picked_v);
			}
			__syncthreads();
			if (threadIdx.x < 32)
			{
				getFinalCostReadShared(t_min_cost, t_picked_cost, t_picked_v, b_min_cost, b_picked_cost, b_picked_v);
				getFinalCostCombineSmall(t_min_cost, t_picked_cost, t_picked_v);
			}
		}

		template <class SC>
		__device__ __forceinline__ void getFinalCostWriteTemp(volatile SC *o_min_cost, volatile SC *o_picked_cost, volatile SC *o_picked_v, SC t_min_cost, SC t_picked_cost, SC t_picked_v)
		{
			if (threadIdx.x == 0)
			{
				o_min_cost[blockIdx.x] = t_min_cost;
				o_picked_cost[blockIdx.x] = t_picked_cost;
				o_picked_v[blockIdx.x] = t_picked_v;
			}
		}

		template <class SC>
		__device__ __forceinline__ void getFinalCostReadTemp(SC &t_min_cost, SC &t_picked_cost, SC &t_picked_v, volatile SC *o_min_cost, volatile SC *o_picked_cost, volatile SC *o_picked_v, SC max)
		{
			if (threadIdx.x < gridDim.x)
			{
				t_min_cost = o_min_cost[threadIdx.x];
				t_picked_cost = o_picked_cost[threadIdx.x];
				t_picked_v = o_picked_v[threadIdx.x];
			}
			else
			{
				t_min_cost = max;
				t_picked_cost = max;
				t_picked_v = max;
			}
		}

		template <class SC>
		__device__ __forceinline__ void getFinalCostReadTempLarge(SC &t_min_cost, SC &t_picked_cost, SC &t_picked_v, volatile SC *o_min_cost, volatile SC *o_picked_cost, volatile SC *o_picked_v, SC max)
		{
			if (threadIdx.x < gridDim.x)
			{
				t_min_cost = o_min_cost[threadIdx.x];
				t_picked_cost = o_picked_cost[threadIdx.x];
				t_picked_v = o_picked_v[threadIdx.x];
				// read additional values
				for (int i = threadIdx.x + blockDim.x; i < gridDim.x; i += blockDim.x)
				{
					SC c_min_cost = o_min_cost[i];
					SC c_picked_cost = o_picked_cost[i];
					SC c_picked_v = o_picked_v[i];
					if (c_min_cost < t_min_cost) t_min_cost = c_min_cost;
					if (c_picked_cost < t_picked_cost) t_picked_cost = c_picked_cost;
					if (c_picked_v < t_picked_v) t_picked_v = c_picked_v;
				}
			}
			else
			{
				t_min_cost = max;
				t_picked_cost = max;
				t_picked_v = max;
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void getFinalCostWrite(MS *s, SC t_min_cost, SC t_picked_cost, SC t_picked_v)
		{
			if (threadIdx.x == 0)
			{
				s->min = t_min_cost;
				s->picked = t_picked_cost;
				s->v_jmin = t_picked_v;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getFinalCostSmall_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_picked_cost, volatile SC *o_picked_v, TC *tt, SC *v, SC max, int j_picked, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_picked_cost, t_picked_v;

			getFinalCostRead(t_min_cost, t_picked_cost, t_picked_v, j, tt, v, j_picked, max, size);
			getFinalCostCombineSmall(t_min_cost, t_picked_cost, t_picked_v);
			getFinalCostWriteTemp(o_min_cost, o_picked_cost, o_picked_v, t_min_cost, t_picked_cost, t_picked_v);

			if (semaphoreWarp(semaphore))
			{
				getFinalCostReadTemp(t_min_cost, t_picked_cost, t_picked_v, o_min_cost, o_picked_cost, o_picked_v, max);
				getFinalCostCombineSmall(t_min_cost, t_picked_cost, t_picked_v);
				getFinalCostWrite(s, t_min_cost, t_picked_cost, t_picked_v);
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getFinalCostMedium_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_picked_cost, volatile SC *o_picked_v, TC *tt, SC *v, SC max, int j_picked, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_picked_cost[8], b_picked_v[8];

			SC t_min_cost, t_picked_cost, t_picked_v;

			getFinalCostRead(t_min_cost, t_picked_cost, t_picked_v, j, tt, v, j_picked, max, size);
			getFinalCostCombineMedium(t_min_cost, t_picked_cost, t_picked_v, b_min_cost, b_picked_cost, b_picked_v);
			getFinalCostWriteTemp(o_min_cost, o_picked_cost, o_picked_v, t_min_cost, t_picked_cost, t_picked_v);

			if (semaphoreBlock(semaphore))
			{
				getFinalCostReadTemp(t_min_cost, t_picked_cost, t_picked_v, o_min_cost, o_picked_cost, o_picked_v, max);
				getFinalCostCombineMedium(t_min_cost, t_picked_cost, t_picked_v, b_min_cost, b_picked_cost, b_picked_v);
				getFinalCostWrite(s, t_min_cost, t_picked_cost, t_picked_v);
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getFinalCostLarge_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_picked_cost, volatile SC *o_picked_v, TC *tt, SC *v, SC max, int j_picked, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_picked_cost[32], b_picked_v[32];

			SC t_min_cost, t_picked_cost, t_picked_v;

			getFinalCostRead(t_min_cost, t_picked_cost, t_picked_v, j, tt, v, j_picked, max, size);
			getFinalCostCombineLarge(t_min_cost, t_picked_cost, t_picked_v, b_min_cost, b_picked_cost, b_picked_v);
			getFinalCostWriteTemp(o_min_cost, o_picked_cost, o_picked_v, t_min_cost, t_picked_cost, t_picked_v);

			if (semaphoreBlock(semaphore))
			{
				getFinalCostReadTempLarge(t_min_cost, t_picked_cost, t_picked_v, o_min_cost, o_picked_cost, o_picked_v, max);
				getFinalCostCombineLarge(t_min_cost, t_picked_cost, t_picked_v, b_min_cost, b_picked_cost, b_picked_v);
				getFinalCostWrite(s, t_min_cost, t_picked_cost, t_picked_v);
			}
		}

		template <class MS, class SC>
		__global__ void getFinalCostSmallVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_picked_cost, volatile SC *o_picked_v, SC *v, SC max, int j_picked, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_picked_cost, t_picked_v;

			getFinalCostRead(t_min_cost, t_picked_cost, t_picked_v, j, v, j_picked, max, size);
			getFinalCostCombineSmall(t_min_cost, t_picked_cost, t_picked_v);
			getFinalCostWriteTemp(o_min_cost, o_picked_cost, o_picked_v, t_min_cost, t_picked_cost, t_picked_v);

			if (semaphoreWarp(semaphore))
			{
				getFinalCostReadTemp(t_min_cost, t_picked_cost, t_picked_v, o_min_cost, o_picked_cost, o_picked_v, max);
				getFinalCostCombineSmall(t_min_cost, t_picked_cost, t_picked_v);
				getFinalCostWrite(s, t_min_cost, t_picked_cost, t_picked_v);
			}
		}

		template <class MS, class SC>
		__global__ void getFinalCostMediumVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_picked_cost, volatile SC *o_picked_v, SC *v, SC max, int j_picked, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_picked_cost[8], b_picked_v[8];

			SC t_min_cost, t_picked_cost, t_picked_v;

			getFinalCostRead(t_min_cost, t_picked_cost, t_picked_v, j, v, j_picked, max, size);
			getFinalCostCombineMedium(t_min_cost, t_picked_cost, t_picked_v, b_min_cost, b_picked_cost, b_picked_v);
			getFinalCostWriteTemp(o_min_cost, o_picked_cost, o_picked_v, t_min_cost, t_picked_cost, t_picked_v);

			if (semaphoreBlock(semaphore))
			{
				getFinalCostReadTemp(t_min_cost, t_picked_cost, t_picked_v, o_min_cost, o_picked_cost, o_picked_v, max);
				getFinalCostCombineMedium(t_min_cost, t_picked_cost, t_picked_v, b_min_cost, b_picked_cost, b_picked_v);
				getFinalCostWrite(s, t_min_cost, t_picked_cost, t_picked_v);
			}
		}

		template <class MS, class SC>
		__global__ void getFinalCostLargeVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile SC *o_picked_cost, volatile SC *o_picked_v, SC *v, SC max, int j_picked, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_picked_cost[32], b_picked_v[32];

			SC t_min_cost, t_picked_cost, t_picked_v;

			getFinalCostRead(t_min_cost, t_picked_cost, t_picked_v, j, v, j_picked, max, size);
			getFinalCostCombineLarge(t_min_cost, t_picked_cost, t_picked_v, b_min_cost, b_picked_cost, b_picked_v);
			getFinalCostWriteTemp(o_min_cost, o_picked_cost, o_picked_v, t_min_cost, t_picked_cost, t_picked_v);

			if (semaphoreBlock(semaphore))
			{
				getFinalCostReadTempLarge(t_min_cost, t_picked_cost, t_picked_v, o_min_cost, o_picked_cost, o_picked_v, max);
				getFinalCostCombineLarge(t_min_cost, t_picked_cost, t_picked_v, b_min_cost, b_picked_cost, b_picked_v);
				getFinalCostWrite(s, t_min_cost, t_picked_cost, t_picked_v);
			}
		}
	}
}
