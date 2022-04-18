#pragma once

namespace lap
{
	namespace cuda
	{
		template <class SC, class TC>
		__device__ __forceinline__ void getMinimalCostRead(SC &t_min_cost, int &t_jmin, SC &t_min_cost_real, int j, TC *tt, SC *v, int *taken, SC max, int size, int dim2)
		{
			t_min_cost_real = max;
			t_min_cost = max;
			t_jmin = dim2;

			if (j < size)
			{
				SC t_cost = (SC)tt[j] - v[j];
				if (taken[j] == 0)
				{
					t_min_cost = t_cost;
					t_jmin = j;
				}
				t_min_cost_real = t_cost;
			}
		}

		template <class SC, class TC>
		__device__ __forceinline__ void getMinimalCostRead(SC &t_min_cost, int &t_jmin, SC &t_min_cost_real, int j, TC *tt, SC *v, int *taken, int last_taken, SC max, int size, int dim2)
		{
			t_min_cost_real = max;
			t_min_cost = max;
			t_jmin = dim2;

			if (j < size)
			{
				SC t_cost = (SC)tt[j] - v[j];
				if (j == last_taken)
				{
					taken[j] = 1;
				}
				else if (taken[j] == 0)
				{
					t_min_cost = t_cost;
					t_jmin = j;
				}
				t_min_cost_real = t_cost;
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinimalCostRead(SC &t_min_cost, int &t_jmin, SC &t_min_cost_real, int j, SC *v, int *taken, SC max, int size, int dim2)
		{
			t_min_cost_real = max;
			t_min_cost = max;
			t_jmin = dim2;

			if (j < size)
			{
				SC t_cost = -v[j];
				if (taken[j] == 0)
				{
					t_min_cost = t_cost;
					t_jmin = j;
				}
				t_min_cost_real = t_cost;
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinimalCostRead(SC &t_min_cost, int &t_jmin, SC &t_min_cost_real, int j, SC *v, int *taken, int last_taken, SC max, int size, int dim2)
		{
			t_min_cost_real = max;
			t_min_cost = max;
			t_jmin = dim2;

			if (j < size)
			{
				SC t_cost = -v[j];
				if (j == last_taken)
				{
					taken[j] = 1;
				}
				else if (taken[j] == 0)
				{
					t_min_cost = t_cost;
					t_jmin = j;
				}
				t_min_cost_real = t_cost;
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinimalCostCombineSmall(SC &t_min_cost, int &t_jmin, SC &t_min_cost_real)
		{
			minWarpIndex(t_min_cost, t_jmin);
			minWarp(t_min_cost_real);
		}

		template <class SC>
		__device__ __forceinline__ void getMinimalCostCombineTiny(SC &t_min_cost, int &t_jmin, SC &t_min_cost_real)
		{
			minWarpIndex8(t_min_cost, t_jmin);
			minWarp8(t_min_cost_real);
		}

		template <class SC>
		__device__ __forceinline__ void getMinimalCostWriteShared(SC *b_min_cost, int *b_jmin, SC *b_min_cost_real, SC t_min_cost, int t_jmin, SC t_min_cost_real)
		{
			int bidx = threadIdx.x >> 5;
			b_min_cost[bidx] = t_min_cost;
			b_min_cost_real[bidx] = t_min_cost_real;
			b_jmin[bidx] = t_jmin;
		}

		template <class SC>
		__device__ __forceinline__ void getMinimalCostReadShared(SC &t_min_cost, int &t_jmin, SC &t_min_cost_real, SC *b_min_cost, int *b_jmin, SC *b_min_cost_real)
		{
			t_min_cost = b_min_cost[threadIdx.x];
			t_min_cost_real = b_min_cost_real[threadIdx.x];
			t_jmin = b_jmin[threadIdx.x];
		}

		template <class SC>
		__device__ __forceinline__ void getMinimalCostCombineMedium(SC &t_min_cost, int &t_jmin, SC &t_min_cost_real, SC *b_min_cost, int *b_jmin, SC *b_min_cost_real)
		{
			getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
			if ((threadIdx.x & 0x1f) == 0)
			{
				getMinimalCostWriteShared(b_min_cost, b_jmin, b_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);
			}
			__syncthreads();
			if (threadIdx.x < 8)
			{
				getMinimalCostReadShared(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostCombineTiny(t_min_cost, t_jmin, t_min_cost_real);
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinimalCostCombineLarge(SC &t_min_cost, int &t_jmin, SC &t_min_cost_real, SC *b_min_cost, int *b_jmin, SC *b_min_cost_real)
		{
			getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
			if ((threadIdx.x & 0x1f) == 0)
			{
				getMinimalCostWriteShared(b_min_cost, b_jmin, b_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);
			}
			__syncthreads();
			if (threadIdx.x < 32)
			{
				getMinimalCostReadShared(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinimalCostWriteTemp(volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC t_min_cost, int t_jmin, SC t_min_cost_real)
		{
			if (threadIdx.x == 0)
			{
				o_min_cost[blockIdx.x] = t_min_cost;
				o_jmin[blockIdx.x] = t_jmin;
				o_min_cost_real[blockIdx.x] = t_min_cost_real;
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinimalCostReadTemp(SC &t_min_cost, int &t_jmin, SC &t_min_cost_real, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC max, int dim2)
		{
			if (threadIdx.x < gridDim.x)
			{
				t_min_cost = o_min_cost[threadIdx.x];
				t_jmin = o_jmin[threadIdx.x];
				t_min_cost_real = o_min_cost_real[threadIdx.x];
			}
			else
			{
				t_min_cost = max;
				t_jmin = dim2;
				t_min_cost_real = max;
			}
		}

		template <class SC>
		__device__ __forceinline__ void getMinimalCostReadTempLarge(SC &t_min_cost, int &t_jmin, SC &t_min_cost_real, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC max, int dim2)
		{
			if (threadIdx.x < gridDim.x)
			{
				t_min_cost = o_min_cost[threadIdx.x];
				t_jmin = o_jmin[threadIdx.x];
				t_min_cost_real = o_min_cost_real[threadIdx.x];
				// read additional values
				for (int i = threadIdx.x + blockDim.x; i < gridDim.x; i += blockDim.x)
				{
					SC c_min_cost = o_min_cost[i];
					int c_jmin = o_jmin[i];
					SC c_min_cost_real = o_min_cost_real[i];
					if ((c_min_cost < t_min_cost) || ((c_min_cost == t_min_cost) && (c_jmin < t_jmin)))
					{
						t_jmin = c_jmin;
						t_min_cost = c_min_cost;
					}
					if (c_min_cost_real < t_min_cost_real) t_min_cost_real = c_min_cost_real;
				}
			}
			else
			{
				t_min_cost = max;
				t_jmin = dim2;
				t_min_cost_real = max;
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void getMinimalCostWrite(MS *s, SC t_min_cost, int t_jmin, SC t_min_cost_real, SC *v, SC max, int dim2)
		{
			if (threadIdx.x == 0)
			{
				s->min = t_min_cost_real;
				s->picked = t_min_cost;
				s->jmin = t_jmin;
				if (t_jmin < dim2)
				{
					s->v_jmin = v[t_jmin];
				}
				else
				{
					s->v_jmin = max;
				}
			}
		}

		template <class MS, class SC>
		__device__ __forceinline__ void getMinimalCostWrite(MS *s, SC t_min_cost, int t_jmin, int start, SC t_min_cost_real, SC *v, SC max, int dim2)
		{
			if (threadIdx.x == 0)
			{
				s->min = t_min_cost_real;
				s->picked = t_min_cost;
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
		__global__ void getMinimalCostSmall_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, TC *tt, SC *v, int *picked, SC max, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, tt, v, picked, max, size, dim2);
			getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreWarp(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinimalCostMedium_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, TC *tt, SC *v, int *picked, SC max, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_min_cost_real[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, tt, v, picked, max, size, dim2);
			getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinimalCostLarge_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, TC *tt, SC *v, int *picked, SC max, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_min_cost_real[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, tt, v, picked, max, size, dim2);
			getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTempLarge(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinimalCostSmallSmall_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, TC *tt, SC *v, int *picked, volatile MS *s_old, SC max, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			int last_picked = getLastPickedSmall(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, tt, v, picked, last_picked, max, size, dim2);
			getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreWarp(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinimalCostSmallMedium_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, TC *tt, SC *v, int *picked, volatile MS *s_old, SC max, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_min_cost_real[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			int last_picked = getLastPicked(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, tt, v, picked, last_picked, max, size, dim2);
			getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinimalCostSmallLarge_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, TC *tt, SC *v, int *picked, volatile MS *s_old, SC max, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_min_cost_real[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			int last_picked = getLastPicked(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, tt, v, picked, last_picked, max, size, dim2);
			getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTempLarge(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinimalCostLargeSmall_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, TC *tt, SC *v, int *picked, volatile MS *s_old, SC max, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			int last_picked = getLastPickedLarge(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, tt, v, picked, last_picked, max, size, dim2);
			getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreWarp(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinimalCostLargeMedium_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, TC *tt, SC *v, int *picked, volatile MS *s_old, SC max, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_min_cost_real[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			int last_picked = getLastPickedLarge(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, tt, v, picked, last_picked, max, size, dim2);
			getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinimalCostLargeLarge_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, TC *tt, SC *v, int *picked, volatile MS *s_old, SC max, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_min_cost_real[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			int last_picked = getLastPickedLarge(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, tt, v, picked, last_picked, max, size, dim2);
			getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTempLarge(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinimalCostSingleSmall_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, TC *tt, SC *v, int *picked, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, tt, v, picked, max, size, dim2);
			getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreWarp(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0) picked[t_jmin] = 1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinimalCostSingleMedium_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, TC *tt, SC *v, int *picked, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_min_cost_real[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, tt, v, picked, max, size, dim2);
			getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0) picked[t_jmin] = 1;
			}
		}

		template <class MS, class SC, class TC>
		__global__ void getMinimalCostSingleLarge_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, TC *tt, SC *v, int *picked, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_min_cost_real[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, tt, v, picked, max, size, dim2);
			getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTempLarge(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0) picked[t_jmin] = 1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinimalCostSmallVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC *v, int *picked, SC max, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, v, picked, max, size, dim2);
			getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreWarp(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinimalCostMediumVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC *v, int *picked, SC max, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_min_cost_real[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, v, picked, max, size, dim2);
			getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinimalCostLargeVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC *v, int *picked, SC max, int start, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_min_cost_real[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, v, picked, max, size, dim2);
			getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTempLarge(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinimalCostSmallSmallVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC *v, int *picked, volatile MS *s_old, SC max, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			int last_picked = getLastPickedSmall(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, v, picked, last_picked, max, size, dim2);
			getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreWarp(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinimalCostSmallMediumVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC *v, int *picked, volatile MS *s_old, SC max, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_min_cost_real[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			int last_picked = getLastPicked(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, v, picked, last_picked, max, size, dim2);
			getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinimalCostSmallLargeVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC *v, int *picked, volatile MS *s_old, SC max, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_min_cost_real[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			int last_picked = getLastPicked(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, v, picked, last_picked, max, size, dim2);
			getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTempLarge(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinimalCostLargeSmallVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC *v, int *picked, volatile MS *s_old, SC max, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			int last_picked = getLastPickedLarge(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, v, picked, last_picked, max, size, dim2);
			getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreWarp(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinimalCostLargeMediumVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC *v, int *picked, volatile MS *s_old, SC max, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_min_cost_real[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			int last_picked = getLastPickedLarge(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, v, picked, last_picked, max, size, dim2);
			getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinimalCostLargeLargeVirtual_kernel(MS *s, MS *s2, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC *v, int *picked, volatile MS *s_old, SC max, int start, int size, int dim2, int devices)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_min_cost_real[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			int last_picked = getLastPickedLarge(s2, s_old, semaphore + 1, max, start, size, dim2, devices);

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, v, picked, last_picked, max, size, dim2);
			getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTempLarge(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, start, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0)  s2->jmin = -1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinimalCostSingleSmallVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC *v, int *picked, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, v, picked, max, size, dim2);
			getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreWarp(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineSmall(t_min_cost, t_jmin, t_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0) picked[t_jmin] = 1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinimalCostSingleMediumVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC *v, int *picked, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[8], b_min_cost_real[8];
			__shared__ int b_jmin[8];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, v, picked, max, size, dim2);
			getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTemp(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineMedium(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0) picked[t_jmin] = 1;
			}
		}

		template <class MS, class SC>
		__global__ void getMinimalCostSingleLargeVirtual_kernel(MS *s, unsigned int *semaphore, volatile SC *o_min_cost, volatile int *o_jmin, volatile SC *o_min_cost_real, SC *v, int *picked, SC max, int size, int dim2)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			__shared__ SC b_min_cost[32], b_min_cost_real[32];
			__shared__ int b_jmin[32];

			SC t_min_cost, t_min_cost_real;
			int t_jmin;

			getMinimalCostRead(t_min_cost, t_jmin, t_min_cost_real, j, v, picked, max, size, dim2);
			getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
			getMinimalCostWriteTemp(o_min_cost, o_jmin, o_min_cost_real, t_min_cost, t_jmin, t_min_cost_real);

			if (semaphoreBlock(semaphore))
			{
				getMinimalCostReadTempLarge(t_min_cost, t_jmin, t_min_cost_real, o_min_cost, o_jmin, o_min_cost_real, max, dim2);
				getMinimalCostCombineLarge(t_min_cost, t_jmin, t_min_cost_real, b_min_cost, b_jmin, b_min_cost_real);
				getMinimalCostWrite(s, t_min_cost, t_jmin, t_min_cost_real, v, max, dim2);
				if (threadIdx.x == 0) picked[t_jmin] = 1;
			}
		}
	}
}
