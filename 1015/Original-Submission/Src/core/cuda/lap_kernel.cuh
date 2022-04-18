#pragma once

// collection of all CUDA kernel required for the solver

namespace lap
{
	namespace cuda
	{
		template <class C>
		__device__ __forceinline__ void minWarp(C &value)
		{
			C value2 = __shfl_xor_sync(0xffffffff, value, 1, 32);
			if (value2 < value) value = value2;
			value2 = __shfl_xor_sync(0xffffffff, value, 2, 32);
			if (value2 < value) value = value2;
			value2 = __shfl_xor_sync(0xffffffff, value, 4, 32);
			if (value2 < value) value = value2;
			value2 = __shfl_xor_sync(0xffffffff, value, 8, 32);
			if (value2 < value) value = value2;
			value2 = __shfl_xor_sync(0xffffffff, value, 16, 32);
			if (value2 < value) value = value2;
		}

		template <class C>
		__device__ __forceinline__ void maxWarp(C &value)
		{
			C value2 = __shfl_xor_sync(0xffffffff, value, 1, 32);
			if (value2 > value) value = value2;
			value2 = __shfl_xor_sync(0xffffffff, value, 2, 32);
			if (value2 > value) value = value2;
			value2 = __shfl_xor_sync(0xffffffff, value, 4, 32);
			if (value2 > value) value = value2;
			value2 = __shfl_xor_sync(0xffffffff, value, 8, 32);
			if (value2 > value) value = value2;
			value2 = __shfl_xor_sync(0xffffffff, value, 16, 32);
			if (value2 > value) value = value2;
		}

		template <class C>
		__device__ __forceinline__ void minWarp8(C &value)
		{
			C value2 = __shfl_xor_sync(0xff, value, 1, 32);
			if (value2 < value) value = value2;
			value2 = __shfl_xor_sync(0xff, value, 2, 32);
			if (value2 < value) value = value2;
			value2 = __shfl_xor_sync(0xff, value, 4, 32);
			if (value2 < value) value = value2;
		}

		template <class C>
		__device__ __forceinline__ void maxWarp8(C &value)
		{
			C value2 = __shfl_xor_sync(0xff, value, 1, 32);
			if (value2 > value) value = value2;
			value2 = __shfl_xor_sync(0xff, value, 2, 32);
			if (value2 > value) value = value2;
			value2 = __shfl_xor_sync(0xff, value, 4, 32);
			if (value2 > value) value = value2;
		}

		template <class C>
		__device__ __forceinline__ void minWarpIndex(C &value, int &index)
		{
			C old_val = value;
			minWarp(value);
			if (old_val != value) index = 0x7fffffff;
			minWarp(index);
		}

		template <class C>
		__device__ __forceinline__ void minWarpIndex8(C &value, int &index)
		{
			C old_val = value;
			minWarp8(value);
			if (old_val != value) index = 0x7fffffff;
			minWarp8(index);
		}

		template <class C>
		__device__ __forceinline__ void minWarpIndex(C &value, int &index, int &old)
		{
			C old_val = value;
			minWarp(value);
			bool active = ((old < 0) && (old_val == value));
			int mask = __ballot_sync(0xffffffff, active);
			if (mask == 0)
			{
				active = (old_val == value);
				mask = __ballot_sync(0xffffffff, active);
			}
			int old_index = index;
			if (!active)
			{
				index = 0x7fffffff;
			}
			minWarp(index);
			int first = __ffs(__ballot_sync(0xffffffff, old_index == index)) - 1;
			old = __shfl_sync(0xffffffff, old, first, 32);
		}

		template <class C>
		__device__ __forceinline__ void minWarpIndex8(C &value, int &index, int &old)
		{
			C old_val = value;
			minWarp8(value);
			bool active = ((old < 0) && (old_val == value));
			int mask = __ballot_sync(0xff, active);
			if (mask == 0)
			{
				active = (old_val == value);
				mask = __ballot_sync(0xff, active);
			}
			int old_index = index;
			if (!active)
			{
				index = 0x7fffffff;
			}
			minWarp8(index);
			int first = __ffs(__ballot_sync(0xff, old_index == index)) - 1;
			old = __shfl_sync(0xff, old, first, 32);
		}

		__device__ __forceinline__ bool semaphoreWarp(unsigned int *semaphore)
		{
			__threadfence();
			int sem;
			if (threadIdx.x == 0) sem = atomicInc(semaphore, gridDim.x - 1);
			sem = __shfl_sync(0xffffffff, sem, 0, 32);
			return (sem == gridDim.x - 1);
		}

		__device__ __forceinline__ bool semaphoreBlock(unsigned int *semaphore)
		{
			__threadfence();
			__shared__ int sem;
			__syncthreads();
			if (threadIdx.x == 0) sem = atomicInc(semaphore, gridDim.x - 1);
			__syncthreads();
			return (sem == gridDim.x - 1);
		}

		template <class MS, class SC>
		__device__ int getLastPickedSmall(volatile MS *s2, volatile MS *s_old, unsigned int *semaphore, SC max, int start, int size, int dim2, int devices)
		{
			int last_picked;

			int sem;
			if (threadIdx.x == 0) sem = atomicInc(semaphore, gridDim.x - 1);
			sem = __shfl_sync(0xffffffff, sem, 0, 32);
			if (sem == 0)
			{
				SC t_picked_cost;
				if (threadIdx.x < devices)
				{
					while ((last_picked = s_old[threadIdx.x].jmin) < 0);
				}
				else
				{
					last_picked = dim2;
				}

				__threadfence_system();
				__syncwarp(0xffffffff);

				if (threadIdx.x < devices)
				{
					t_picked_cost = s_old[threadIdx.x].picked;
				}
				else
				{
					t_picked_cost = max;
				}

				minWarpIndex(t_picked_cost, last_picked);

				if (threadIdx.x == 0)
				{
					s2->jmin = last_picked;
				}
			}
			else
			{
				if (threadIdx.x == 0) while ((last_picked = s2->jmin) < 0) {}
				last_picked = __shfl_sync(0xffffffff, last_picked, 0, 32);
			}

			return last_picked - start;
		}

		template <class MS, class SC>
		__device__ int getLastPicked(volatile MS *s2, volatile MS *s_old, unsigned int *semaphore, SC max, int start, int size, int dim2, int devices)
		{
			int last_picked;

			__shared__ int b_last_picked;
			if (threadIdx.x < 32)
			{
				int sem;
				if (threadIdx.x == 0) sem = atomicInc(semaphore, gridDim.x - 1);
				sem = __shfl_sync(0xffffffff, sem, 0, 32);
				if (sem == 0)
				{
					SC t_picked_cost;
					if (threadIdx.x < devices)
					{
						while ((last_picked = s_old[threadIdx.x].jmin) < 0);
					}
					else
					{
						last_picked = dim2;
					}

					__threadfence_system();
					__syncwarp(0xffffffff);

					if (threadIdx.x < devices)
					{
						t_picked_cost = s_old[threadIdx.x].picked;
					}
					else
					{
						t_picked_cost = max;
					}

					minWarpIndex(t_picked_cost, last_picked);

					if (threadIdx.x == 0)
					{
						s2->jmin = b_last_picked = last_picked;
					}
				}
				else
				{
					if (threadIdx.x == 0)
					{
						while ((last_picked = s2->jmin) < 0) {}
						b_last_picked = last_picked;
					}
				}
			}
			__syncthreads();
			return b_last_picked - start;
		}

		template <class MS, class SC>
		__device__ int getLastPickedLarge(volatile MS *s2, volatile MS *s_old, unsigned int *semaphore, SC max, int start, int size, int dim2, int devices)
		{
			int last_picked = -1;

			__shared__ int b_last_picked;
			if (threadIdx.x < 32)
			{
				int sem;
				if (threadIdx.x == 0) sem = atomicInc(semaphore, gridDim.x - 1);
				sem = __shfl_sync(0xffffffff, sem, 0, 32);
				if (sem == 0)
				{
					while (last_picked < 0)
					{
						last_picked = dim2;
						for (int ii = threadIdx.x; ii < devices; ii += 32)
						{
							int c_last_picked = s_old[ii].jmin;
							if (c_last_picked < last_picked) last_picked = c_last_picked;
						}
					}

					__threadfence_system();
					__syncwarp(0xffffffff);

					last_picked = dim2;
					SC t_picked_cost = max;
					for (int ii = threadIdx.x; ii < devices; ii += 32)
					{
						int c_last_picked = s_old[ii].jmin;
						SC c_picked_cost = s_old[ii].picked;
						if ((c_picked_cost < t_picked_cost) || ((c_picked_cost == t_picked_cost) && (c_last_picked < last_picked)))
						{
							last_picked = c_last_picked;
							t_picked_cost = c_picked_cost;
						}
					}
					minWarpIndex(t_picked_cost, last_picked);

					if (threadIdx.x == 0)
					{
						s2->jmin = b_last_picked = last_picked;
					}
				}
				else
				{
					if (threadIdx.x == 0)
					{
						while ((last_picked = s2->jmin) < 0) {}
						b_last_picked = last_picked;
					}
				}
			}
			__syncthreads();
			return b_last_picked - start;
		}

		template <class MS, class SC>
		__global__ void findMaxSmall_kernel(MS *s, SC *max, SC min, int size)
		{
			int j = threadIdx.x;

			SC v_max = min;

			while (j < size)
			{
				SC c_max = max[j];
				if (c_max > v_max) v_max = c_max;
				j += blockDim.x;
			}
			maxWarp(v_max);
			if (threadIdx.x == 0) s->max = v_max;
		}

		template <class MS, class SC>
		__global__ void findMaxMedium_kernel(MS *s, SC *max, SC min, int size)
		{
			// 256 threads in 8 warps
			__shared__ SC b_max[8];

			int j = threadIdx.x;

			SC v_max = min;

			while (j < size)
			{
				SC c_max = max[j];
				if (c_max > v_max) v_max = c_max;
				j += blockDim.x;
			}
			maxWarp(v_max);
			if ((threadIdx.x & 0x1f) == 0) b_max[threadIdx.x >> 5] = v_max;
			__syncthreads();
			if (threadIdx.x >= 8) return;
			v_max = b_max[threadIdx.x];
			maxWarp8(v_max);
			if (threadIdx.x == 0) s->max = v_max;
		}

		template <class MS, class SC>
		__global__ void findMaxLarge_kernel(MS *s, SC *max, SC min, int size)
		{
			// 1024 threads in 32 warps
			__shared__ SC b_max[8];

			int j = threadIdx.x;

			SC v_max = min;

			while (j < size)
			{
				SC c_max = max[j];
				if (c_max > v_max) v_max = c_max;
				j += blockDim.x;
			}
			maxWarp(v_max);
			if ((threadIdx.x & 0x1f) == 0) b_max[threadIdx.x >> 5] = v_max;
			__syncthreads();
			if (threadIdx.x >= 32) return;
			v_max = b_max[threadIdx.x];
			maxWarp(v_max);
			if (threadIdx.x == 0) s->max = v_max;
		}

		template <class MS, class SC>
		__global__ void findMaxSmall_kernel(MS *s, SC *max, int *active, SC min, int size)
		{
			int j = threadIdx.x;

			SC v_max = min;

			while (j < size)
			{
				if (active[j] >= 0)
				{
					SC c_max = max[j];
					if (c_max > v_max) v_max = c_max;
					j += blockDim.x;
				}
			}
			maxWarp(v_max);
			if (threadIdx.x == 0) s->max = v_max;
		}

		template <class MS, class SC>
		__global__ void findMaxMedium_kernel(MS *s, SC *max, int *active, SC min, int size)
		{
			// 256 threads in 8 warps
			__shared__ SC b_max[8];

			int j = threadIdx.x;

			SC v_max = min;

			while (j < size)
			{
				if (active[j] >= 0)
				{
					SC c_max = max[j];
					if (c_max > v_max) v_max = c_max;
					j += blockDim.x;
				}
			}
			maxWarp(v_max);
			if ((threadIdx.x & 0x1f) == 0) b_max[threadIdx.x >> 5] = v_max;
			__syncthreads();
			if (threadIdx.x >= 8) return;
			v_max = b_max[threadIdx.x];
			maxWarp8(v_max);
			if (threadIdx.x == 0) s->max = v_max;
		}

		template <class MS, class SC>
		__global__ void findMaxLarge_kernel(MS *s, SC *max, int *active, SC min, int size)
		{
			// 1024 threads in 32 warps
			__shared__ SC b_max[8];

			int j = threadIdx.x;

			SC v_max = min;

			while (j < size)
			{
				if (active[j] >= 0)
				{
					SC c_max = max[j];
					if (c_max > v_max) v_max = c_max;
					j += blockDim.x;
				}
			}
			maxWarp(v_max);
			if ((threadIdx.x & 0x1f) == 0) b_max[threadIdx.x >> 5] = v_max;
			__syncthreads();
			if (threadIdx.x >= 32) return;
			v_max = b_max[threadIdx.x];
			maxWarp(v_max);
			if (threadIdx.x == 0) s->max = v_max;
		}

		template <class SC, class MS>
		__global__ void subtractMaximum_kernel(SC *v, MS *max_struct, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j >= size) return;

			v[j] -= max_struct->max;
		}

		template <class SC>
		__global__ void subtractMaximum_kernel(SC *v, SC max, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j >= size) return;

			v[j] -= max;
		}

		template <class SC, class MS>
		__global__ void subtractMaximumLimited_kernel(SC *v, MS *max_struct, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j >= size) return;

			SC tmp = v[j] - max_struct->max;
			if (tmp < SC(0)) tmp = SC(0);
			v[j] = tmp;
		}

		template <class SC>
		__global__ void subtractMaximumLimited_kernel(SC *v, SC max, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j >= size) return;

			SC tmp = v[j] - max;
			if (tmp < SC(0)) tmp = SC(0);
			v[j] = tmp;
		}

		template <class SC>
		__global__ void interpolateV_kernel(SC *v, SC *v2, double ratio2, int size)
		{
			int i = threadIdx.x + blockIdx.x * blockDim.x;

			if (i >= size) return;

			v[i] = (SC)((double)v2[i] * ratio2 + (double)v[i] * (1.0 - ratio2));
		}
	}
}

// including additional kernel (pre-conditioning)
#include "lap_update_v.cuh"
#include "lap_get_min_max_best.cuh"
#include "lap_get_min_second_best.cuh"
#include "lap_get_minimal_cost.cuh"
#include "lap_get_final_cost.cuh"
#include "lap_update_estimated_v.cuh"

// include addtional kernel (solver)
#include "lap_initialize_search.cuh"
#include "lap_continue_search.cuh"
#include "lap_update_column.cuh"

// include kernel for cost function
#include "lap_cost.cuh"
