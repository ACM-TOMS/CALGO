#pragma once

namespace lap
{
	namespace cuda
	{
		template <class SC, class TC>
		__device__ __forceinline__ void updateEstimatedVFirst(SC *min_v, TC *tt, int *picked, SC min_cost, int jmin, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j >= size) return;

			min_v[j] = (SC)tt[j] - min_cost;
			if (jmin == j) picked[j] = 1;
		}

		template <class SC, class TC>
		__device__ __forceinline__ void updateEstimatedVSecond(SC *v, SC *min_v, TC *tt, int *picked, SC min_cost, int jmin, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j >= size) return;

			SC tmp = (SC)tt[j] - min_cost;
			if (tmp < min_v[j])
			{
				v[j] = min_v[j];
				min_v[j] = tmp;
			}
			else v[j] = tmp;
			if (jmin == j) picked[j] = 1;
		}

		template <class SC, class TC>
		__device__ __forceinline__ void updateEstimatedV(SC *v, SC *min_v, TC *tt, int *picked, SC min_cost, int jmin, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j >= size) return;

			SC tmp = (SC)tt[j] - min_cost;
			if (tmp < min_v[j])
			{
				v[j] = min_v[j];
				min_v[j] = tmp;
			}
			else if (tmp < v[j]) v[j] = tmp;
			if (jmin == j) picked[j] = 1;
		}

		template <class SC>
		__device__ __forceinline__ void updateEstimatedVFirst(SC *min_v, int *picked, SC min_cost, int jmin, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j >= size) return;

			min_v[j] = -min_cost;
			if (jmin == j) picked[j] = 1;
		}

		template <class SC>
		__device__ __forceinline__ void updateEstimatedVSecond(SC *v, SC *min_v, int *picked, SC min_cost, int jmin, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j >= size) return;

			SC tmp = -min_cost;
			if (tmp < min_v[j])
			{
				v[j] = min_v[j];
				min_v[j] = tmp;
			}
			else v[j] = tmp;
			if (jmin == j) picked[j] = 1;
		}

		template <class SC>
		__device__ __forceinline__ void updateEstimatedV(SC *v, SC *min_v, int *picked, SC min_cost, int jmin, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;

			if (j >= size) return;

			SC tmp = -min_cost;
			if (tmp < min_v[j])
			{
				v[j] = min_v[j];
				min_v[j] = tmp;
			}
			else if (tmp < v[j]) v[j] = tmp;
			if (jmin == j) picked[j] = 1;
		}

		template <class MS, class SC>
		__device__ __forceinline__ void updateEstimateVGetMin(int &jmin, SC &min_cost, volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, volatile MS *s, int start, int dim2, SC max, int devices)
		{
			__shared__ int b_jmin;
			__shared__ SC b_min_cost;

			if (threadIdx.x < 32)
			{
				int sem;
				if (threadIdx.x == 0) sem = atomicInc(semaphore, gridDim.x - 1);
				sem = __shfl_sync(0xffffffff, sem, 0, 32);
				if (sem == 0)
				{
					if (threadIdx.x < devices) while (data_valid[threadIdx.x] == 0);
					__threadfence_system();
					__syncwarp(0xffffffff);

					int t_jmin;
					SC t_min_cost;
					SC t_picked_cost;
					if (threadIdx.x < devices)
					{
						t_jmin = s[threadIdx.x].jmin;
						t_picked_cost = s[threadIdx.x].picked;
						t_min_cost = s[threadIdx.x].min;
					}
					else
					{
						t_jmin = dim2;
						t_picked_cost = max;
						t_min_cost = max;
					}
					minWarpIndex(t_picked_cost, t_jmin);
					minWarp(t_min_cost);
					if (threadIdx.x == 0)
					{
						s2->min = b_min_cost = t_min_cost;
						__threadfence();
						s2->jmin = b_jmin = t_jmin;
					}
				}
				else
				{
					if (threadIdx.x == 0)
					{
						int t_jmin;
						while ((t_jmin = s2->jmin) < 0) {}
						__threadfence();
						b_jmin = t_jmin;
						b_min_cost = s2->min;
					}
				}
			}
			__syncthreads();
			jmin = b_jmin - start;
			min_cost = b_min_cost;
		}

		template <class MS, class SC>
		__device__ __forceinline__ void updateEstimateVGetMinLarge(int &jmin, SC &min_cost, volatile MS *s2, unsigned int *semaphore, volatile int *data_valid, volatile MS *s, int start, int dim2, SC max, int devices)
		{
			__shared__ int b_jmin;
			__shared__ SC b_min_cost;

			if (threadIdx.x < 32)
			{
				int sem;
				if (threadIdx.x == 0) sem = atomicInc(semaphore, gridDim.x - 1);
				sem = __shfl_sync(0xffffffff, sem, 0, 32);
				if (sem == 0)
				{
					bool is_valid = false;
					do
					{
						is_valid = true;
						for (int t = threadIdx.x; t < devices; t += 32) is_valid &= (data_valid[t] != 0);
					} while (!is_valid);

					__threadfence_system();
					__syncwarp(0xffffffff);

					int t_jmin = dim2;
					SC t_min_cost = max;
					SC t_picked_cost = max;

					for (int t = threadIdx.x; t < devices; t += 32)
					{
						int c_jmin = s[t].jmin;
						SC c_min_cost = s[t].min;
						SC c_picked_cost = s[t].picked;
						if ((c_picked_cost < t_picked_cost) || ((c_picked_cost == t_picked_cost) && (c_jmin < t_jmin)))
						{
							t_jmin = c_jmin;
							t_picked_cost = c_picked_cost;
						}
						if (c_min_cost < t_min_cost) t_min_cost = c_min_cost;
					}
					minWarpIndex(t_picked_cost, t_jmin);
					minWarp(t_min_cost);
					if (threadIdx.x == 0)
					{
						s2->min = b_min_cost = t_min_cost;
						__threadfence();
						s2->jmin = b_jmin = t_jmin;
					}
				}
				else
				{
					if (threadIdx.x == 0)
					{
						int t_jmin;
						while ((t_jmin = s2->jmin) < 0) {}
						__threadfence();
						b_jmin = t_jmin;
						b_min_cost = s2->min;
					}
				}
			}

			__syncthreads();
			jmin = b_jmin - start;
			min_cost = b_min_cost;
		}

		template <class SC, class TC>
		__global__ void updateEstimatedV_kernel(int i, SC* v, SC* min_v, TC* tt, int* picked, SC* min_cost, int* jmin, int size)
		{
			if (i == 0) updateEstimatedVFirst(min_v, tt, picked, *min_cost, *jmin, size);
			else if (i == 1) updateEstimatedVSecond(v, min_v, tt, picked, *min_cost, *jmin, size);
			else updateEstimatedV(v, min_v, tt, picked, *min_cost, *jmin, size);
		}

		template <class MS, class SC, class TC>
		__global__ void updateEstimatedVSmall_kernel(int i, SC* v, SC* min_v, volatile MS* s2, unsigned int* semaphore, TC* tt, int* picked, volatile int* data_valid, MS* s, int start, int size, int dim2, SC max, int devices)
		{
			int jmin;
			SC min_cost;
			updateEstimateVGetMin(jmin, min_cost, s2, semaphore, data_valid, s, start, dim2, max, devices);
			if (i == 0) updateEstimatedVFirst(min_v, tt, picked, min_cost, jmin, size);
			else if (i == 1) updateEstimatedVSecond(v, min_v, tt, picked, min_cost, jmin, size);
			else updateEstimatedV(v, min_v, tt, picked, min_cost, jmin, size);
		}

		template <class MS, class SC, class TC>
		__global__ void updateEstimatedVLarge_kernel(int i, SC* v, SC* min_v, volatile MS* s2, unsigned int* semaphore, TC* tt, int* picked, volatile int* data_valid, MS* s, int start, int size, int dim2, SC max, int devices)
		{
			int jmin;
			SC min_cost;
			updateEstimateVGetMinLarge(jmin, min_cost, s2, semaphore, data_valid, s, start, dim2, max, devices);
			if (i == 0) updateEstimatedVFirst(min_v, tt, picked, min_cost, jmin, size);
			else if (i == 1) updateEstimatedVSecond(v, min_v, tt, picked, min_cost, jmin, size);
			else updateEstimatedV(v, min_v, tt, picked, min_cost, jmin, size);
		}

		template <class SC>
		__global__ void updateEstimatedVVirtual_kernel(int i, SC* v, SC* min_v, int* picked, SC* min_cost, int* jmin, int size)
		{
			if (i == 0) updateEstimatedVFirst(min_v, picked, *min_cost, *jmin, size);
			else if (i == 1) updateEstimatedVSecond(v, min_v, picked, *min_cost, *jmin, size);
			else updateEstimatedV(v, min_v, picked, *min_cost, *jmin, size);
		}

		template <class MS, class SC>
		__global__ void updateEstimatedVSmallVirtual_kernel(int i, SC* v, SC* min_v, volatile MS* s2, unsigned int* semaphore, int* picked, volatile int* data_valid, MS* s, int start, int size, int dim2, SC max, int devices)
		{
			int jmin;
			SC min_cost;
			updateEstimateVGetMin(jmin, min_cost, s2, semaphore, data_valid, s, start, dim2, max, devices);
			if (i == 0) updateEstimatedVFirst(min_v, picked, min_cost, jmin, size);
			else if (i == 1) updateEstimatedVSecond(v, min_v, picked, min_cost, jmin, size);
			else updateEstimatedV(v, min_v, picked, min_cost, jmin, size);
		}

		template <class MS, class SC>
		__global__ void updateEstimatedVLargeVirtual_kernel(int i, SC* v, SC* min_v, volatile MS* s2, unsigned int* semaphore, int* picked, volatile int* data_valid, MS* s, int start, int size, int dim2, SC max, int devices)
		{
			int jmin;
			SC min_cost;
			updateEstimateVGetMinLarge(jmin, min_cost, s2, semaphore, data_valid, s, start, dim2, max, devices);
			if (i == 0) updateEstimatedVFirst(min_v, picked, min_cost, jmin, size);
			else if (i == 1) updateEstimatedVSecond(v, min_v, picked, min_cost, jmin, size);
			else updateEstimatedV(v, min_v, picked, min_cost, jmin, size);
		}
	}
}
