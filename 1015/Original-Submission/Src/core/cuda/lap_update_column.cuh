#pragma once

namespace lap
{
	namespace cuda
	{
		template <class SC>
		__global__ void updateColumnPrices_kernel(char *colactive, int jmin, SC min, SC *v, SC *d, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;
			if (j >= size) return;

			if ((j == jmin) || (colactive[j] == 0))
			{
				SC dlt = min - d[j];
				v[j] -= dlt;
			}
		}

		template <class SC, class TC>
		__global__ void updateColumnPricesEpsilon_kernel(char *colactive, int jmin, SC min, SC *v, SC *d, SC *total_d, SC *total_eps, TC eps, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;
			if (j >= size) return;

			if ((j == jmin) || (colactive[j] == 0))
			{
				SC dlt = min - d[j];
				total_d[j] += dlt;
				v[j] -= dlt + eps;
				total_eps[j] += eps;
			}
		}

		template <class SC>
		__global__ void updateColumnPricesCopy_kernel(char *colactive, int jmin, SC min, SC *v, SC *d, int *dst, int *src, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;
			if (j >= size) return;

			dst[j] = src[j];
			if ((j == jmin) || (colactive[j] == 0))
			{
				SC dlt = min - d[j];
				v[j] -= dlt;
			}
		}

		template <class SC, class TC>
		__global__ void updateColumnPricesEpsilonCopy_kernel(char *colactive, int jmin, SC min, SC *v, SC *d, SC *total_d, SC *total_eps, TC eps, int *dst, int *src, int size)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;
			if (j >= size) return;

			dst[j] = src[j];
			if ((j == jmin) || (colactive[j] == 0))
			{
				SC dlt = min - d[j];
				total_d[j] += dlt;
				v[j] -= dlt + eps;
				total_eps[j] += eps;
			}
		}

		template <class SC>
		__global__ void updateColumnPricesFast_kernel(char *colactive, int jmin, SC min, SC *v, SC *d, int size, int *colsol, int csol)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;
			if (j >= size) return;
			if (j == 0) *colsol = csol;

			if ((j == jmin) || (colactive[j] == 0))
			{
				SC dlt = min - d[j];
				v[j] -= dlt;
			}
		}

		template <class SC, class TC>
		__global__ void updateColumnPricesEpsilonFast_kernel(char *colactive, int jmin, SC min, SC *v, SC *d, SC *total_d, SC *total_eps, TC eps, int size, int *colsol, int csol)
		{
			int j = threadIdx.x + blockIdx.x * blockDim.x;
			if (j >= size) return;
			if (j == 0) *colsol = csol;

			if ((j == jmin) || (colactive[j] == 0))
			{
				SC dlt = min - d[j];
				total_d[j] += dlt;
				v[j] -= dlt + eps;
				total_eps[j] += eps;
			}
		}
	}
}
