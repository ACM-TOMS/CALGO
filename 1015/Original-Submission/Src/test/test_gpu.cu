#define LAP_CUDA
// required for multiple devices
#define LAP_CUDA_OPENMP
#define LAP_QUIET
//#define LAP_DISPLAY_EVALUATED
//#define LAP_DEBUG
//#define LAP_NO_MEM_DEBUG
//#define LAP_ROWS_SCANNED
// should only be enabled for testing purposes
//#define LAP_CUDA_ALLOW_WDDM
//#define LAP_CUDA_COMPARE_CPU
#define LAP_MINIMIZE_V

//#define RANDOM_SEED 1234

#include "../lap.h"

#include <random>
#include <string>
#include "test_options.h"
#include "image.h"
#include <iomanip>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>

template <class C> void testRandom(long long min_tab, long long max_tab, long long max_memory, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent);
template <class C> void testSanity(long long min_tab, long long max_tab, long long max_memory, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent);
template <class C> void testSanityCached(long long min_cached, long long max_cached, long long max_memory, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent);
template <class C> void testGeometric(long long min_tab, long long max_tab, long long max_memory, int runs, bool epsilon, bool disjoint, std::string name_C, std::vector<int> &devs, bool silent);
template <class C> void testGeometricCached(long long min_cached, long long max_cached, long long max_memory, int runs, bool epsilon, bool disjoint, std::string name_C, std::vector<int> &devs, bool silent);
template <class C> void testRandomLowRank(long long min_tab, long long max_tab, long long max_memory, long long min_rank, long long max_rank, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent);
template <class C> void testRandomLowRankCached(long long min_cached, long long max_cached, long long max_memory, long long min_rank, long long max_rank, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent);
template <class C> void testImages(std::vector<std::string> &images, long long max_memory, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent);
template <class C> void testInteger(long long min_tab, long long max_tab, long long max_memory, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent);


int main(int argc, char* argv[])
{
	Options opt;
	int r = opt.parseOptions(argc, argv);
	if (r != 0) return r;

	if (opt.use_double)
	{
		if (opt.use_single)
		{
			if (opt.run_sanity) testSanity<double>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, false, std::string("double"), opt.devices, opt.silent);
			if (opt.run_sanity_cached) testSanityCached<double>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.runs, false, std::string("double"), opt.devices, opt.silent);
			if (opt.run_random) testRandom<double>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, false, std::string("double"), opt.devices, opt.silent);
			if (opt.run_random_low_rank) testRandomLowRank<double>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.lap_min_rank, opt.lap_max_rank, opt.runs, false, std::string("double"), opt.devices, opt.silent);
			if (opt.run_random_low_rank_cached) testRandomLowRankCached<double>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.lap_min_rank, opt.lap_max_rank, opt.runs, false, std::string("double"), opt.devices, opt.silent);
			if (opt.run_geometric) testGeometric<double>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, false, false, std::string("double"), opt.devices, opt.silent);
			if (opt.run_geometric_disjoint) testGeometric<double>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, false, true, std::string("double"), opt.devices, opt.silent);
			if (opt.run_geometric_cached) testGeometricCached<double>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.runs, false, false, std::string("double"), opt.devices, opt.silent);
			if (opt.run_geometric_disjoint_cached) testGeometricCached<double>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.runs, false, true, std::string("double"), opt.devices, opt.silent);
			if (opt.images.size() > 1) testImages<double>(opt.images, opt.lap_max_memory, opt.runs, false, std::string("double"), opt.devices, opt.silent);
		}
		if (opt.use_epsilon)
		{
			if (opt.run_sanity) testSanity<double>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, true, std::string("double"), opt.devices, opt.silent);
			if (opt.run_sanity_cached) testSanityCached<double>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.runs, true, std::string("double"), opt.devices, opt.silent);
			if (opt.run_random) testRandom<double>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, true, std::string("double"), opt.devices, opt.silent);
			if (opt.run_random_low_rank) testRandomLowRank<double>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.lap_min_rank, opt.lap_max_rank, opt.runs, true, std::string("double"), opt.devices, opt.silent);
			if (opt.run_random_low_rank_cached) testRandomLowRankCached<double>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.lap_min_rank, opt.lap_max_rank, opt.runs, true, std::string("double"), opt.devices, opt.silent);
			if (opt.run_geometric) testGeometric<double>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, true, false, std::string("double"), opt.devices, opt.silent);
			if (opt.run_geometric_disjoint) testGeometric<double>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, true, true, std::string("double"), opt.devices, opt.silent);
			if (opt.run_geometric_cached) testGeometricCached<double>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.runs, true, false, std::string("double"), opt.devices, opt.silent);
			if (opt.run_geometric_disjoint_cached) testGeometricCached<double>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.runs, true, true, std::string("double"), opt.devices, opt.silent);
			if (opt.images.size() > 1) testImages<double>(opt.images, opt.lap_max_memory, opt.runs, true, std::string("double"), opt.devices, opt.silent);
		}
	}
	if (opt.use_float)
	{
		if (opt.use_single)
		{
			if (opt.run_sanity) testSanity<float>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, false, std::string("float"), opt.devices, opt.silent);
			if (opt.run_sanity_cached) testSanityCached<float>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.runs, false, std::string("float"), opt.devices, opt.silent);
			if (opt.run_random) testRandom<float>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, false, std::string("float"), opt.devices, opt.silent);
			if (opt.run_random_low_rank) testRandomLowRank<float>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.lap_min_rank, opt.lap_max_rank, opt.runs, false, std::string("float"), opt.devices, opt.silent);
			if (opt.run_random_low_rank_cached) testRandomLowRankCached<float>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.lap_min_rank, opt.lap_max_rank, opt.runs, false, std::string("float"), opt.devices, opt.silent);
			if (opt.run_geometric) testGeometric<float>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, false, false, std::string("float"), opt.devices, opt.silent);
			if (opt.run_geometric_disjoint) testGeometric<float>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, false, true, std::string("float"), opt.devices, opt.silent);
			if (opt.run_geometric_cached) testGeometricCached<float>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.runs, false, false, std::string("float"), opt.devices, opt.silent);
			if (opt.run_geometric_disjoint_cached) testGeometricCached<float>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.runs, false, true, std::string("float"), opt.devices, opt.silent);
			if (opt.images.size() > 1) testImages<float>(opt.images, opt.lap_max_memory, opt.runs, false, std::string("float"), opt.devices, opt.silent);
		}
		if (opt.use_epsilon)
		{
			if (opt.run_sanity) testSanity<float>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, true, std::string("float"), opt.devices, opt.silent);
			if (opt.run_sanity_cached) testSanityCached<float>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.runs, true, std::string("float"), opt.devices, opt.silent);
			if (opt.run_random) testRandom<float>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, true, std::string("float"), opt.devices, opt.silent);
			if (opt.run_random_low_rank) testRandomLowRank<float>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.lap_min_rank, opt.lap_max_rank, opt.runs, true, std::string("float"), opt.devices, opt.silent);
			if (opt.run_random_low_rank_cached) testRandomLowRankCached<float>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.lap_min_rank, opt.lap_max_rank, opt.runs, true, std::string("float"), opt.devices, opt.silent);
			if (opt.run_geometric) testGeometric<float>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, true, false, std::string("float"), opt.devices, opt.silent);
			if (opt.run_geometric_disjoint) testGeometric<float>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, true, true, std::string("float"), opt.devices, opt.silent);
			if (opt.run_geometric_cached) testGeometricCached<float>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.runs, true, false, std::string("float"), opt.devices, opt.silent);
			if (opt.run_geometric_disjoint_cached) testGeometricCached<float>(opt.lap_min_cached, opt.lap_max_cached, opt.lap_max_memory, opt.runs, true, true, std::string("float"), opt.devices, opt.silent);
			if (opt.images.size() > 1) testImages<float>(opt.images, opt.lap_max_memory, opt.runs, true, std::string("float"), opt.devices, opt.silent);
		}
	}
	if (opt.run_integer)
	{
		if (opt.use_double)
		{
			if (opt.use_single) testInteger<double>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, false, std::string("double"), opt.devices, opt.silent);
			if (opt.use_epsilon) testInteger<double>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, true, std::string("double"), opt.devices, opt.silent);
		}
		if (opt.use_float)
		{
			if (opt.use_single) testInteger<float>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, false, std::string("float"), opt.devices, opt.silent);
			if (opt.use_epsilon) testInteger<float>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, true, std::string("float"), opt.devices, opt.silent);
		}
		if (opt.use_int)
		{
			if (opt.use_single) testInteger<int>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, false, std::string("int"), opt.devices, opt.silent);
			if (opt.use_epsilon) testInteger<int>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, true, std::string("int"), opt.devices, opt.silent);
		}
		if (opt.use_long)
		{
			if (opt.use_single) testInteger<long long>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, false, std::string("long long"), opt.devices, opt.silent);
			if (opt.use_epsilon) testInteger<long long>(opt.lap_min_tab, opt.lap_max_tab, opt.lap_max_memory, opt.runs, true, std::string("long long"), opt.devices, opt.silent);
		}
	}

#ifndef LAP_QUIET
	lap::allocationLogger.destroy();
#endif

	checkCudaErrors(cudaProfilerStop());
	return 0;
}

template <class SC, class TC, class CF, class STATE, class TP>
void solveCachingCUDA(TP &start_time, int N1, int N2, CF &get_cost, STATE *state, lap::cuda::Worksharing &ws, long long max_memory, int *rowsol, bool epsilon)
{
	int devices = (int)ws.device.size();

	lap::cuda::SimpleCostFunction<TC, CF, STATE> costFunction(get_cost, state, devices);

	// different cache size, so always use SLRU
	lap::cuda::CachingIterator<SC, TC, decltype(costFunction), lap::CacheSLRU> iterator(N1, N2, max_memory / sizeof(TC), costFunction, ws);

	// pre-load cache
	for (int t = 0; t < devices; t++)
	{
		int rows = std::min(N1, iterator.getCache(t).getEntries());
		for (int i = 0; i < rows; i++)
		{
			int idx;
			iterator.getCache(t).find(idx, i);
		}
		checkCudaErrors(cudaSetDevice(ws.device[t]));
		iterator.fillRows(t, rows);
	}
	for (int t = 0; t < devices; t++)
	{
		checkCudaErrors(cudaSetDevice(ws.device[t]));
		checkCudaErrors(cudaDeviceSynchronize());
	}

	lap::displayTime(start_time, "setup complete", std::cout);

	lap::cuda::solve<SC, TC>(N1, N2, costFunction, iterator, rowsol, epsilon);

	std::stringstream ss;
	ss << "cost = " << std::setprecision(std::numeric_limits<SC>::max_digits10) << lap::cuda::cost<SC, TC>(N1, N2, costFunction, rowsol, ws.stream[0]);
	lap::displayTime(start_time, ss.str().c_str(), std::cout);
}

template <class SC, class TC, class CF, class STATE, class TP>
void solveDirectCUDA(TP& start_time, int N1, int N2, CF& get_cost, STATE* state, lap::cuda::Worksharing& ws, long long max_memory, int* rowsol, bool epsilon)
{
	int devices = (int)ws.device.size();

	lap::cuda::SimpleCostFunction<TC, CF, STATE> costFunction(get_cost, state, devices);
	lap::cuda::GPUTableCost<TC> costMatrix(N1, N2, costFunction, ws);

	lap::cuda::DirectIterator<SC, TC, decltype(costMatrix)> iterator(N1, N2, costMatrix, ws);

	lap::displayTime(start_time, "setup complete", std::cout);

	lap::cuda::solve<SC, TC>(N1, N2, costMatrix, iterator, rowsol, epsilon);

	std::stringstream ss;
	ss << "cost = " << std::setprecision(std::numeric_limits<SC>::max_digits10) << lap::cuda::cost<SC, TC>(N1, N2, costFunction, rowsol, ws.stream[0]);
	lap::displayTime(start_time, ss.str().c_str(), std::cout);
}

template <class SC, class TC, class CF, class TP>
void solveCachingTableCUDA(TP &start_time, int N1, int N2, CF &get_cost_cpu, lap::cuda::Worksharing &ws, long long max_memory, int *rowsol, bool epsilon, bool sequential, bool pinned)
{
	lap::cuda::CpuCostFunction<TC, decltype(get_cost_cpu)> cpuCostFunction(get_cost_cpu, sequential);
	lap::cuda::CPUTableCost<TC> costMatrix(N1, N2, cpuCostFunction, ws, pinned);

	int devices = (int)ws.device.size();

	// different cache size, so always use SLRU
	lap::cuda::CachingIterator<SC, TC, decltype(costMatrix), lap::CacheSLRU> iterator(N1, N2, max_memory / sizeof(TC), costMatrix, ws);

	// pre-load cache
	for (int t = 0; t < devices; t++)
	{
		int rows = std::min(N1, iterator.getCache(t).getEntries());
		for (int i = 0; i < rows; i++)
		{
			int idx;
			iterator.getCache(t).find(idx, i);
		}
		checkCudaErrors(cudaSetDevice(ws.device[t]));
		iterator.fillRows(t, rows);
	}
	for (int t = 0; t < devices; t++)
	{
		checkCudaErrors(cudaSetDevice(ws.device[t]));
		checkCudaErrors(cudaDeviceSynchronize());
	}

	lap::displayTime(start_time, "setup complete", std::cout);

	lap::cuda::solve<SC, TC>(N1, N2, costMatrix, iterator, rowsol, epsilon);

	std::stringstream ss;
	ss << "cost = " << std::setprecision(std::numeric_limits<SC>::max_digits10) << lap::cost<SC>(N1, N2, costMatrix, rowsol);
	lap::displayTime(start_time, ss.str().c_str(), std::cout);
}

template <class SC, class TC, class CF, class TP>
void solveDirectTableCUDA(TP& start_time, int N1, int N2, CF& get_cost_cpu, lap::cuda::Worksharing& ws, long long max_memory, int* rowsol, bool epsilon, bool sequential, bool pinned)
{
	lap::cuda::CpuCostFunction<TC, decltype(get_cost_cpu)> cpuCostFunction(get_cost_cpu, sequential);
	lap::cuda::CPUTableCost<TC> hostMatrix(N1, N2, cpuCostFunction, ws, pinned);
	lap::cuda::GPUTableCost<TC> costMatrix(N1, N2, hostMatrix, ws);

	int devices = (int)ws.device.size();

	lap::cuda::DirectIterator<SC, TC, decltype(costMatrix)> iterator(N1, N2, costMatrix, ws);

	lap::displayTime(start_time, "setup complete", std::cout);

	lap::cuda::solve<SC, TC>(N1, N2, costMatrix, iterator, rowsol, epsilon);

	std::stringstream ss;
	ss << "cost = " << std::setprecision(std::numeric_limits<SC>::max_digits10) << lap::cost<SC>(N1, N2, hostMatrix, rowsol);
	lap::displayTime(start_time, ss.str().c_str(), std::cout);
}

template <class SC, class TC, class CF, class STATE, class TP>
void solveCUDA(TP& start_time, int N1, int N2, CF& get_cost_cpu, STATE* state, lap::cuda::Worksharing& ws, long long max_memory, int* rowsol, bool epsilon, bool silent)
{
	bool useTable = true;
	int devices = (int)ws.device.size();
	for (int t = 0; t < devices; t++)
	{
		long long required = (long long)N1 * (long long)(ws.part[t].second - ws.part[t].first) * sizeof(TC);
		if (required > max_memory) useTable = false;
	}
	if (useTable)
	{
		if (!silent) lap::displayTime(start_time, "Solver using GPU table", std::cout);
		solveDirectCUDA<SC, TC, CF, STATE, TP>(start_time, N1, N2, get_cost_cpu, state, ws, max_memory, rowsol, epsilon);
	}
	else
	{
		if (!silent) lap::displayTime(start_time, "Solver using GPU caching", std::cout);
		solveCachingCUDA<SC, TC, CF, STATE, TP>(start_time, N1, N2, get_cost_cpu, state, ws, max_memory, rowsol, epsilon);
	}
}

template <class SC, class TC, class CF, class TP>
void solveTableCUDA(TP& start_time, int N1, int N2, CF& get_cost_cpu, lap::cuda::Worksharing& ws, long long max_memory, int* rowsol, bool epsilon, bool sequential, bool pinned, bool silent)
{
	bool useTable = true;
	int devices = (int)ws.device.size();
	for (int t = 0; t < devices; t++)
	{
		long long required = (long long)N1 * (long long)(ws.part[t].second - ws.part[t].first) * sizeof(TC);
		if (required > max_memory) useTable = false;
	}
	if (useTable)
	{
		if (!silent) lap::displayTime(start_time, "Solver using GPU table", std::cout);
		solveDirectTableCUDA<SC, TC, CF, TP>(start_time, N1, N2, get_cost_cpu, ws, max_memory, rowsol, epsilon, sequential, pinned);
	}
	else
	{
		if (!silent)
		{
			if (pinned) lap::displayTime(start_time, "Solver using GPU caching of pinned CPU table", std::cout);
			else lap::displayTime(start_time, "Solver using GPU caching of pageable CPU table", std::cout);
		}
		solveCachingTableCUDA<SC, TC, CF, TP>(start_time, N1, N2, get_cost_cpu, ws, max_memory, rowsol, epsilon, sequential, pinned);
	}
}

// needs to be declared outside of a function
template <class C>
struct GeometricState
{
	C *tab_s;
	C *tab_t;
};

template <class C>
void testGeometricCached(long long min_cached, long long max_cached, long long max_memory, int runs, bool epsilon, bool disjoint, std::string name_C, std::vector<int> &devs, bool silent)
{
	for (long long NN = min_cached * min_cached; NN <= max_cached * max_cached; NN <<= 1)
	{
		for (int r = 0; r < runs; r++)
		{
			int N = (int)floor(sqrt((double)NN));

			int max_devices = 1;
			while ((size_t)max_devices * max_memory / sizeof(C) < (size_t)N * (size_t)N) max_devices++;

			std::cout << "Geometric";
			if (disjoint) std::cout << " Disjoint";
			std::cout << " R^2<" << name_C << "> " << N << "x" << N << " (" << (double)max_memory / 1073741824.0 << "GB / GPU)";
			if (epsilon) std::cout << " with epsilon scaling";
			std::cout << std::endl;

			auto start_time = std::chrono::high_resolution_clock::now();

			std::uniform_real_distribution<C> distribution(0.0, 1.0);
#ifdef RANDOM_SEED
			std::mt19937_64 generator(RANDOM_SEED);
#else
			std::random_device rd;
			std::mt19937_64 generator(rd());
#endif

			C *tab_s = new C[2 * N];
			C *tab_t = new C[2 * N];
			for (int i = 0; i < N; i++)
			{
				tab_s[i] = distribution(generator);
				tab_t[i] = distribution(generator);
				tab_s[i + N] = distribution(generator);
				tab_t[i + N] = distribution(generator);
			}

			// order of coordinates is different, first all x then all y
			if (disjoint)
			{
				for (int i = 0; i < N; i++)
				{
					if ((i << 1) < N)
					{
						tab_t[i] += C(1.0);
					}
					else
					{
						tab_s[i] += C(1.0);
						tab_s[i + N] += C(1.0);
						tab_t[i + N] += C(1.0);
					}
				}
			}

			lap::cuda::Worksharing ws(N, 256, devs, max_devices, silent);
			int num_enabled = (int)ws.device.size();

			typedef GeometricState<C> State;

			State *d_state = new State[num_enabled];

			for (int i = 0; i < num_enabled; i++)
			{
				d_state[i].tab_s = 0;
				d_state[i].tab_t = 0;
			}

			for (int i = 0; i < num_enabled; i++)
			{
				checkCudaErrors(cudaSetDevice(ws.device[i]));
				lapAllocDevice(d_state[i].tab_s, 2 * N, __FILE__, __LINE__);
				lapAllocDevice(d_state[i].tab_t, 2 * N, __FILE__, __LINE__);
				checkCudaErrors(cudaMemcpy(d_state[i].tab_s, tab_s, 2 * N * sizeof(C), cudaMemcpyHostToDevice));
				checkCudaErrors(cudaMemcpy(d_state[i].tab_t, tab_t, 2 * N * sizeof(C), cudaMemcpyHostToDevice));
			}

			int *rowsol = new int[N];

			// cost function
			auto get_cost = [N] __device__(int x, int y, State &state)
			{
				float d0 = state.tab_s[x] - state.tab_t[y];
				float d1 = state.tab_s[x + N] - state.tab_t[y + N];
				return d0 * d0 + d1 * d1;
			};

			solveCUDA<C, C>(start_time, N, N, get_cost, d_state, ws, max_memory, rowsol, epsilon, silent);

			for (int i = 0; i < num_enabled; i++)
			{
				checkCudaErrors(cudaSetDevice(ws.device[i]));
				lapFreeDevice(d_state[i].tab_s);
				lapFreeDevice(d_state[i].tab_t);
			}

			delete[] rowsol;
			delete[] tab_s;
			delete[] tab_t;
			delete[] d_state;
		}
	}
}

template <class C, class GETCOST, class STATE>
__global__
void getGTCost_sanity_kernel(C *cost, GETCOST getcost, STATE state, int N)
{
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	if (x >= N) return;

	cost[x] = getcost(x, x, state);
}

// needs to be declared outside of a function
template <class C>
struct SanityState
{
	C *vec;
};

template <class C>
void testSanityCached(long long min_cached, long long max_cached, long long max_memory, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent)
{
	for (long long NN = min_cached * min_cached; NN <= max_cached * max_cached; NN <<= 1)
	{
		for (int r = 0; r < runs; r++)
		{
			int N = (int)floor(sqrt((double)NN));

			int max_devices = 1;
			while ((size_t)max_devices * max_memory / sizeof(C) < (size_t)N * (size_t)N) max_devices++;

			std::cout << "Sanity<" << name_C << "> " << N << "x" << N << " (" << (double)max_memory / 1073741824.0 << "GB / GPU)";
			if (epsilon) std::cout << " with epsilon scaling";
			std::cout << std::endl;

			auto start_time = std::chrono::high_resolution_clock::now();

			std::uniform_real_distribution<C> distribution(0.0, 1.0);
#ifdef RANDOM_SEED
			std::mt19937_64 generator(RANDOM_SEED);
#else
			std::random_device rd;
			std::mt19937_64 generator(rd());
#endif

			C *vec = new C[N << 1];

			for (long long i = 0; i < N << 1; i++) vec[i] = distribution(generator);

			lap::cuda::Worksharing ws(N, 256, devs, max_devices, silent);
			int num_enabled = (int)ws.device.size();

			typedef SanityState<C> State;
			State *d_state = new State[num_enabled];

			for (int i = 0; i < num_enabled; i++)
			{
				d_state[i].vec = 0;
			}

			for (int i = 0; i < num_enabled; i++)
			{
				checkCudaErrors(cudaSetDevice(ws.device[i]));
				lapAllocDevice(d_state[i].vec, 2 * N, __FILE__, __LINE__);
				checkCudaErrors(cudaMemcpy(d_state[i].vec, vec, 2 * N * sizeof(C), cudaMemcpyHostToDevice));
			}

			int *rowsol = new int[N];

			// cost function
			auto get_cost = [N] __device__(int x, int y, State &state)
			{
				C r = state.vec[x] + state.vec[y + N];
				if (x != y) r += C(0.1);

				return r;
			};

			solveCUDA<C, C>(start_time, N, N, get_cost, d_state, ws, max_memory, rowsol, epsilon, silent);

			bool passed = true;
			for (long long i = 0; (passed) && (i < N); i++)
			{
				passed &= (rowsol[i] == i);
			}
			std::stringstream ss;
			if (passed) ss << "test passed: ";
			else ss << "test failed: ";
			{
				// set device back to 0
				checkCudaErrors(cudaSetDevice(ws.device[0]));
				C my_cost(0);
				C *row = new C[N];
				// calculate costs directly
				{
					C *d_row;
					lapAllocDevice(d_row, N, __FILE__, __LINE__);
					dim3 block_size, grid_size;
					block_size.x = 256;
					grid_size.x = (N + block_size.x - 1) / block_size.x;
					getGTCost_sanity_kernel<<<grid_size, block_size>>>(d_row, get_cost, d_state[0], N);
					checkCudaErrors(cudaMemcpy(row, d_row, N * sizeof(C), cudaMemcpyDeviceToHost));
					lapFreeDevice(d_row);
				}
				for (int i = 0; i < N; i++) my_cost += row[i];
				delete[] row;
				ss << "ground truth cost = " << std::setprecision(std::numeric_limits<C>::max_digits10) << my_cost;
			}
			lap::displayTime(start_time, ss.str().c_str(), std::cout);

			for (int i = 0; i < num_enabled; i++)
			{
				checkCudaErrors(cudaSetDevice(ws.device[i]));
				lapFreeDevice(d_state[i].vec);
			}

			delete[] rowsol;
			delete[] vec;
			delete[] d_state;
		}
	}
}

template <class C>
void testRandomLowRankCached(long long min_cached, long long max_cached, long long max_memory, long long min_rank, long long max_rank, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent)
{
	for (long long rank = min_rank; rank <= max_rank; rank <<= 1)
	{
		for (long long NN = min_cached * min_cached; NN <= max_cached * max_cached; NN <<= 1)
		{
			for (int r = 0; r < runs; r++)
			{
				int N = (int)floor(sqrt((double)NN));

				int max_devices = 1;
				while ((size_t)max_devices * max_memory / sizeof(C) < (size_t)N * (size_t)N) max_devices++;

				std::cout << "RandomLowRank<" << name_C << "> " << N << "x" << N << " (" << (double)max_memory / 1073741824.0 << "GB / GPU)";
				if (epsilon) std::cout << " with epsilon scaling";
				std::cout << std::endl;

				auto start_time = std::chrono::high_resolution_clock::now();

				std::uniform_real_distribution<C> distribution(0.0, 1.0);
#ifdef RANDOM_SEED
				std::mt19937_64 generator(RANDOM_SEED);
#else
				std::random_device rd;
				std::mt19937_64 generator(rd());
#endif

				// The following matrix will have at most the seletcted rank.
				C *vec = new C[N * rank];
				for (long long i = 0; i < rank; i++)
				{
					for (long long j = 0; j < N; j++) vec[i * N + j] = distribution(generator);
				}

				lap::cuda::Worksharing ws(N, 256, devs, max_devices, silent);
				int num_enabled = (int)ws.device.size();

				typedef SanityState<C> State;
				State *d_state = new State[num_enabled];

				for (int i = 0; i < num_enabled; i++)
				{
					d_state[i].vec = 0;
				}

				for (int i = 0; i < num_enabled; i++)
				{
					checkCudaErrors(cudaSetDevice(ws.device[i]));
					lapAllocDevice(d_state[i].vec, N * rank, __FILE__, __LINE__);
					checkCudaErrors(cudaMemcpy(d_state[i].vec, vec, N * rank * sizeof(C), cudaMemcpyHostToDevice));
				}

				int *rowsol = new int[N];

				// cost function
				auto get_cost = [rank, N] __device__(int x, int y, State &state)
				{
					C sum(0);
#pragma unroll(8)
					for (long long k = 0; k < rank; k++)
					{
						sum += state.vec[k * N + x] * state.vec[k * N + y];
					}
					sum /= C(rank);

					return sum;
				};

				solveCUDA<C, C>(start_time, N, N, get_cost, d_state, ws, max_memory, rowsol, epsilon, silent);

				for (int i = 0; i < num_enabled; i++)
				{
					checkCudaErrors(cudaSetDevice(ws.device[i]));
					lapFreeDevice(d_state[i].vec);
				}

				delete[] rowsol;
				delete[] vec;
				delete[] d_state;
			}
		}
	}
}

template <class C>
void testInteger(long long min_tab, long long max_tab, long long max_memory, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent)
{
	// random costs (directly supply cost matrix)
	//for (int range = 0; range < 3; range++)
	int range = 2;
	{
		for (long long NN = min_tab * min_tab; NN <= max_tab * max_tab; NN <<= 1)
		{
			for (int r = 0; r < runs; r++)
			{
				int N = (int)floor(sqrt((double)NN));

				int max_devices = 1;
				while ((size_t)max_devices * max_memory / sizeof(C) < (size_t)N * (size_t)N) max_devices++;

				std::cout << "Integer";
				std::cout << "<" << name_C << " ";
				if (range == 0) std::cout << "1/10n";
				else if (range == 1) std::cout << "n";
				else std::cout << "10n";
				std::cout << "> " << N << "x" << N << " table";
				if (epsilon) std::cout << " with epsilon scaling";
				std::cout << std::endl;

				int n;
				if (range == 0) n = N / 10;
				else if (range == 1) n = N;
				else n = 10 * N;
				std::uniform_int_distribution<int> distribution(0, n);
#ifdef RANDOM_SEED
				std::mt19937_64 generator(RANDOM_SEED);
#else
				std::random_device rd;
				std::mt19937_64 generator(rd());
#endif

				auto start_time = std::chrono::high_resolution_clock::now();

				auto get_cost = [&distribution, &generator](int x, int y) -> int
				{
					return distribution(generator);
				};

				lap::cuda::Worksharing ws(N, 256, devs, max_devices, silent);

				int *rowsol = new int[N];

				solveTableCUDA<C, int>(start_time, N, N, get_cost, ws, max_memory, rowsol, epsilon, true, N < max_tab, silent);

				delete[] rowsol;
			}
		}
	}
}

template <class C> void testRandom(long long min_tab, long long max_tab, long long max_memory, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent)
{
	// random costs (directly supply cost matrix)
	for (long long NN = min_tab * min_tab; NN <= max_tab * max_tab; NN <<= 1)
	{
		for (int r = 0; r < runs; r++)
		{
			int N = (int)floor(sqrt((double)NN));

			int max_devices = 1;
			while ((size_t)max_devices * max_memory / sizeof(C) < (size_t)N * (size_t)N) max_devices++;

			std::cout << "Random";
			std::cout << "<" << name_C << "> " << N << "x" << N << " table";
			if (epsilon) std::cout << " with epsilon scaling";
			std::cout << std::endl;

			std::uniform_real_distribution<C> distribution(0.0, 1.0);
#ifdef RANDOM_SEED
			std::mt19937_64 generator(RANDOM_SEED);
#else
			std::random_device rd;
			std::mt19937_64 generator(rd());
#endif

			auto start_time = std::chrono::high_resolution_clock::now();

			int *rowsol = new int[N];

			auto get_cost = [&distribution, &generator](int x, int y) -> C
			{
				return distribution(generator);
			};

			lap::cuda::Worksharing ws(N, 256, devs, max_devices, silent);

			solveTableCUDA<C, C>(start_time, N, N, get_cost, ws, max_memory, rowsol, epsilon, true, N < max_tab, silent);

			delete[] rowsol;
		}
	}
}

template <class C> void testSanity(long long min_tab, long long max_tab, long long max_memory, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent)
{
	// random costs (directly supply cost matrix)
	for (long long NN = min_tab * min_tab; NN <= max_tab * max_tab; NN <<= 1)
	{
		for (int r = 0; r < runs; r++)
		{
			int N = (int)floor(sqrt((double)NN));

			int max_devices = 1;
			while ((size_t)max_devices * max_memory / sizeof(C) < (size_t)N * (size_t)N) max_devices++;

			std::cout << "Sanity";
			std::cout << "<" << name_C << "> " << N << "x" << N << " table";
			if (epsilon) std::cout << " with epsilon scaling";
			std::cout << std::endl;

			std::uniform_real_distribution<C> distribution(0.0, 1.0);
#ifdef RANDOM_SEED
			std::mt19937_64 generator(RANDOM_SEED);
#else
			std::random_device rd;
			std::mt19937_64 generator(rd());
#endif

			auto start_time = std::chrono::high_resolution_clock::now();

			int *rowsol = new int[N];

			C *vec = new C[N << 1];

			for (long long i = 0; i < N << 1; i++) vec[i] = distribution(generator);

			// cost functions
			auto get_cost = [&vec, &N](int x, int y) -> C
			{
				C r = vec[x] + vec[y + N];
				if (x == y) return r;
				else return r + C(0.1);
			};


			lap::cuda::Worksharing ws(N, 256, devs, max_devices, silent);

			solveTableCUDA<C, C>(start_time, N, N, get_cost, ws, max_memory, rowsol, epsilon, false, N < max_tab, silent);

			bool passed = true;
			for (long long i = 0; (passed) && (i < N); i++)
			{
				passed &= (rowsol[i] == i);
			}
			std::stringstream ss;
			if (passed) ss << "test passed: ";
			else ss << "test failed: ";
			C real_cost(0);
			for (int i = 0; i < N; i++) real_cost += get_cost(i, i);
			ss << "ground truth cost = " << std::setprecision(std::numeric_limits<C>::max_digits10) << real_cost;
			lap::displayTime(start_time, ss.str().c_str(), std::cout);

			delete[] vec;
			delete[] rowsol;
		}
	}
}

template <class C> void testGeometric(long long min_tab, long long max_tab, long long max_memory, int runs, bool epsilon, bool disjoint, std::string name_C, std::vector<int> &devs, bool silent)
{
	// geometric costs in R^2 (supply function for calculating cost matrix)
	for (long long NN = min_tab * min_tab; NN <= max_tab * max_tab; NN <<= 1)
	{
		for (int r = 0; r < runs; r++)
		{
			int N = (int)floor(sqrt((double)NN));

			int max_devices = 1;
			while ((size_t)max_devices * max_memory / sizeof(C) < (size_t)N * (size_t)N) max_devices++;

			std::cout << "Geometric";
			if (disjoint) std::cout << " Disjoint";
			std::cout << " R^2<" << name_C << "> " << N << "x" << N << " table";
			if (epsilon) std::cout << " with epsilon scaling";
			std::cout << std::endl;

			auto start_time = std::chrono::high_resolution_clock::now();

			std::uniform_real_distribution<C> distribution(0.0, 1.0);
#ifdef RANDOM_SEED
			std::mt19937_64 generator(RANDOM_SEED);
#else
			std::random_device rd;
			std::mt19937_64 generator(rd());
#endif

			C *tab_s = new C[2 * N];
			C *tab_t = new C[2 * N];

			for (int i = 0; i < 2 * N; i++)
			{
				tab_s[i] = distribution(generator);
				tab_t[i] = distribution(generator);
			}

			if (disjoint)
			{
				for (int i = 0; i < 2 * N; i += 2)
				{
					if (i < N)
					{
						tab_t[i] += C(1);
					}
					else
					{
						tab_s[i] += C(1);
						tab_s[i + 1] += C(1);
						tab_t[i + 1] += C(1);
					}
				}
			}

			// cost function
			auto get_cost = [&tab_s, &tab_t](int x, int y) -> C
			{
				int xx = x + x;
				int yy = y + y;
				C d0 = tab_s[xx] - tab_t[yy];
				C d1 = tab_s[xx + 1] - tab_t[yy + 1];
				return d0 * d0 + d1 * d1;
			};

			int *rowsol = new int[N];

			lap::cuda::Worksharing ws(N, 256, devs, max_devices, silent);

			solveTableCUDA<C, C>(start_time, N, N, get_cost, ws, max_memory, rowsol, epsilon, false, N < max_tab, silent);

			delete[] tab_s;
			delete[] tab_t;
			delete[] rowsol;
		}
	}
}

template <class C> void testRandomLowRank(long long min_tab, long long max_tab, long long max_memory, long long min_rank, long long max_rank, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent)
{
	// random costs (directly supply cost matrix)
	for (long long rank = min_rank; rank <= max_rank; rank <<= 1)
	{
		for (long long NN = min_tab * min_tab; NN <= max_tab * max_tab; NN <<= 1)
		{
			for (int r = 0; r < runs; r++)
			{
				int N = (int)floor(sqrt((double)NN));

				int max_devices = 1;
				while ((size_t)max_devices * max_memory / sizeof(C) < (size_t)N * (size_t)N) max_devices++;

				std::cout << "RandomLowRank<" << name_C << "> " << N << "x" << N << " table rank = " << rank;
				if (epsilon) std::cout << " with epsilon scaling";
				std::cout << std::endl;

				auto start_time = std::chrono::high_resolution_clock::now();

				std::uniform_real_distribution<C> distribution(0.0, 1.0);
#ifdef RANDOM_SEED
				std::mt19937_64 generator(RANDOM_SEED);
#else
				std::random_device rd;
				std::mt19937_64 generator(rd());
#endif

				// The following matrix will have at most the seletcted rank.
				C *vec = new C[N * rank];
				C max_vec;
				C min_vec;
				for (long long i = 0; i < rank; i++)
				{
					for (long long j = 0; j < N; j++) vec[i * N + j] = distribution(generator);
					max_vec = vec[i * N];
					for (long long j = 1; j < N; j++) max_vec = std::max(max_vec, vec[i * N + j]);
					min_vec = vec[i * N];
					for (long long j = 1; j < N; j++) min_vec = std::min(min_vec, vec[i * N + j]);
				}

				// cost function
				auto get_cost = [&vec, &N, &rank, &max_vec](int x, int y) -> C
				{
					C sum(0);
					for (long long k = 0; k < rank; k++)
					{
						sum += vec[k * N + x] * vec[k * N + y];
					}
					return sum / C(rank);
				};

				int *rowsol = new int[N];

				lap::cuda::Worksharing ws(N, 256, devs, max_devices, silent);

				solveTableCUDA<C, C>(start_time, N, N, get_cost, ws, max_memory, rowsol, epsilon, false, N < max_tab, silent);

				delete[] vec;
				delete[] rowsol;
			}
		}
	}
}

// needs to be declared outside of a function
struct ImagesState
{
	unsigned char *c00;
	unsigned char *c01;
	unsigned char *c02;
	unsigned char *c10;
	unsigned char *c11;
	unsigned char *c12;
};

template <class C> void testImages(std::vector<std::string> &images, long long max_memory, int runs, bool epsilon, std::string name_C, std::vector<int> &devs, bool silent)
{
	std::cout << "Comparing images ";
	if (epsilon) std::cout << " with epsilon scaling";
	std::cout << std::endl;
	for (unsigned int a = 0; a < images.size() - 1; a++)
	{
		for (unsigned int b = a + 1; b < images.size(); b++)
		{
			PPMImage img_a(images[a]);
			PPMImage img_b(images[b]);
			std::cout << "Comparing image \"" << images[a] << "\" (" << img_a.width << "x" << img_a.height << ") with image \"" << images[b] << "\" (" << img_b.width << "x" << img_b.height << ")." << std::endl;
			for (int r = 0; r < runs; r++)
			{
				auto start_time = std::chrono::high_resolution_clock::now();

				int N1 = std::min(img_a.width * img_a.height, img_b.width * img_b.height);
				int N2 = std::max(img_a.width * img_a.height, img_b.width * img_b.height);

				long long max_memory_local = max_memory - (N1 * 3 + N2 + 3);

				int max_devices = 1;
				while ((size_t)max_devices * max_memory / sizeof(C) < (size_t)N1 * (size_t)N2) max_devices++;

				lap::cuda::Worksharing ws(N2, 256, devs, max_devices, silent);
				int num_devices = (int)ws.device.size();
				typedef ImagesState State;
				// make sure img[0] is at most as large as img[1]
				PPMImage *img[2];
				img[0] = new PPMImage[num_devices];
				img[1] = new PPMImage[num_devices];
				// rearrange data for GPU
				int size_a = img_a.width * img_a.height;
				unsigned char *buf_a = new unsigned char[size_a * 3];
				int size_b = img_b.width * img_b.height;
				unsigned char *buf_b = new unsigned char[size_b * 3];
				for (int y = 0; y < img_a.height; y++)
				{
					for (int x = 0; x < img_a.width; x++)
					{
						int off = x + y * img_a.width;
						buf_a[off] = img_a.raw[off * 3];
						buf_a[off + size_a] = img_a.raw[off * 3 + 1];
						buf_a[off + 2 * size_a] = img_a.raw[off * 3 + 2];
					}
				}
				for (int y = 0; y < img_b.height; y++)
				{
					for (int x = 0; x < img_b.width; x++)
					{
						int off = x + y * img_b.width;
						buf_b[off] = img_b.raw[off * 3];
						buf_b[off + size_b] = img_b.raw[off * 3 + 1];
						buf_b[off + 2 * size_b] = img_b.raw[off * 3 + 2];
					}
				}
				for (int t = 0; t < num_devices; t++)
				{
					checkCudaErrors(cudaSetDevice(ws.device[t]));
					if (img_a.width * img_a.height < img_b.width * img_b.height)
					{
						img[0][t].width = img_a.width;
						img[0][t].height = img_a.height;
						img[0][t].max_val = img_a.max_val;
						lapAllocDevice(img[0][t].raw, img[0][t].width * img[0][t].height * 3, __FILE__, __LINE__);
						img[1][t].width = img_b.width;
						img[1][t].height = img_b.height;
						img[1][t].max_val = img_b.max_val;
						lapAllocDevice(img[1][t].raw, img[1][t].width * img[1][t].height * 3, __FILE__, __LINE__);

						checkCudaErrors(cudaMemcpyAsync(img[0][t].raw, buf_a, img[0][t].width * img[0][t].height * 3, cudaMemcpyHostToDevice));
						checkCudaErrors(cudaMemcpyAsync(img[1][t].raw, buf_b, img[1][t].width * img[1][t].height * 3, cudaMemcpyHostToDevice));
					}
					else
					{
						img[0][t].width = img_b.width;
						img[0][t].height = img_b.height;
						img[0][t].max_val = img_b.max_val;
						lapAllocDevice(img[0][t].raw, img[0][t].width * img[0][t].height * 3, __FILE__, __LINE__);
						img[1][t].width = img_a.width;
						img[1][t].height = img_a.height;
						img[1][t].max_val = img_a.max_val;
						lapAllocDevice(img[1][t].raw, img[1][t].width * img[1][t].height * 3, __FILE__, __LINE__);

						checkCudaErrors(cudaMemcpyAsync(img[1][t].raw, buf_a, img[1][t].width * img[1][t].height * 3, cudaMemcpyHostToDevice));
						checkCudaErrors(cudaMemcpyAsync(img[0][t].raw, buf_b, img[0][t].width * img[0][t].height * 3, cudaMemcpyHostToDevice));
					}
					checkCudaErrors(cudaDeviceSynchronize());
				}
				delete[] buf_a;
				delete[] buf_b;

				// setup state and other arguments
				typedef ImagesState State;
				State *d_state = new State[num_devices];
				int w0 = img[0][0].width;
				int h0 = img[0][0].height;
				int mval0 = img[0][0].max_val;
				int w1 = img[1][0].width;
				int h1 = img[1][0].height;
				int mval1 = img[1][0].max_val;
				int size_0 = img[0][0].width * img[0][0].height;
				int size_1 = img[1][0].width * img[1][0].height;
				for (int t = 0; t < num_devices; t++)
				{
					d_state[t].c00 = img[0][t].raw;
					d_state[t].c01 = img[0][t].raw + size_0;
					d_state[t].c02 = img[0][t].raw + 2 * size_0;
					d_state[t].c10 = img[1][t].raw;
					d_state[t].c11 = img[1][t].raw + size_1;
					d_state[t].c12 = img[1][t].raw + 2 * size_1;
				}

				auto get_cost = [w0, h0, mval0, w1, h1, mval1] __device__(int x, int y, State &state)
				{
					C r = C(state.c00[x]) / C(mval0) - C(state.c10[y]) / C(mval1);
					C g = C(state.c01[x]) / C(mval0) - C(state.c11[y]) / C(mval1);
					C b = C(state.c02[x]) / C(mval0) - C(state.c12[y]) / C(mval1);
					C u = C(x % w0) / C(w0 - 1) - C(y % w1) / C(w1 - 1);
					C v = C(x / w0) / C(h0 - 1) - C(y / w1) / C(h1 - 1);
					return r * r + g * g + b * b + u * u + v * v;
				};

				int *rowsol = new int[N2];

				solveCUDA<C, C>(start_time, N1, N2, get_cost, d_state, ws, max_memory_local, rowsol, epsilon, silent);

				for (int t = 0; t < num_devices; t++)
				{
					lapFreeDevice(img[0][t].raw);
					lapFreeDevice(img[1][t].raw);
				}
				delete[] rowsol;
				delete[] img[0];
				delete[] img[1];
				delete[] d_state;
			}
		}
	}
}
