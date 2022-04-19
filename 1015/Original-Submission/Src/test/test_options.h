#pragma once

#include <string.h>
#include <stdlib.h>
#include <vector>
#include <string>

class Options
{
public:
	long long lap_min_tab;
	long long lap_max_tab;
	long long lap_min_rank;
	long long lap_max_rank;
	long long lap_min_cached;
	long long lap_max_cached;
	long long lap_max_memory;

	bool use_double;
	bool use_float;
	bool use_int;
	bool use_long;
	bool use_single;
	bool use_epsilon;
	bool use_omp;

	bool run_sanity;
	bool run_sanity_cached;
	bool run_random;
	bool run_random_low_rank;
	bool run_random_low_rank_cached;
	bool run_geometric;
	bool run_geometric_cached;
	bool run_geometric_disjoint;
	bool run_geometric_disjoint_cached;
	bool run_integer;
	std::vector<std::string> images;
	std::vector<int> devices;
	bool silent;

	int runs;
public:
	Options()
	{
		lap_min_tab = lap_max_tab = lap_min_cached = lap_max_cached = lap_min_rank = lap_max_rank = lap_max_memory = 0ll;
		use_double = use_float = use_int = use_long = use_single = use_epsilon = use_omp = false;
		run_sanity = run_random = run_geometric = run_geometric_disjoint = run_random_low_rank = run_integer = false;
		run_sanity_cached = run_geometric_cached = run_geometric_disjoint_cached = run_random_low_rank_cached = false;
		runs = 1;
		silent = false;
	}
public:
	void setDefaultSize()
	{
		lap_min_tab = 1000ll;
		lap_max_tab = 64000ll;
		lap_min_cached = 64000ll;
		lap_max_cached = 256000ll;
		lap_max_memory = 4ll * 1024ll * 1024ll * 1024ll;
	}
	void setDefault()
	{
		setDefaultSize();
		use_double = true;
		use_single = true;
		use_epsilon = true;
		use_omp = true;
		run_random = true;
		run_geometric = true;
		run_geometric_cached = true;
		run_geometric_disjoint = true;
		run_geometric_disjoint_cached = true;
		run_integer = false;
	}

	int parseOptions(int argc, char* argv[])
	{
		if (argc == 1)
		{
			setDefault();
		}
		for (int i = 1; i < argc; i++)
		{
			if (!strcmp(argv[i], "-default"))
			{
				setDefault();
			}
			else if (!strcmp(argv[i], "-default_size"))
			{
				setDefaultSize();
			}
			else if (!strcmp(argv[i], "-table_min"))
			{
				lap_min_tab = atoll(argv[++i]);
			}
			else if (!strcmp(argv[i], "-table_max"))
			{
				lap_max_tab = atoll(argv[++i]);
			}
			else if (!strcmp(argv[i], "-rank_min"))
			{
				lap_min_rank = atoll(argv[++i]);
			}
			else if (!strcmp(argv[i], "-rank_max"))
			{
				lap_max_rank = atoll(argv[++i]);
			}
			else if (!strcmp(argv[i], "-cached_min"))
			{
				lap_min_cached = atoll(argv[++i]);
			}
			else if (!strcmp(argv[i], "-cached_max"))
			{
				lap_max_cached = atoll(argv[++i]);
			}
			else if (!strcmp(argv[i], "-memory"))
			{
				lap_max_memory = atoll(argv[++i]);
			}
			else if (!strcmp(argv[i], "-double"))
			{
				use_double = true;
			}
			else if (!strcmp(argv[i], "-float"))
			{
				use_float = true;
			}
			else if (!strcmp(argv[i], "-long"))
			{
				use_long = true;
			}
			else if (!strcmp(argv[i], "-int"))
			{
				use_int = true;
			}
			else if (!strcmp(argv[i], "-single"))
			{
				use_single = true;
			}
			else if (!strcmp(argv[i], "-epsilon"))
			{
				use_epsilon = true;
			}
			else if (!strcmp(argv[i], "-sanity"))
			{
				run_sanity = true;
			}
			else if (!strcmp(argv[i], "-sanity_cached"))
			{
				run_sanity_cached = true;
			}
			else if (!strcmp(argv[i], "-random"))
			{
				run_random = true;
			}
			else if (!strcmp(argv[i], "-random_low_rank"))
			{
				run_random_low_rank = true;
			}
			else if (!strcmp(argv[i], "-random_low_rank_cached"))
			{
				run_random_low_rank_cached = true;
			}
			else if (!strcmp(argv[i], "-geometric"))
			{
				run_geometric = true;
			}
			else if (!strcmp(argv[i], "-geometric_cached"))
			{
				run_geometric_cached = true;
			}
			else if (!strcmp(argv[i], "-geometric_disjoint"))
			{
				run_geometric_disjoint = true;
			}
			else if (!strcmp(argv[i], "-geometric_disjoint_cached"))
			{
				run_geometric_disjoint_cached = true;
			}
			else if (!strcmp(argv[i], "-integer"))
			{
				run_integer = true;
			}
			else if (!strcmp(argv[i], "-silent"))
			{
				silent = true;
			}
			else if (!strcmp(argv[i], "-device"))
			{
				devices.push_back(atoi(argv[++i]));
			}
			else if (!strcmp(argv[i], "-omp"))
			{
#ifdef _OPENMP
				use_omp = true;
#else
				std::cout << "OpenMP not enabled." << std::endl;
#endif
			}
			else if (!strcmp(argv[i], "-img"))
			{
				images.push_back(argv[++i]);
			}
			else if (!strcmp(argv[i], "-runs"))
			{
				runs = atoi(argv[++i]);
			}
			else
			{
				std::cout << "Unkown option: " << argv[i] << std::endl;
				return -1;
			}
		}
		return 0;
	}
};
