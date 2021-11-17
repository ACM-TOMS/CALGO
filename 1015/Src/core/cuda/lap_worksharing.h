#pragma once
#include <utility>
#include <vector>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>

namespace lap
{
	namespace cuda
	{
		class Worksharing
		{
		public:
			std::pair<int, int> *part;
			std::vector<int> device;
			std::vector<cudaStream_t> stream;
			std::vector<int> sm_count;
			std::vector<int> threads_per_sm;
			bool silent;
		public:
			Worksharing(int size, int multiple, std::vector<int> &devs, int max_devices, bool silent) : silent(silent)
			{
				max_devices = std::min(max_devices, (size + multiple - 1) / multiple);
				int device_count;

				cudaDeviceProp deviceProp;

				if (devs.empty())
				{
					cudaGetDeviceCount(&device_count);

#ifndef LAP_CUDA_ALLOW_WDDM
					bool allow_wddm = false;
					bool done_searching = false;

					while (!done_searching)
					{
						for (int current_device = 0; ((current_device < device_count) && ((int)device.size() < max_devices)); current_device++)
						{
							cudaGetDeviceProperties(&deviceProp, current_device);



							// If this GPU is not running on Compute Mode prohibited, then we can add it to the list
							if (deviceProp.computeMode != cudaComputeModeProhibited)
							{
								if ((allow_wddm) || (deviceProp.tccDriver))
								{
									if (!silent) lapInfo << "Adding device " << current_device << ": " << deviceProp.name << std::endl;
									if (!silent) lapInfo << "  pciBusID = " << deviceProp.pciBusID << " pciDeviceID = " << deviceProp.pciDeviceID << " pciDomainID = " << deviceProp.pciDomainID << std::endl;
									device.push_back(current_device);
									sm_count.push_back(deviceProp.multiProcessorCount);
									threads_per_sm.push_back(deviceProp.maxThreadsPerMultiProcessor);
								}
							}
						}
						if (device.empty())
						{
							if (allow_wddm) done_searching = true;
							else allow_wddm = true;
						}
						else
						{
							done_searching = true;
						}
					}
#else
					for (int current_device = 0; ((current_device < device_count) && ((int)device.size() < max_devices)); current_device++)
					{
						cudaGetDeviceProperties(&deviceProp, current_device);

						// If this GPU is not running on Compute Mode prohibited, then we can add it to the list
						if (deviceProp.computeMode != cudaComputeModeProhibited)
						{
							if (!silent) lapInfo << "Adding device " << current_device << ": " << deviceProp.name << std::endl;
							device.push_back(current_device);
							sm_count.push_back(deviceProp.multiProcessorCount);
							threads_per_sm.push_back(deviceProp.maxThreadsPerMultiProcessor);
						}
					}
#endif
				}
				else
				{
					device_count = (int)devs.size();
					for (int i = 0; ((i < device_count) && ((int)device.size() < max_devices)); i++)
					{
						int current_device = devs[i];
						cudaGetDeviceProperties(&deviceProp, current_device);

						// If this GPU is not running on Compute Mode prohibited, then we can add it to the list
						if (deviceProp.computeMode != cudaComputeModeProhibited)
						{
							if (!silent) lapInfo << "Adding device " << current_device << ": " << deviceProp.name << std::endl;
							device.push_back(current_device);
							sm_count.push_back(deviceProp.multiProcessorCount);
							threads_per_sm.push_back(deviceProp.maxThreadsPerMultiProcessor);
						}
					}
				}

				if (device.size() == 0)
				{
					std::cout << "No suitable CUDA device found." << std::endl;
					exit(-1);
				}

#ifndef LAP_CUDA_OPENMP
				// single threaded code supports only a single device (TODO: select the fastest device)
				while (device.size() > 1) device.pop_back();
#endif

				int devices = (int)device.size();
				lapAlloc(part, devices, __FILE__, __LINE__);
				for (int p = 0; p < devices; p++)
				{
					long long x0 = (long long)p * (long long)size;
					x0 += devices * multiple - 1;
					x0 /= devices * multiple;
					part[p].first = (int)(multiple * x0);
					if (p + 1 != devices)
					{
						long long x1 = ((long long)p + 1ll) * (long long)size;
						x1 += devices * multiple - 1;
						x1 /= devices * multiple;
						part[p].second = (int)(multiple * x1);
					}
					else part[p].second = size;
				}
				stream.resize(devices);
				for (int t = 0; t < devices; t++)
				{
					cudaSetDevice(device[t]);
					checkCudaErrors(cudaStreamCreateWithFlags(&stream[t], cudaStreamNonBlocking));
				}
			}
			~Worksharing()
			{
				if (part != 0) lapFree(part);
				int devices = (int)device.size();
				for (int t = 0; t < devices; t++)
				{
					cudaSetDevice(device[t]);
					cudaStreamDestroy(stream[t]);
				}
			}
			int find(int x)
			{
				int devices = (int)device.size();
				for (int t = 0; t < devices; t++)
				{
					if ((x >= part[t].first) && (x < part[t].second)) return t;
				}
				return -1;
			}
			bool peerAccess()
			{
				int devices = (int)device.size();
				bool can_access = true;
				for (int d0 = 0; d0 < devices - 1; d0++)
				{
					for (int d1 = d0 + 1; d1 < devices; d1++)
					{
						if (device[d0] != device[d1])
						{
							int canAccessPeer;
							bool access = true;
							cudaDeviceCanAccessPeer(&canAccessPeer, device[d0], device[d1]);
							if (canAccessPeer == 0) access = false;
							cudaDeviceCanAccessPeer(&canAccessPeer, device[d1], device[d0]);
							if (canAccessPeer == 0) access = false;
							if (!silent) lapInfo << "Peer access between device " << device[d0] << " and device " << device[d1] << (access?"":" not") << " possible." << std::endl;
							can_access &= access;
						}
					}
				}
				if (can_access)
				{
					for (int d0 = 0; d0 < devices - 1; d0++)
					{
						for (int d1 = d0 + 1; d1 < devices; d1++)
						{
							if (device[d0] != device[d1])
							{
								cudaSetDevice(device[d0]);
								cudaDeviceEnablePeerAccess(device[d1], 0);
								cudaSetDevice(device[d1]);
								cudaDeviceEnablePeerAccess(device[d0], 0);
							}
						}
					}
				}
				return can_access;
			}
		};
	}
}
