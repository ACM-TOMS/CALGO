#pragma once

#include <algorithm>

namespace lap
{
	template <class SC, class TC, class CF>
	class DirectIterator
	{
	protected:
		int dim, dim2;
	public:
		CF &costfunc;
	public:
		DirectIterator(int dim, int dim2, CF &costfunc) : dim(dim), dim2(dim2), costfunc(costfunc) {}
		~DirectIterator() {}

		void getHitMiss(long long &hit, long long &miss) { hit = miss = 0; }

		__forceinline const TC *getRow(int i) { return costfunc.getRow(i); }
	};
}
