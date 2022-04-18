#pragma once

#include "lap_worksharing.h"

namespace lap
{
	namespace cuda
	{
		template <class SC, class TC, class CF>
		class DirectIterator
		{
		protected:
			int dim, dim2;
		public:
			CF &costfunc;
			Worksharing &ws;

		public:
			DirectIterator(int dim, int dim2, CF &costfunc, Worksharing &ws) : dim(dim), dim2(dim2), costfunc(costfunc), ws(ws) {}
			~DirectIterator() {}

			void getHitMiss(long long &hit, long long &miss) { hit = miss = 0; }

			__forceinline const TC *getRow(int t, int i, bool async) { return costfunc.getRow(t, i); }
			__forceinline bool checkRow(int t, int i) { return true; }
		};
	}
}
