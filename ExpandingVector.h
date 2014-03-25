#ifndef __EXPANDING_VECTOR__
#define __EXPANDING_VECTOR__
#include "FFTWVector.h"
template<class T>
class ExpandingVector : public std::vector<T, FFTW_allocator<T> >
{
	typedef size_t size_type;
	typedef T value_type;
public:
	void resize (size_type n, value_type val = value_type())
	{
		bool needs_to_resize =	n > this->size() || this->size() > 2 * n;

		if(needs_to_resize)
		{
			double log_2 = std::ceil(std::log2((double)n));
			std::vector<T, FFTW_allocator<T> >::resize(1 << (int)log_2);
		}
	}
	void expand(size_type n)
	{
		if(n > this->size())
		{
			this->resize(n);
		}
	}
};
#endif	//__EXPANDING_VECTOR__
