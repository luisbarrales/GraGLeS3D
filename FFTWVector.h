#ifndef __FFTW_VECTOR__
#define __FFTW_VECTOR__
#include <fftw3.h>
#include <vector>
/*!
 * \cond HIDDEN_SYMBOLS
 */
template<class T>
class FFTW_allocator
{
public:
	typedef size_t size_type;
	typedef T* pointer;
	typedef const T* const_pointer;
	typedef T& reference;
	typedef const T& const_reference;
	typedef T value_type;

	FFTW_allocator()
	{
	}

	pointer allocate(size_type n, const void * = 0)
	{
		pointer t = (pointer) fftw_malloc(n * sizeof(T));
		return t;
	}

	void deallocate(void* p, size_type)
	{
		if (p)
		{
			fftw_free(p);
		}
	}

	pointer address(reference x) const
	{
		return &x;
	}
	const_pointer address(const_reference x) const
	{
		return &x;
	}
	void construct(pointer p, const_reference val)
	{
		new ((pointer) p) T(val);
	}
	void construct(pointer p)
	{
		new ((pointer) p) T();
	}
	void destroy(pointer p)
	{
		p->~T();
	}

	size_type max_size() const
	{
		return size_type(-1);
	}

	template<class U>
	struct rebind
	{
		typedef FFTW_allocator<U> other;
	};

	template<class U>
	FFTW_allocator(const FFTW_allocator<U>&)
	{
	}

	template<class U>
	FFTW_allocator& operator=(const FFTW_allocator<U>&)
	{
		return *this;
	}
};
/*!
 * \endcond
 */

/*!
 * \class FFTW_Utils
 * \brief Structure that stores utility objects for working with FFTW.
 *	\n Example usage is <b>FFTW_Utils</b>< <em>Desired Type</em> ><b>::FFTW_Vector</b> <em>Desired Name</em>
 */
template<class T>
struct FFTW_Utils
{
	typedef std::vector<T, FFTW_allocator<T> > FFTW_Vector; /*!< Class for FFTW memory aligned vector. */
};

#endif	//__FFTW_VECTOR__
