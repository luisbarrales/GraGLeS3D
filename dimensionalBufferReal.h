#ifndef		__DIMENSIONAL_BUFFER_REAL__
#define		__DIMENSIONAL_BUFFER_REAL__

#include "dimensionalBuffer.h"

class DimensionalBufferReal : public DimensionalBuffer<float>
{
public:

	DimensionalBufferReal() : DimensionalBuffer<float>()
		{}
	DimensionalBufferReal(unsigned int upperLeftX, unsigned int upperLeftY,
					  unsigned int lowerRightX, unsigned int lowerRightY) :
						  DimensionalBuffer(upperLeftX, upperLeftY, lowerRightX, lowerRightY)
	{}

	void clampValues(double minimumValue, double maximumValue)
	{
		for(int i=0; i<m_values.size(); i++)
		{
			if ( m_values[i] < minimumValue )
			{
				m_values[i] = minimumValue; continue;
			}
			if (m_values[i] > maximumValue )
			{
				m_values[i] = maximumValue; continue;
			}
		}
	}
};
#endif		//__DIMENSIONAL_BUFFER_REAL__
