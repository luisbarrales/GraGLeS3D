#ifndef __POOLED_DIMENSIONAL_BUFFER_REAL__
#define __POOLED_DIMENSIONAL_BUFFER_REAL__
class PooledDimensionalBufferReal
{
public:
	PooledDimensionalBufferReal(char* pool, unsigned int size,
			unsigned int upperLeftX, unsigned int upperLeftY,
			unsigned int lowerRightX, unsigned int lowerRightY) :
				m_xMin(upperLeftX), m_xMax(lowerRightX), m_yMin(upperLeftY), m_yMax(lowerRightY),
				m_pool(pool), m_poolSize(size)
	{
	}
	float getValueAt(unsigned int row, unsigned int column)
	{
		float* pointer = (float*) m_pool;
		return pointer[(row - m_yMin) * (m_xMax - m_xMin) + (column - m_xMin)];
	}

	void setValueAt(unsigned int row, unsigned int column, float value)
	{
		float* pointer = (float*) m_pool;
		pointer[(row - m_yMin) * (m_xMax - m_xMin) + (column - m_xMin)] = value;
	}
	void clearValues(float value)
	{
		for(int i = 0; i < m_poolSize / sizeof(float); i++)
		{
			float* pointer = (float*)m_pool;
			m_pool[i] = value;
		}
	}

private:
	char*	m_pool;
	int 	m_poolSize;
	int 	m_xMin;
	int 	m_xMax;
	int 	m_yMin;
	int 	m_yMax;
};
#endif 	//__POOLED_DIMENSIONAL_BUFFER_REAL__