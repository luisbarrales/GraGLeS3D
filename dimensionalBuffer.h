#ifndef		__DIMENSIONAL_BUFFER__
#define		__DIMENSIONAL_BUFFER__

#include <vector>

/*!
 * \class DimensionalBuffer
 * \brief Class that encapsulates a dimensional buffer
 *
 * Dimensional buffer for ease in management of dimensional data. The data exists in space, where
 * the minimum x and y coordinates are 0 and the maximal coordinate value is bounded by the size
 * of signed 32 bit integer.
 */
template<class T>
class	DimensionalBuffer
{
public:
	/*!
	 * \brief Basic constructor. Just initializes the object with some default values.
	 */
	DimensionalBuffer() :
			m_xMin(0), m_xMax(1), m_yMin(0), m_yMax(1)
	{}
	/*!
	* \brief Constructor, which receives the coordinates of the managed region.
	*
	* This constructor takes the boundaries of the managed region as input, and properly
	* initializes the internal memory.
    */
	DimensionalBuffer(unsigned int upperLeftX, unsigned int upperLeftY,
					  unsigned int lowerRightX, unsigned int lowerRightY) :
		m_xMin(upperLeftX), m_yMin(upperLeftY), m_xMax(lowerRightX),
		m_yMax(lowerRightY)
	{
		m_values.resize((m_xMax - m_xMin + 1) * (m_yMax - m_yMin + 1));
	}
	/*!
	* \brief Default destructor.
    */
	~DimensionalBuffer()
	{}
	/*!
	* \brief Method that returns the value at the given coordinates.
	*
	* This method retrieves the value at the specified coordinates. If such a value does not exist
	* this method throws an out of bound exception.
	* \param row the y coordinate of the element.
	* \param column the x coordinate of the element.
	*/
	T getValueAt(unsigned int row, unsigned int column) const
	{
		//Will throw exception if accessed out of bound.
		//TODO: Analyze performance and replace with [] if needed.
		return m_values.at((row - m_yMin) * (m_xMax - m_xMin) + (column - m_xMin));
	}
	/*!
	* \brief Method that sets the value at the given coordinates.
	*
	* This method sets the value at the specified coordinates. If such a value does not exist
	* this method throws an out of bound exception.
	* \param row the y coordinate of the element.
	* \param column the x coordinate of the element.
	* \param value the value to be set.
	*/
	void setValueAt(unsigned int row, unsigned int column, T value)
	{
		//Will throw exception if accessed out of bound.
		//TODO: Analyze performance and replace with [] if needed.
		m_values.at((row - m_yMin) * (m_xMax - m_xMin) + (column - m_xMin)) = value;
	}
	/*!
	* \brief This method resizes the dimensions.
	*
	* This method resizes the dimensions and properly manages the internal data.
	* \param upperLeftX the desired new minimal x coordinate.
	* \param upperLeftY the desired new minimal y coordinate.
	* \param lowerRightX the desired new maximal x coordinate.
	* \param lowerRightY the desired new maximal y coordinate.
	*/
	void resize(unsigned int upperLeftX, unsigned int upperLeftY,
			  	unsigned int lowerRightX, unsigned int lowerRightY)
	{
		m_xMin = upperLeftX;
		m_xMax = lowerRightX;
		m_yMin = upperLeftY;
		m_yMax = lowerRightY;

		m_values.resize((m_xMax - m_xMin) * (m_yMax - m_yMin));
	}

	/*!
	* \brief This method resizes the dimensions.
	*
	* Resizes the current area to a square area and manages the internal memory.
	* \param maximumLength The maximal value for the x or y coordinates.
	*/
	void resizeToSquare(unsigned int maximumLength)
	{
		int grid_size = maximumLength;
		unsigned int height = m_yMax - m_yMin;
		unsigned int width  = m_xMax - m_xMin;
		//First resize the rectangle to a square
		if(width == height)
			return;
		else if (width > height)
		{
			int diff = width - height;
			m_yMin -= diff/2;
			m_yMax += diff/2 + diff%2;
		}
		else
		{
			int diff = height - width;
			m_xMin -= diff/2;
			m_xMax += diff/2 + diff%2;
		}
		//Now move the square in the bounds if it has left them
		if( m_xMin < 0 )
		{
			int delta = -(int)m_xMin;
			m_xMax += delta;
			m_xMin += delta;
		}
		else if (m_xMax > grid_size)
		{
			int delta = grid_size - m_xMax;
			m_xMax += delta;
			m_xMin += delta;
		}
		if( m_yMin < 0 )
		{
			int delta = -(int)m_yMin;
			m_yMax += delta;
			m_yMin += delta;
		}
		else if (m_yMax > grid_size)
		{
			int delta = grid_size - m_yMax;
			m_yMax += delta;
			m_yMin += delta;
		}

		m_values.resize((m_xMax-m_xMin) * (m_yMax-m_yMin));
	}
	/*!
	* \brief This method fills the area with the provided values.
	*/
	void clearValues(T value)
	{
		std::fill(m_values.begin(), m_values.end(), value);
	}
	/*!
	* \brief This method restricts all values to the range of [minimumValue ; maximumValue]
	* \param minimumValue The minimum value.
	* \param maximumValue The maximum value.
	*/
	void clampValues(T minimumValue, T maximumValue)
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
	inline int getMinX() const
	{return m_xMin;}
	inline int getMaxX() const
	{return m_xMax;}
	inline int getMinY() const
	{return m_yMin;}
	inline int getMaxY() const
	{return m_yMax;}
	inline T* getRawData()
	{return &m_values[0];}

private:

	int 	m_xMin;
	int 	m_xMax;
	int 	m_yMin;
	int 	m_yMax;
	std::vector<T>	m_values;
};
#endif		//__DISTANCE_BUFFER__
