#ifndef		__MARCHING_SQUARES_ALGORITHM__
#define		__MARCHING_SQUARES_ALGORITHM__

#include "dimensionalBuffer.h"

struct SPoint
{
	SPoint() : x(-1), y(-1)
	{}
	SPoint(double _x, double _y) : x(_x), y(_y)
	{}
	SPoint(const SPoint& other) : x(other.x), y(other.y)
	{}
	double x;
	double y;
};

enum E_CORNER_VALUES
{
	TOP_LEFT = 1,
	TOP_RIGHT = 1<<1,
	BOTTOM_LEFT = 1<<2,
	BOTTOM_RIGHT = 1<<3
};

enum E_MOVEMENT_DIRECTIONS
{
	UP,
	DOWN,
	LEFT,
	RIGHT,
	NO_DIRECTION
};

class MarchingSquaresAlgorithm
{
public:
	MarchingSquaresAlgorithm(DimensionalBuffer<double>& distance_buffer);
	
	virtual ~MarchingSquaresAlgorithm(){}
	
	virtual bool	generateContour(std::vector<SPoint>& output);
	virtual bool	isInside(int row, int column);
	virtual SPoint	generatePoint(E_MOVEMENT_DIRECTIONS dir);
			void	insertPoint(std::vector<SPoint>& output, SPoint p);
private:

	DimensionalBuffer<double>& m_DistanceBuffer;
	int m_top;
	int m_bottom;
	int m_left;
	int m_right;

};
#endif		__MARCHING_SQUARES_ALGORITHM__
