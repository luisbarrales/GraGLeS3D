#ifndef		__MARCHING_SQUARES_ALGORITHM__
#define		__MARCHING_SQUARES_ALGORITHM__

#include "dimensionalBuffer.h"
#include "ggLS.h"

/*!
 * \struct SPoint
 * \brief Structure used to represent a two dimensional point
 *
 * The point represented by this structure has coordinates of type double. No operators
 * are overloaded for this structure.
 */
struct SPoint
{
	SPoint() : x(-1), y(-1),energy(0)
	{}
	SPoint(double _x, double _y, double _energy) : x(_x), y(_y), energy(_energy)
	{}
	SPoint(const SPoint& other) : x(other.x), y(other.y), energy(other.energy) 
	{}
	double x;
	double y;
	double energy;
};

/*!
 * \enum E_CORNER_VALUES
 * \brief Enumeration that is used to mark the edges of a square that are inside the area
 * of the object whose contour is being traced by the marching squares algorithm.
 */
enum E_CORNER_VALUES
{
	TOP_LEFT = 1,
	TOP_RIGHT = 1<<1,
	BOTTOM_LEFT = 1<<2,
	BOTTOM_RIGHT = 1<<3
};
/*!
 * \enum E_MOVEMENT_DIRECTIONS
 * \brief Enumeration that describes the directions in which the marching square shall advance.
 */
enum E_MOVEMENT_DIRECTIONS
{
	UP,
	DOWN,
	LEFT,
	RIGHT,
	NO_DIRECTION
};
/*!
 * \class MarchingSquaresAlgorithm
 * \brief Class that implements the marching squares algorithm.
 *
 * This class handles 2 dimensional patterns. It uses the standard marching squares algorithm in
 * clock-wise direction. Currently it needs to be improved for ambiguity points. The class has
 * virtual methods that can be overridden to customize the way the algorithm works. The basic
 * implementation is suited to the Levelset project.
 */
class MarchingSquaresAlgorithm
{
public:
	/*!
		 * \brief Basic constructor. Requires a dimensional that defines the shape of the object.
	*/
	MarchingSquaresAlgorithm(DimensionalBuffer<varprecision>& distance_buffer);
	/*!
		 * \brief Basic destructor. Virtual so that the class can be easily inherited.
	*/
	virtual ~MarchingSquaresAlgorithm(){}
	/*!
		 * \brief This method generates the contour as a vector of SPoint objects. In this vector
		 * the first and the last point will be the same. Although this method can be overridden
		 * The core of the algorithm resides here.
	*/
	virtual bool	generateContour(std::vector<SPoint>& output);
	/*!
		 * \brief This method is used to determine whether a point at the specified row and column
		 * is actually inside the object. Can be overridden to customize the contour generation process.
	*/
	virtual bool	isInside(int row, int column);
	/*!
		 * \brief This method generates a point by the current movement direction.
		 * Can be overridden to customize how points on the contour are actually generated.
	*/
	virtual SPoint	generatePoint(E_MOVEMENT_DIRECTIONS dir);
	/*!
		 * \brief This method inserts a point in the output by first checking if the point
		 * is not already there. Prevents point duplication.
	*/
			void	insertPoint(std::vector<SPoint>& output, SPoint p);
private:

	DimensionalBuffer<varprecision>& m_DistanceBuffer;
	int m_top;
	int m_bottom;
	int m_left;
	int m_right;

};
#endif		//__MARCHING_SQUARES_ALGORITHM__
