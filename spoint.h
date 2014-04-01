#ifndef __SPOINT__
#define	 __SPOINT__
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

#endif	//__SPOINT__
