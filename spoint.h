#ifndef __SPOINT__
#define	 __SPOINT__
#include <cmath>
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
	double squaredDistanceTo(const SPoint& other) const
	{
		return (x-other.x)*(x-other.x) + (y-other.y)*(y-other.y);
	}
	bool operator==(const SPoint &other) const
	{
		return (x == other.x) && (y == other.y);
	}
	bool operator<(const SPoint& other) const
	{
		return y<other.y || (!(other.y<y) && x<other.x);
	}
	SPoint operator+(const SPoint& other)
	{
		SPoint result(0,0,0);
		result.x = this->x + other.x;
		result.y = this->y + other.y;
		return result;
	}
	SPoint operator-(const SPoint& other)
	{
		SPoint result(0,0,0);
		result.x = this->x - other.x;
		result.y = this->y - other.y;
		return result;
	}
	SPoint operator*(const double other)
	{
		SPoint result;
		result.x = this->x * other;
		result.y = this->y * other;
		return result;
	}
	double dot(const SPoint& other) const
	{
		return x*other.x + y*other.y;
	}
	double len() const
	{
		return sqrt(x*x + y*y);
	}
	double x;
	double y;
	double energy;
};

#endif	//__SPOINT__
