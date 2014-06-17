#ifndef		__CONTOUR_SECTOR__
#define		__CONTOUR_SECTOR__
#include "junction.h"
#include <vector>
enum E_CONTOUR_TYPE
{
	E_INTERPOLATION_SECTOR,
	E_CONSTANT_JUNCTION_SECTOR,
	E_INVALID_TYPE_SECTOR
};
class ContourSector
{
public:
	ContourSector(GrainJunction* initialJunction);
	bool 		isPointInside(int i, int j);
	double		getWeigth(int i, int j);
	bool		mergeWith(GrainJunction* junction);
	bool		isPointWithinSectoinRadiuses(const SPoint& point) const;
	void 		setLeftContourPoint(int ID);
	void 		setRightContourPoint(int ID);
	static		double		INNER_CIRCLE_RADIUS;
public:
	E_CONTOUR_TYPE				m_sectorType;
	int 						m_leftContourPointID;
	int 						m_rightContourPointID;
	std::vector<GrainJunction*>	m_junctions;
};

#endif		//__CONTOUR_SECTOR__