#ifndef		__CONTOUR_SECTOR__
#define		__CONTOUR_SECTOR__
#include "junction.h"
#include <vector>
enum E_CONTOUR_TYPE
{
	E_CONSTANT_LINE_SECTOR,
	E_INTERPOLATION_SECTOR,
	E_CONSTANT_JUNCTION_SECTOR,
	E_INVALID_TYPE_SECTOR
};
class ContourSector
{
public:
	ContourSector();
	bool 		isPointInside(int i, int j);
	double		getWeigth(int i, int j);
private:
	E_CONTOUR_TYPE			m_sectorType;
	int 					m_leftContourPointID;
	int 					m_rightContourPointID;
	std::vector<GrainJunction*>	m_junctions;
};

#endif		//__CONTOUR_SECTOR__
