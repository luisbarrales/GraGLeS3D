#ifndef		__CONTOUR_SECTOR__
#define		__CONTOUR_SECTOR__
#include "junction.h"
#include <fstream>
#include <vector>
using namespace std;
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
	void		debugPrintSector(vector<SPoint>& contour_grain, ofstream& ofs);
	static		double		INNER_CIRCLE_RADIUS;

private:
	E_CONTOUR_TYPE				m_sectorType;
	int 						m_leftContourPointID;
	int 						m_rightContourPointID;
	std::vector<GrainJunction*>	m_junctions;
};

#endif		//__CONTOUR_SECTOR__
