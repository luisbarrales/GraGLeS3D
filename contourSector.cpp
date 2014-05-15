#include "contourSector.h"
#include "spoint.h"

double	ContourSector::INNER_CIRCLE_RADIUS = 0.0;

ContourSector::ContourSector(GrainJunction* initialJunction) :
	m_sectorType(E_CONSTANT_JUNCTION_SECTOR),
	m_leftContourPointID(-1),
	m_rightContourPointID(-1)
{
	m_junctions.push_back(initialJunction);
}

bool ContourSector::mergeWith(GrainJunction* junction)
{
	for(int i=0; i<m_junctions.size(); i++)
	{
		if(m_junctions[i]->coordinates.squaredDistanceTo(junction->coordinates) <= INNER_CIRCLE_RADIUS*INNER_CIRCLE_RADIUS)
		{
			m_junctions.push_back(junction);
			return true;
		}
	}
	return false;
}

bool ContourSector::isPointWithinSectoinRadiuses(const SPoint& point) const
{
	for(int i=0; i<m_junctions.size(); i++)
	{
		if(m_junctions[i]->coordinates.squaredDistanceTo(point) <= INNER_CIRCLE_RADIUS*INNER_CIRCLE_RADIUS)
		{
			return true;
		}
	}
	return false;
}

void ContourSector::setLeftContourPoint(int ID)
{
	m_leftContourPointID = ID;
}
void ContourSector::setRightContourPoint(int ID)
{
	m_rightContourPointID = ID;
}
