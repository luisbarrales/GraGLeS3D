#include "contourSector.h"
#include "grahamScan.h"
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
		if(m_junctions[i]->coordinates.squaredDistanceTo(junction->coordinates) <= 4*INNER_CIRCLE_RADIUS*INNER_CIRCLE_RADIUS)
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

void ContourSector::debugPrintSector(vector<SPoint>& contour_grain, ofstream& ofs)
{
	int realContourSize = contour_grain.size()-1;
	int left = m_leftContourPointID;
	int lefter = (left-1+realContourSize)%realContourSize;
	int right = m_rightContourPointID;
	int righter = (right-1+realContourSize)%realContourSize;
	//NEW_CODE
	int width = right < left ? left - right + 1: left + realContourSize - right + 1;
	vector<SPoint> temp_pts;
	temp_pts.resize(2*width);
	for(int i=right, cnt=0; cnt< width; i = (i+1+realContourSize)%realContourSize, cnt++)
	{
		int prev_pnt = (i-1+realContourSize)%realContourSize;
		SPoint v1;
		v1.x = -(contour_grain[i] - contour_grain[prev_pnt]).y;
		v1.y = (contour_grain[i] - contour_grain[prev_pnt]).x;
		v1 = v1*(1/v1.len());
		temp_pts[cnt] = contour_grain[i] + v1*2.0;
		temp_pts[temp_pts.size()-cnt-1] = contour_grain[i] + v1*(-2.0);
	}
	GrahamScan scanner(temp_pts);
	scanner.generateCovnexHull(temp_pts);
	for(int i=0; i<temp_pts.size(); i++)
	{
		ofs<<temp_pts[i].x<<" "<<temp_pts[i].y<<"\n";
	}
	ofs<<"\n";
}
