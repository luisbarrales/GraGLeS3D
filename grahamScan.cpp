#include "grahamScan.h"
#include <algorithm>
#include <set>

const SPoint* GrahamScan::m_compareReference = NULL;

GrahamScan::GrahamScan(voro::voronoicell_neighbor& voro_cell, unsigned int cellID, double* partPos)
{
	vector<double> vv;
	double x1[2], x2[2];
	set<SPoint>	pointset;
	//Iterate over the voro++ structure and produce unique points on the contour.
	voro_cell.vertices(partPos[3*(cellID-1)],partPos[3*(cellID-1)+1],partPos[3*(cellID-1)+2],vv);
	for (int ii = 0; ii < voro_cell.p; ii++)
	{
		for (int jj = 0; jj < voro_cell.nu[ii]; jj++)
		{

			int k = voro_cell.ed[ii][jj];
			x1[0] = vv[3 * ii];
			x1[1] = vv[3 * ii + 1];
			x2[0] = vv[3 * k];
			x2[1] = vv[3 * k + 1];
			pointset.insert(SPoint(x1[1], x1[0], 0));
			pointset.insert(SPoint(x2[1], x2[0], 0));
		}
	}
	std::copy(pointset.begin(), pointset.end(), std::back_inserter(m_sortedPoints));
	m_compareReference = &(*pointset.begin());
	vector<SPoint>::iterator start = m_sortedPoints.begin(); start++;
	sort(start, m_sortedPoints.end(), GrahamScan::compareByOri);
}
GrahamScan::GrahamScan(vector<SPoint>& point_cloud)
{
	std::copy(point_cloud.begin(), point_cloud.end(), std::back_inserter(m_sortedPoints));
	sort(m_sortedPoints.begin(), m_sortedPoints.end());
	m_compareReference = &m_sortedPoints[0];
	vector<SPoint>::iterator start = m_sortedPoints.begin(); start++;
	sort(start, m_sortedPoints.end(), GrahamScan::compareByOri);
}
bool GrahamScan::compareByOri(const SPoint& lhs, const SPoint& rhs)
{
	double val = (lhs.y - m_compareReference->y) * (rhs.x - lhs.x) -
	              (lhs.x - m_compareReference->x) * (rhs.y - lhs.y);
	if(val == 0)
	{
		return m_compareReference->squaredDistanceTo(lhs) < m_compareReference->squaredDistanceTo(rhs);
	}
	else
	{
		return val < 0;
	}
}
void GrahamScan::generateCovnexHull(vector<SPoint>& output_hull)
{
	output_hull.clear();
	output_hull.push_back(m_sortedPoints[0]);
	output_hull.push_back(m_sortedPoints[1]);
	output_hull.push_back(m_sortedPoints[2]);
	m_compareReference = &output_hull[1];
	int curr_top = 2;
	for(int i=3; i<m_sortedPoints.size(); i++)
	{
		while(compareByOri(output_hull[curr_top], m_sortedPoints[i]) == false)
		{
			curr_top--;
			output_hull.pop_back();
			m_compareReference = &output_hull[curr_top-1];
		}
		output_hull.push_back(m_sortedPoints[i]);
		curr_top++;
		m_compareReference = &output_hull[curr_top-1];
	}
	output_hull.push_back(output_hull[0]);
}
