#include "grahamScan.h"
#include <algorithm>
#include <set>

const SPoint* GrahamScan::m_compareReference = NULL;

GrahamScan::GrahamScan(voro::voronoicell_neighbor& voro_cell, unsigned int cellID, double* partPos) :
		m_cellSampler(voro_cell),
		m_cellID(cellID),
		m_partPos(partPos)
{}
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
	vector<double> vv;
	double x1[2], x2[2];
	set<SPoint>	pointset;
	//Iterate over the voro++ structure and produce unique points on the contour.
	m_cellSampler.vertices(m_partPos[3*(m_cellID-1)],m_partPos[3*(m_cellID-1)+1],m_partPos[3*(m_cellID-1)+2],vv);
	for (int ii = 0; ii < m_cellSampler.p; ii++)
	{
		for (int jj = 0; jj < m_cellSampler.nu[ii]; jj++)
		{

			int k = m_cellSampler.ed[ii][jj];
			x1[0] = vv[3 * ii];
			x1[1] = vv[3 * ii + 1];
			x2[0] = vv[3 * k];
			x2[1] = vv[3 * k + 1];
			pointset.insert(SPoint(x1[1], x1[0], 0));
			pointset.insert(SPoint(x2[1], x2[0], 0));
		}
	}
	vector<SPoint> sorted_points;
	std::copy(pointset.begin(), pointset.end(), std::back_inserter(sorted_points));
	m_compareReference = &(*pointset.begin());
	sort(sorted_points.begin()++, sorted_points.end(), GrahamScan::compareByOri);

	output_hull.push_back(sorted_points[0]);
	output_hull.push_back(sorted_points[1]);
	output_hull.push_back(sorted_points[2]);
	m_compareReference = &output_hull[1];
	int curr_top = 2;
	for(int i=3; i<sorted_points.size(); i++)
	{
		while(compareByOri(output_hull[curr_top], sorted_points[i])==false)
		{
			curr_top--;
			m_compareReference = &output_hull[curr_top-1];
		}
		output_hull.push_back(sorted_points[i]);
		curr_top++;
		m_compareReference = &output_hull[curr_top-1];
	}
	output_hull.push_back(output_hull[0]);
}
