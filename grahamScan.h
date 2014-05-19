#ifndef		__GRAHAM_SCAN_ALGORITHM__
#define		__GRAHAM_SCAN_ALGORITHM__
#include "voro++/include/voro++/voro++.hh"
#include "spoint.h"
#include <vector>
using namespace std;
class GrahamScan
{
public:
	GrahamScan(voro::voronoicell_neighbor& voro_cell, unsigned int cellID, double* partPos);
	GrahamScan(vector<SPoint>& point_cloud);
	void generateCovnexHull(vector<SPoint>& output_hull);
private:
	vector<SPoint>				m_sortedPoints;
	static const SPoint* 		m_compareReference;

	static bool compareByOri(const SPoint& lhs, const SPoint& rhs);
};
#endif		//__GRAHAM_SCAN_ALGORITHM__
