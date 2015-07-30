#ifndef		__MARCHING_CUBES_ALGORITHM__
#define		__MARCHING_CUBES_ALGORITHM__

#include "dimensionalBufferReal.h"
#include "triangle.h"
#include <vector>
using namespace std;

class LSbox;

struct GridCell
{
	Vector3d	points[8];
	double		values[8];
};

class MarchingCubesAlgorithm
{
public:
	MarchingCubesAlgorithm(DimensionalBufferReal& distance_buffer, LSbox* current_grain);
	~MarchingCubesAlgorithm(){}
	bool generateHull(vector<Triangle>& triangles);
	bool isInside(int row, int column, int depth);
	unsigned int getIdentifiedNeighborCount() {return m_distinctGrains.size();}
private:

	bool polygonoizeCube(GridCell& inputData, vector<Triangle>& triangles);
	Vector3d VertexInterp(Vector3d p1,Vector3d p2,double valp1,double valp2);
	int generateAdditionalInformation();
	inline void	insertDistinctInteger(int intger, int* array, int& elem_count) const;
	void recordNewGrains(int* grainIDs, int length);

	DimensionalBufferReal& 		m_DistanceBuffer;
	LSbox*						m_currentGrain;

	Vector3i					m_leftBottomFront;
	vector<unsigned int>		m_distinctGrains;
};

#endif		//__MARCHING_CUBES_ALGORITHM__
