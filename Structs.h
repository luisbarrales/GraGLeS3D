class myQuaternion;


struct VolEvolution {
	double dA;
	int nVertex;

	VolEvolution(double da, int nr) :
			dA(da), nVertex(nr) {
	}
};

struct Face {
	double area;
	unsigned int grainA;
	unsigned int grainB;
	Face(double _area, unsigned int _grainA, unsigned int _grainB) :
			area(_area), grainA(_grainA), grainB(_grainB) {
	}
};

struct TextureData {
	unsigned int id;
	unsigned int NeighbourCount;
	unsigned int intersectsBoundaryGrain;
	double volume;
	double surfaceArea;
	double GBEnergy;
	double BulkEnergy;
	double euler[3];
	double barycenter[3];
	TextureData(unsigned int _id, unsigned int _NeighbourCount,
			bool _intersectsBoundaryGrain, double _volume, double _dA,
			double _surfaceArea, double _GBEnergy, double _BulkEnergy,
			myQuaternion *ori) :
			id(_id), NeighbourCount(_NeighbourCount), intersectsBoundaryGrain(
					_intersectsBoundaryGrain), volume(_volume), surfaceArea(
					_surfaceArea), GBEnergy(_GBEnergy), BulkEnergy(_BulkEnergy) {
		double* newori = ori->Quaternion2Euler();
		*euler = *newori;
		delete [] newori;

	}
	~TextureData() {
		delete[] euler;
	}
};
