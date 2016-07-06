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
	bool intersectsBoundaryGrain;
	double volume;
	double dA;
	double perimeter;
	double GBEnergy;
	double BulkEnergy;
	double* euler;
	TextureData(unsigned int _id, unsigned int _NeighbourCount,
			bool _intersectsBoundaryGrain, double _volume, double _dA,
			double _perimeter, double _GBEnergy, double _BulkEnergy,
			myQuaternion *ori) :
			id(_id), NeighbourCount(_NeighbourCount), intersectsBoundaryGrain(
					_intersectsBoundaryGrain), volume(_volume), dA(_dA), perimeter(
					_perimeter), GBEnergy(_GBEnergy), BulkEnergy(_BulkEnergy) {
		euler = new double(3);
		euler = ori->Quaternion2Euler();

	}
	~TextureData() {
		delete[] euler;
	}
};
