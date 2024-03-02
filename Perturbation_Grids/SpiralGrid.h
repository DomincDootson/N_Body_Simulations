#ifndef SPIRALGRID
#define SPIRALGRID

#include <Eigen/Dense>

#include "PerturbationGrid.h"
#include "../Bodies/Bodies.h"

#include "../../Spiral2D/Spiral2D.h"

class SpiralGrid : public PerturbationGrid
{
public:
	template <class T>
	SpiralGrid(const T & bf, const double rMaxGrid, const int nGrid, const Spiral2D & spiral) : 
	PerturbationGrid(bf, rMaxGrid, nGrid),
	m_spiral{spiral}
	{}

	~SpiralGrid() {}


	bool updateGridNow(const int nBodyIndex, const int totalNsteps) const {return true;} // This is updating based on the 
	void updateGrid() {} // We need to fill in this function

	Eigen::VectorXcd operator()(int timeIndex) const {return m_spiral(timeIndex);}

	double maxDensity() const {return m_spiral.maxDensity();}
	double density(const double r, const double phi) const {return m_spiral.density(r, phi);}
	double density(const double r) const {return m_spiral.density(r);}

private: 
	Spiral2D m_spiral;
};

#endif