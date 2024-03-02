#ifndef CHOOSENPERTURBATIONGRID
#define CHOOSENPERTURBATIONGRID

#include "PerturbationGrid.h"
#include "../../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../../Volterra_Solver/ExpansionCoeff.h"


class ChoosenPerturbationGrid : public PerturbationGrid
{
public:
	template <class T>
	ChoosenPerturbationGrid(const T & bf, const double rMaxGrid, const int nGrid, const std::string & filename, const int numbTimeStep) : 
	PerturbationGrid(bf, rMaxGrid, nGrid),
	m_perturbationCoeff(filename, numbTimeStep)
	{}
	~ChoosenPerturbationGrid() {}

	bool updateGridNow(const int nBodyIndex, const int totalNsteps) const; 
	void updateGrid(); 

private:
	ExpansionCoeff m_perturbationCoeff;
};


bool ChoosenPerturbationGrid::updateGridNow(const int nBodyIndex, const int totalNsteps) const {
	
	if (nBodyIndex % (totalNsteps/m_perturbationCoeff.nTimeStep()) == 0  && m_perturbationIndex < m_perturbationCoeff.nTimeStep()){
		return true;
	}
	return false;
}


void ChoosenPerturbationGrid::updateGrid() {
	Eigen::ArrayXXcd potential = Eigen::ArrayXXcd::Zero(m_potentialArray.rows(), m_potentialArray.cols());
	for (int i = 0; i < v_individualPotentials.size(); ++i){ potential += m_perturbationCoeff(m_perturbationIndex)[i] * v_individualPotentials[i];}
	if (m_fourierHarmonic == 0){m_potentialArray = potential.real();}
	else {m_potentialArray = 2 * potential.real();}
	takeTimeStep();
}

#endif