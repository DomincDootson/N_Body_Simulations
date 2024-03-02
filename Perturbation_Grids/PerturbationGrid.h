#ifndef PERTURBATIONGRID
#define PERTURBATIONGRID

#include <valarray>

#include "../../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../Bodies/Bodies.h"

class PerturbationGrid
{
public:
	template <class T>
	PerturbationGrid(const T & bf, const double rMaxGrid, const int nGrid) : m_potentialArray(nGrid, nGrid),
	m_spacing{2*rMaxGrid/ ((double) nGrid-1)}, m_rMaxGrid{rMaxGrid}, m_centre{0.5 * (nGrid-1)},
	v_individualPotentials{bf.individualPotential(nGrid, m_rMaxGrid)},
	m_perturbationIndex{0}, m_fourierHarmonic{bf.fourierHarmonic()}
	{} 

	~PerturbationGrid() {}

	virtual bool updateGridNow(const int nBodyIndex, const int totalNsteps) const = 0; 
	virtual void updateGrid() = 0; 

	std::valarray<double> perturbationAccels(const Bodies & ptle) const;

	void saveArray(const std::string & filename) const {
		std::ofstream out(filename); std::cout << "Saving potential to: " << filename << '\n';
		for (int i =0; i < m_potentialArray.rows(); ++i)
			for (int j = 0; j < m_potentialArray.cols(); ++j){
				if (j == m_potentialArray.cols()-1) {out << m_potentialArray(i,j) <<'\n';}
				else {out << m_potentialArray(i,j) << ',';}
			}
	}

	double potential(const double xPosition, const double yPosition) const;
	std::valarray<double> potential(const Bodies & ptle) const;
	
	void saveAccel(const std::string & filename) const { 
		std::ofstream out(filename);

		for (int i = 1; i < m_potentialArray.cols()-1; i++) {
			for (int j =1; j < m_potentialArray.rows()-1; j++) {
				if (j == m_potentialArray.cols()-2) {out << xAccel(i,j) <<'\n';}
				else {out << xAccel(i,j) << ',';}
			}
		}
	}
	

//protected:
		
	Eigen::ArrayXXd m_potentialArray;
	const double m_spacing, m_rMaxGrid, m_centre;
	const std::vector<Eigen::ArrayXXcd> v_individualPotentials; 
	const int m_fourierHarmonic; 

	double m_perturbationIndex;

	bool checkOnGrid(const double xPosition, const double yPosition) const;

	double pos2Index(const double position) const {return position/m_spacing + m_centre;} // 0 < i,j < m_numbgridspace-1 !! Note not strictly

	double xAccel(const int i, const int j) const {return -(1/(2*m_spacing)) * (m_potentialArray(i,j+1) - m_potentialArray(i,j-1));} 
	double xAccel(const double xPosition, const double yPosition) const;

	double yAccel(const int i, const int j) const {return -(1/(2*m_spacing)) * (m_potentialArray(i+1,j) - m_potentialArray(i-1,j));} 
	double yAccel(const double xPosition, const double yPosition) const;
	

	void takeTimeStep() {m_perturbationIndex += 1;}


};


// Functions that get accels //

int floorInt(const double index) {return index/1;}
int ceilInt(const double index) {return (index+1)/1;}

bool PerturbationGrid::checkOnGrid(const double xPosition, const double yPosition) const{
	if ((xPosition*xPosition + yPosition*yPosition) < (m_rMaxGrid-2*m_spacing)*(m_rMaxGrid-2*m_spacing)){
		return true;}
	return false;
}


double PerturbationGrid::xAccel(const double xPosition, const double yPosition) const 
{
	if (!checkOnGrid(xPosition, yPosition)) {return 0;}
	double j{pos2Index(xPosition)}, i{pos2Index(yPosition)};
	int i0{floorInt(i)}, j0{floorInt(j)}, i1{ceilInt(i)}, j1{ceilInt(j)};
	double accel_i0_j0{xAccel(i0, j0)}, accel_i1_j0{xAccel(i1, j0)}, accel_i0_j1{xAccel(i0, j1)}, accel_i1_j1{xAccel(i1, j1)};
	return accel_i0_j0*(i1-i)*(j1-j) + accel_i1_j0*(i-i0)*(j1-j) + accel_i0_j1*(i1-i)*(j-j0) + accel_i1_j1*(i-i0)*(j-j0); //linear interpolation 
}

double PerturbationGrid::yAccel(const double xPosition, const double yPosition) const 
{
	if (!checkOnGrid(xPosition, yPosition)) {return 0;}
	double j{pos2Index(xPosition)}, i{pos2Index(yPosition)};
	int i0{floorInt(i)}, j0{floorInt(j)}, i1{ceilInt(i)}, j1{ceilInt(j)};
	double accel_i0_j0{yAccel(i0, j0)}, accel_i1_j0{yAccel(i1, j0)}, accel_i0_j1{yAccel(i0, j1)}, accel_i1_j1{yAccel(i1, j1)};

	return accel_i0_j0*(i1-i)*(j1-j) + accel_i1_j0*(i-i0)*(j1-j) + accel_i0_j1*(i1-i)*(j-j0) + accel_i1_j1*(i-i0)*(j-j0); //linear interpolation 
}


std::valarray<double> PerturbationGrid::perturbationAccels(const Bodies & ptle) const 
{
	std::valarray<double> accels(2*ptle.n);
	for (int nParticle = 0; nParticle < ptle.n; ++nParticle){
		accels[2*nParticle]   = xAccel(ptle.xy[2*nParticle], ptle.xy[2*nParticle+1]);
		accels[2*nParticle+1] = yAccel(ptle.xy[2*nParticle], ptle.xy[2*nParticle+1]);
	}
	return accels;
}
 
double PerturbationGrid::potential(const double xPosition, const double yPosition) const 
{
	if (!checkOnGrid(xPosition, yPosition)) {return 0;}
	double j{pos2Index(xPosition)}, i{pos2Index(yPosition)};
	int i0{floorInt(i)}, j0{floorInt(j)}, i1{ceilInt(i)}, j1{ceilInt(j)};
	return m_potentialArray(i0,j0) * (i1-i) * (j1-j) + m_potentialArray(i1, j0) * (i-i0) * (j1-j) 
	+ m_potentialArray(i0, j1) * (i1-i) * (j-j0) + m_potentialArray(i1, j1) * (i-i0) * (j-j0); 
}

std::valarray<double> PerturbationGrid::potential(const Bodies & ptle) const {
	std::valarray<double> pot(ptle.n);

	for (int i = 0; i < ptle.n; ++i) {
		pot[i] = potential(ptle.xy[2*i], ptle.xy[2*i+1]);
	}
	return pot; 
}


#endif