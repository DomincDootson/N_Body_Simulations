#ifndef NBODYSPIRAL
#define NBODYSPIRAL

#include <Eigen/Dense>
#include <algorithm>
#include <random>

#include "NBody.h"

#include "../Perturbation_Grids/SpiralGrid.h"

template <class Tbf>
class NBodySpiral : public NBody<Tbf>
{
public:
	NBodySpiral(const int nParticles, const int numbTimeSteps, const double timesStep, const Tbf & bf, const Spiral2D & spiral, const double sigma = 0.35) : 
	NBody<Tbf>(nParticles, numbTimeSteps, timesStep, bf, 0.4, sigma),
	m_spiralGrid(bf, 20, 300, spiral)
	{this->m_foreground = spiralSamplingPosition();}
	~NBodySpiral() {}

	void testParticleEvolution(const std::string & filename);
	void nBodyEvolution(const std::string & filename);

private:
	SpiralGrid m_spiralGrid; 
	std::vector<double> v_density; 

	double mValue() const; 
	double ratio(double r) const {return 1+m_spiralGrid.density(r)/this->m_DF.densityFromVector(r);}


	Bodies spiralSamplingPosition(); 
	bool toKeep(double r, double phi, double comparisionValue, double maxValue) const;
};


template <class Tbf>
void NBodySpiral<Tbf>::testParticleEvolution(const std::string & filename) { 
	std::ofstream out(filename);
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		if (time % this->m_skip == 0) {this->outputCoefficents(out); std::cout << "Fraction of test particle: " << time/((double) this->m_numbTimeSteps) << '\n';}
	    this->backgroundParticleEvolution(false);
		this->foregroundParticleEvolution(false);
	}
	out.close();
}

template <class Tbf>
void NBodySpiral<Tbf>::nBodyEvolution(const std::string & filename) { 
	std::ofstream out(filename);
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		if (time % this->m_skip == 0) {this->outputCoefficents(out); std::cout << "Fraction of NBody: " << time/((double) this->m_numbTimeSteps) << '\n';
		std::cout << this->m_foreground.xy[0] << " " << this->m_background.xy[0] << '\n';
	}
	    this->backgroundParticleEvolution(true);
		this->foregroundParticleEvolution(true);	
	}
	out.close();
}













/* Sampling Functions */ 
/* ------------------ */ 

template <class Tbf>
bool NBodySpiral<Tbf>::toKeep(double r, double phi, double comparisionValue, double maxValue) const {
	double normlisedDensity =(this->m_DF.densityFromVector(r)+m_spiralGrid.density(r, phi))/(maxValue *this->m_DF.densityFromVector(r));
	 
	if (comparisionValue < normlisedDensity) {return true;}
	else if (normlisedDensity != normlisedDensity) {return true;} // Incase the radius is off the grid
	else {return false;}
}

template <class Tbf>
double NBodySpiral<Tbf>::mValue() const {
	double r0{0}, r1{0}, r2{0}, radius{0.20};
	do
	{
		radius +=0.1; 
		r0 = ratio(radius-0.1); r1 = ratio(radius); r2 = ratio(radius+0.1); 
	} while (!(r0 < r1 && r2 < r1));

	return r1; // We assume that the spiral decays quicker than the background 
}


template <class Tbf>
Bodies NBodySpiral<Tbf>::spiralSamplingPosition() {
	
	int nSteps(500); double spacing{0.04};
	this->m_DF.densityVector(nSteps, spacing);

	double maxValue = mValue();

	std::default_random_engine generator; std::uniform_real_distribution<double> distribution(0,1);

	std::valarray<double> angles{this->m_foreground.angle()}, radius{this->m_foreground.radius()};
	std::vector<int> toKeepVector; 

	for (int i = 0; i < this->m_foreground.n; ++i) {
		double r{radius[i]}, phi{angles[i]}, comparisionValue{distribution(generator)}; 
		
		if (toKeep(r, phi, comparisionValue, maxValue)) {toKeepVector.push_back(i);}
	}

	Bodies newPtle(toKeepVector.size());  
	std::cout << toKeepVector.size() << '\n';
	for (int i = 0; i < toKeepVector.size(); ++i) {
		int index{toKeepVector[i]};
		newPtle.xy[2*i] = this->m_foreground.xy[2*index]; newPtle.xy[2*i+1] = this->m_foreground.xy[2*index+1];
		newPtle.vxvy[2*i] = this->m_foreground.vxvy[2*index]; newPtle.vxvy[2*i+1] = this->m_foreground.vxvy[2*index+1];
		newPtle.m[i] = this->m_foreground.m[index] * (((double) this->m_foreground.n)/((double) newPtle.n));
	}
	
	return newPtle;
}


#endif 