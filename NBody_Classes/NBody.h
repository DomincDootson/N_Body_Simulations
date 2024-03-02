#ifndef NBODYCLASS
#define NBODYCLASS

#include <Eigen/Dense>
#include <complex>
#include <string>
#include <chrono>

#include "../Bodies/Bodies.h"
#include "../Box/Box.h"
#include "../Orbit_Sections/OrbitSections.h"

#include "../../DF_Class/Mestel.h"
#include "../../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../../Bar2D/Bar2D.h"

template <class Tbf>
class NBody
{
public:
	NBody(const int nParticles, const int numbTimeSteps, const double timesStep, const Tbf & bf, const double xi = 0.5, const double sigma = 0.35) : 
	 m_DF(),
	 m_diskBox(120, 26.0, 0.18,1), m_m0Box(120, 26.0, 0.18,1),
	 m_basisFunction(bf),
	 m_numbTimeSteps{numbTimeSteps}, m_fourierHarmonic{bf.fourierHarmonic()}, m_skip{100}, m_timeStep{timesStep}, m_xi{xi},
	 m_foreground("Bodies/particleSamples.out", nParticles, xi, sigma), 
	 m_background{m_foreground}
	 {if (m_fourierHarmonic == 0){
	 	m0Grid();}
	 	std::cout << "xi is: " << xi << '\n';
	 	std::cout << "Fourier Harmonic: " << m_fourierHarmonic <<'\n';}

	~NBody() {}

	void cumulativeDistribution(const std::string & filename) const {m_DF.cumulativeDensity(filename);}
	
	void barEvolution(const std::string & filename, Bar2D & bar); // derived class

	void particleSampling() {m_foreground.samplingDF();}

	double totalAngularMomentum() const {return m_background.m[0]*m_background.angularMomentum().sum();}
protected:
	
	Box m_diskBox, m_m0Box;

	Mestel m_DF;

	Tbf m_basisFunction; 

	const int m_numbTimeSteps, m_fourierHarmonic, m_skip;
	const double m_timeStep, m_xi;

	Bodies m_foreground, m_background;

	void outputInfo(int timeIndex, std::ofstream & out);
	void outputInfo(int timeIndex) const {if (timeIndex % m_skip == 0) {std::cout << "Fraction of Run: " << timeIndex/((double) m_numbTimeSteps) << '\n';}}; 
	void outputCoefficents(std::ofstream & out);


	void backgroundParticleEvolution(const bool isSelfConsistent); 

	template <class Tgrid>
	void foregroundParticleEvolution(const bool isSelfConsistent, const Tgrid & perturbationGrid);
	void foregroundParticleEvolution(const bool isSelfConsistent);

	valarray<double>  accelsFromBackground(const Bodies & ptle) const;
	valarray<double>  accelsFromDisk(const Bodies & ptle);
	
	void m0Grid();
};

// Output Functions //
// ---------------- //  

std::string outComplexNumber(const std::complex<double> number){
	std::string outString = std::to_string(real(number));
	if (imag(number) < 0){return outString + '-' +std::to_string(abs(imag(number))) +'j';} 
	else {return outString + '+' + std::to_string(imag(number)) + 'j';}
}

template <class Tbf>
void NBody<Tbf>::outputInfo(int timeIndex, std::ofstream & out) {
	if (timeIndex % m_skip == 0) { 
		outputCoefficents(out); 
		std::cout << "Fraction of Run: " << timeIndex/((double) m_numbTimeSteps) << '\n';
	}
}


template <class Tbf>
void NBody<Tbf>::outputCoefficents(std::ofstream & out) {
	Eigen::VectorXcd coef = m_foreground.responseCoefficents(m_basisFunction) - m_background.responseCoefficents(m_basisFunction);
	for (int i =0; i < coef.size()-1; ++i) {out << outComplexNumber(coef(i)) <<',';}
	out << outComplexNumber(coef(coef.size()-1)) << '\n';
}

// Acceleration Functions //
// ---------------------- // 

template <class Tbf>
valarray<double>  NBody<Tbf>::accelsFromBackground(const Bodies & ptle) const {
	std::valarray<double> accels(2*ptle.n);
	for(int i=0;i<ptle.n;i++) 
    {
    	double x{ptle.xy[2*i]}, y{ptle.xy[2*i+1]}, R2{x*x+y*y}, ax{0}, ay{0}, e2{0.25}; 
    	ax += m_DF.xAccel(x, y); ay += m_DF.yAccel(x, y);
    	//ax += -x/(R2 + e2); ay += -y/(R2 +e2); // Do this to include softening 
    	accels[2*i] = ax; accels[2*i+1] = ay;
    }
	return accels;
}

template <class Tbf>
valarray<double>  NBody<Tbf>::accelsFromDisk(const Bodies & ptle){ 
	m_diskBox.zero();
    m_diskBox.bodies2density_m2(ptle, m_fourierHarmonic, 2400, 720); // Put in values for nPhi and nRing
    m_diskBox.density2pot();
    valarray<double> accels = m_diskBox.pot2accels(ptle);
    
    if (m_fourierHarmonic == 0) {accels -= m_m0Box.pot2accels(ptle);}
   
    return accels; 
}
template <class Tbf>
void NBody<Tbf>::backgroundParticleEvolution(const bool isSelfConsistent) {
	m_background.xy   += m_background.vxvy * m_timeStep * 0.5;
	m_background.vxvy += accelsFromBackground(m_background) * m_timeStep;
	if (isSelfConsistent) { m_background.vxvy += accelsFromDisk(m_background) * m_timeStep;}
	m_background.xy   += m_background.vxvy * m_timeStep * 0.5;	
}

template <class Tbf>
template <class Tgrid>
void NBody<Tbf>::foregroundParticleEvolution(const bool isSelfConsistent, const Tgrid & perturbationGrid) {
	m_foreground.xy   += m_foreground.vxvy * m_timeStep * 0.5;
	std::valarray<double> perturbationAccels = perturbationGrid.perturbationAccels(m_foreground); 
	
	m_foreground.vxvy += (accelsFromBackground(m_foreground) + perturbationAccels) * m_timeStep;
	if (isSelfConsistent) {m_foreground.vxvy += accelsFromDisk(m_foreground) * m_timeStep;}
	m_foreground.xy   += m_foreground.vxvy * m_timeStep * 0.5;
}



template <class Tbf>
void NBody<Tbf>::foregroundParticleEvolution(const bool isSelfConsistent) {
	m_foreground.xy   += m_foreground.vxvy * m_timeStep * 0.5;
	
	m_foreground.vxvy += (accelsFromBackground(m_foreground)) * m_timeStep;
	if (isSelfConsistent) {m_foreground.vxvy += accelsFromDisk(m_foreground) * m_timeStep;}
	m_foreground.xy   += m_foreground.vxvy * m_timeStep * 0.5;
}


template <class Tbf>
void NBody<Tbf>::m0Grid() {
	m_m0Box.zero();
	m_m0Box.bodies2density_m2(m_background, 0, 2400, 720);
	m_m0Box.density2pot();
}

#endif