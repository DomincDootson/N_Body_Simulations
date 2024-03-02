#ifndef NBODYBAR
#define NBODYBAR

#include <Eigen/Dense>

#include "NBody.h"

#include "../Perturbation_Grids/BarGrid.h"
#include "../Orbit_Sections/OrbitSections.h"

template <class Tbf>
class NBodyBar : public NBody<Tbf>
{
public:
	NBodyBar(const int nParticles, const int numbTimeSteps, const double timesStep, const Tbf & bf, const Bar2D & bar, const double sigma = 0.35, const double xi = 0.5) : 
	NBody<Tbf>(nParticles, numbTimeSteps, timesStep, bf, xi, sigma),
	m_barGrid(bf, bf.maxRadius(), 10, bar) // This needs to be updated with the potential size
	{}
	~NBodyBar() {}

	void testParticleEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating); 
	void nBodyEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating);

	void orbitSectionPerturbation(const std::string & filename, const bool isSelfConsistent) {this->orbitSections(filename, isSelfConsistent, m_barGrid);}
	void angularMomentumSections(const std::string & filename, const std::string & bodiesFile = "Bodies/particleSamples.csv"); 
	void countTrappedOrbits(double runTime = 750, const std::string & bodiesFile = "Bodies/particleSamples.csv"); 

	
	double countChangedOrbits(); 
	double outputOrbit(int index, const std::string & filename); 

	void savePhaseSpace(const std::string & filename, bool outPerturbed = true) const; 

private:
	BarGrid m_barGrid; 

	void barUpdate(const double time, const double freelyRotating, const int updateBarEvery);
	void foregroundParticleEvolution(const bool isSelfConsistent, const double time);
};

template <class Tbf>
void NBodyBar<Tbf>::testParticleEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating){ 
	std::ofstream out(diskFile); int updateBarEvery{100}; m_barGrid.updateGrid(0); 
	savePhaseSpace("../Plotting/Nbody_Sormani_Data/Phase_Space_Data/Beginning.csv"); 
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		this->outputInfo(time, out);
	    this->backgroundParticleEvolution(false);
	    
		if (m_barGrid.updateGridNow(time,updateBarEvery)) { barUpdate(time*this->m_timeStep, freelyRotating, updateBarEvery);}
		foregroundParticleEvolution(false, time*this->m_timeStep); // Update this function
		if (time == (this->m_numbTimeSteps/2)) {savePhaseSpace("../Plotting/Nbody_Sormani_Data/Phase_Space_Data/Middle.csv");}
	}
	m_barGrid.saveBarEvolution(barFile);
	out.close();
	savePhaseSpace("../Plotting/Nbody_Sormani_Data/Phase_Space_Data/End.csv");
}


template <class Tbf>
void NBodyBar<Tbf>::nBodyEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating){ 
	std::ofstream out(diskFile); int updateBarEvery{100}; m_barGrid.updateGrid(0);
	savePhaseSpace("../Plotting/Nbody_Sormani_Data/Phase_Space_Data/Beginning.csv"); 
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		this->outputInfo(time, out);
	    this->backgroundParticleEvolution(true);
	    
		if (m_barGrid.updateGridNow(time,updateBarEvery)) { barUpdate(time*this->m_timeStep, freelyRotating, updateBarEvery); }
		foregroundParticleEvolution(true, time*this->m_timeStep); // Update this function
		if (time == (this->m_numbTimeSteps/2)) {savePhaseSpace("../Plotting/Nbody_Sormani_Data/Phase_Space_Data/Middle.csv");}

	}
	out.close();
	m_barGrid.saveBarEvolution(barFile);
	savePhaseSpace("../Plotting/Nbody_Sormani_Data/Phase_Space_Data/End.csv");
}

template <class Tbf>
void NBodyBar<Tbf>::barUpdate(const double time, const double freelyRotating, const int updateBarEvery){
	m_barGrid.driftBar(updateBarEvery * this->m_timeStep, freelyRotating);
	m_barGrid.updateGrid(time);

	Eigen::VectorXcd coef = this->m_foreground.responseCoefficents(this->m_basisFunction) - this->m_background.responseCoefficents(this->m_basisFunction);
	m_barGrid.kickBar(updateBarEvery * this->m_timeStep, coef, freelyRotating, time);	
}


template <class Tbf>
void NBodyBar<Tbf>::angularMomentumSections(const std::string & filename,  const std::string & bodiesFile) 
{
	OrbitSections sectionsClass(96, bodiesFile); int minIndex{0}, index{0}, skip{100}; m_barGrid.updateGrid(); 
	
	do 
	{	
		m_barGrid.driftBar(this->m_timeStep, 0);
		if (m_barGrid.updateGridNow(index, skip)) {m_barGrid.updateGrid();} 

		sectionsClass.driftStep(this->m_timeStep);
		std::valarray<double> accels = (this->accelsFromBackground(sectionsClass.m_ptle)+m_barGrid.perturbationAccels(sectionsClass.m_ptle, index * this->m_timeStep)); // update perturbationAccels

		sectionsClass.m_ptle.vxvy += accels * this->m_timeStep;
		sectionsClass.driftStep(this->m_timeStep);

		sectionsClass.angularMomentumSections(m_barGrid.angle());

		//if (minIndex != sectionsClass.minIndex()) {minIndex = sectionsClass.minIndex(); std::cout << "Min Index: " << minIndex << '\n';}  	
 		index +=1;

	} while (sectionsClass.continueSections()); 

	sectionsClass.outputSections(filename);

}

template <class Tbf>
void NBodyBar<Tbf>::countTrappedOrbits(double runTime, const std::string & bodiesFile) {
	OrbitSections sectionsClass(20, bodiesFile); int minIndex{0}, index{0}, skip{100}; m_barGrid.updateGrid(); 
	sectionsClass.setUpforCounting(); 

	for (double time  = 0; time < runTime; time += this->m_timeStep) {
		m_barGrid.driftBar(this->m_timeStep, 0);
		if (m_barGrid.updateGridNow(index, skip)) {m_barGrid.updateGrid();} 

		sectionsClass.driftStep(this->m_timeStep);
		std::valarray<double> accels = (this->accelsFromBackground(sectionsClass.m_ptle)+m_barGrid.perturbationAccels(sectionsClass.m_ptle, time));

		sectionsClass.m_ptle.vxvy += accels * this->m_timeStep;
		sectionsClass.driftStep(this->m_timeStep);

		if (time > 0.25 * runTime) {sectionsClass.countingSections(m_barGrid.angle());}  	
 		index +=1;
 		if (index % 100000 == 0) {std::cout << "Fraction of time completed: " << time / runTime << '\n';}
	} 	
	std::cout << "Fraction of trapped orbits: " << sectionsClass.fractionOfOrbitsTrapped() << '\n';
	std::cout << "Fraction of orbits around bar: " << sectionsClass.fractionOfOrbitsAroundBar(2) << '\n';  
}



template <class Tbf> 
double NBodyBar<Tbf>::countChangedOrbits() {
	int updateBarEvery{100}; m_barGrid.updateGrid(0); // TRY SMALLER updateBarEvery
	std::valarray<double> hasSwitched(0.0, this->m_foreground.n); //Initialise this to 0 
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		if (m_barGrid.updateGridNow(time,updateBarEvery)) { barUpdate(time*this->m_timeStep, 0, updateBarEvery);
			for (int n =0; n < this->m_foreground.n; ++n) {
				if (hasSwitched[n] != 1) {
					if (this->m_foreground.angularMomentum(n) < 0) {hasSwitched[n] = 1;}
				}
			}
		}
		foregroundParticleEvolution(false, time);
	}
	for (int n =0; n < this->m_foreground.n; ++n) {if (hasSwitched[n] == 1) {std::cout << n <<'\n';}}
	return hasSwitched.sum() / ((double) this->m_foreground.n); 
}

template <class Tbf> 
double NBodyBar<Tbf>::outputOrbit(int index, const std::string & filename) {
	int updateBarEvery{100}; m_barGrid.updateGrid(0); 
	std::ofstream out(filename);
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		if (m_barGrid.updateGridNow(time,updateBarEvery)) { barUpdate(time*this->m_timeStep, 0, updateBarEvery);}

		out << time * (this->m_timeStep) << ','; this->m_foreground.particlePosition(out, index); 
		foregroundParticleEvolution(false, time);
	}
	out.close(); 
}

template<class Tbf>
void NBodyBar<Tbf>::savePhaseSpace(const std::string & filename, bool outPerturbed) const {
	if (outPerturbed) {this->m_foreground.savePhaseSpace(filename, m_barGrid.potential(this->m_foreground, this->m_timeStep * this->m_numbTimeSteps), m_barGrid.patternSpeed());}
	else {this->m_background.savePhaseSpace(filename, 0*m_barGrid.potential(this->m_background, this->m_timeStep * this->m_numbTimeSteps), m_barGrid.patternSpeed());}
	std::cout <<"Phase space saved to: " << filename<< '\n';
}

template<class Tbf>
void NBodyBar<Tbf>::foregroundParticleEvolution(const bool isSelfConsistent, const double time) {
	this->m_foreground.xy   += this->m_foreground.vxvy * this->m_timeStep * 0.5;
	std::valarray<double> perturbationAccels = m_barGrid.perturbationAccels(this->m_foreground, time); 
	
	this->m_foreground.vxvy += (this->accelsFromBackground(this->m_foreground) + perturbationAccels) * this->m_timeStep;
	if (isSelfConsistent) {this->m_foreground.vxvy += this->accelsFromDisk(this->m_foreground) * this->m_timeStep;}
	this->m_foreground.xy   += this->m_foreground.vxvy * this->m_timeStep * 0.5;
} 

#endif 


