#include <iostream>
#include <Eigen/Dense>

#include "NBody_Classes/NBodyPerturbation.h"
#include "NBody_Classes/NBodyBar.h"
#include "NBody_Classes/NBodySpiral.h"


#include "../DF_Class/Mestel.h"
#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"
#include "../Potential_Density_Pair_Classes/GaussianLogBasis.h"

#include "../Bar2D/Bar2D.h"
#include "../Spiral2D/Spiral2D.h"

#include <cmath>

const int NUMBPARTICLES{100000}; //1000000
const int NSTEPS{10000};   //300000
const double TIMESTEP{0.01};


Eigen::VectorXcd gaussianBar(const PotentialDensityPairContainer<GaussianLogBasis> & pd) {
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(pd.maxRadialIndex() + 1);
	for (int i = 0; i < pd.maxRadialIndex() + 1; ++i){
		coeff(i) = 0.001*pd.potential(2.06271, i);
	}
	return - (pd.scriptE()).inverse() * coeff;
}

std::string evolutionFilenameGaussian(int runNumber) {
	return "../Plotting/GaussianTorque/Evolution_" + std::to_string(runNumber) + ".csv";
} 

std::string coefficentFilenameGaussian(int runNumber) {
	return "../Plotting/GaussianTorque/Coefficent_" + std::to_string(runNumber) + ".csv";
} 


std::string evolutionFilenameKalnajs(int runNumber) {
	return "../Plotting/KalnajsTorque/Sormani_Bar/Evolution_" + std::to_string(runNumber) + ".csv";
} 

std::string coefficentFilenameKalnajs(int runNumber) {
	return "../Plotting/KalnajsTorque/Sormani_Bar/Coefficent_" + std::to_string(runNumber) + ".csv";
} 

std::string stemFilename(const std::string & stem, const int runNumber) {
	return "../Plotting/KalnajsTorque/Sormani_Bar/" + stem + "_" + std::to_string(runNumber) + ".csv";
}


void barEvolutionGaussian() 
{
	std::vector<double> params{50, .15, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> pd(params, 50, 2);

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(50+1);
	std::complex<double> unitComplex(0,1);
	

	for (int i =0; i < 5; ++i){
		std::cout << '\n' << '\n' << "Realisation: " << i << '\n';

		coeff = gaussianBar(pd);	
		Bar2D bar(coeff, 0.01, "../Bar2D/barSize.out");

		NBodyBar nbodyBar(1, NSTEPS, TIMESTEP, pd, bar);  
		nbodyBar.testParticleEvolution(coefficentFilenameGaussian(i), evolutionFilenameGaussian(i), 0);
		exit(0); 
		//nbodyBar.nBodyEvolution(coefficentFilenameGaussian(i), evolutionFilenameGaussian(i), 0);
	}
}

void barEvolutionKalnajs(double ep, double omegaP, const std::string & evolutionFile, const bool isSelfConsistent, const double littleSigma)
{
 	std::cout << evolutionFile <<'\n'; 

 	PotentialDensityPairContainer<KalnajsNBasis> pd("../Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_2.dat", 40);
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(40+1); coeff(0) = 0.01; 
	Bar2D bar(coeff, omegaP, "../Bar2D/barSize.out"); // I think we are running this of 0.1 
	bar.sormaniBar(pd, ep, "../Bar2D/Bar_Potentials/Sormani_Large.out");

	NBodyBar nbodyBar(1000000, 10000, 0.01, pd, bar, littleSigma, 0.5);  // 1000000
	if (isSelfConsistent) {nbodyBar.nBodyEvolution("Coeff.csv", evolutionFile, 0);}
	else {nbodyBar.testParticleEvolution("Coeff.csv", evolutionFile, 0);} 
}

void somaniTrappedOrbits(const double littleSigma, const std::string & outFilename)
{
 	int nReal{1}; 
 	PotentialDensityPairContainer<KalnajsNBasis> pd("../Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_2.dat", 70);
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(70+1); coeff(0) = 0.01; 

	std::vector<double> omegaP(10);
	for (int i = 0; i < omegaP.size(); ++i) {omegaP[i] = i * 0.01;}
	std::ofstream out(outFilename); 
	for (auto op : omegaP)
		{
			out<< op <<','; 
			for (int n = 0; n < nReal; ++n) {
				std::cout << op << " " << n <<'\n';
				Bar2D bar(coeff, 0.03, "../Bar2D/barSize.out"); 
				bar.sormaniBar(pd, op, "../Bar2D/Bar_Potentials/Sormani_Medium.out");
			
				NBodyBar nbodyBar(NUMBPARTICLES, NSTEPS, TIMESTEP, pd, bar, littleSigma);  
				if (n < nReal-1) {out << nbodyBar.countChangedOrbits() <<',';}
				else {out<< nbodyBar.countChangedOrbits() <<'\n';}
			}
		}
		out.close(); 
}

void orbitSection(const std::string & sectionBodies, const std::string & outFile) // Orbit sections 
{
	PotentialDensityPairContainer<KalnajsNBasis> pd("../Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_2.dat", 70);
	
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(70+1);
	Bar2D bar(coeff, 0.18, "../Bar2D/barSize.out"); 
	bar.sormaniBar(pd, 0.005, "../Bar2D/Bar_Potentials/Sormani_Large.out");

	NBodyBar nbodyBar(1, 200, 0.0005, pd, bar);
		
	nbodyBar.angularMomentumSections(outFile, sectionBodies);
}

void sormaniBoxOrbit(int particleIndex) {
 	PotentialDensityPairContainer<KalnajsNBasis> pd("../Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_2.dat", 70);
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(70+1); coeff(0) = 0.01; 	
	
	Bar2D bar(coeff, 0.18, "../Bar2D/barSize.out"); 
	bar.sormaniBar(pd, 0.09, "../Bar2D/Bar_Potentials/Sormani_Medium.out");

	NBodyBar nbodyBar(particleIndex+1, NSTEPS, TIMESTEP, pd, bar, 0.35);  
	nbodyBar.outputOrbit(particleIndex, "../Plotting/BoxOrbit_1.csv");
	
}


void calculateDiscAM() {
	std::vector<double> params{4, 20};
 	PotentialDensityPairContainer<KalnajsBasis> pd(params, 10, 2);
 	 NBodyPerturbation nbody(NUMBPARTICLES, NSTEPS, TIMESTEP, pd);
 	 std::cout << "Total angular momentum: " << nbody.totalAngularMomentum() << '\n';
}


/* Outputting Phase Space info */ 
/* --------------------------- */

// Some function that does indidual Run 

void phaseSpaceInfo (const std::string & filename, double ep, double patternSpeed, bool isPerturbed = true) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("../Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_2.dat", 70);
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(70+1); coeff(0) = 0.01; 

	Bar2D bar(coeff, patternSpeed, "../Bar2D/barSize.out"); 
	bar.sormaniBar(pd, ep, "../Bar2D/Bar_Potentials/Sormani_Medium.out");		
	NBodyBar nbodyBar(NUMBPARTICLES, NSTEPS, TIMESTEP, pd, bar, 0.35);  
	nbodyBar.testParticleEvolution("doesntmatter0.csv", "doesntmatter1.csv", 0);

	nbodyBar.savePhaseSpace(filename, isPerturbed); 
}


void testingPhaseSpace() {
	phaseSpaceInfo("Phase_Space_Sormani_0_03_New.csv", 0, 0.03, true);
	phaseSpaceInfo("Phase_Space_Sormani_05_03_New.csv", 0.05, 0.03, true);
	phaseSpaceInfo("Phase_Space_Sormani_05_0_New.csv", 0.05, 0.00, true);
}


/* Some Evolution Functions */ 
/* ------------------------ */ 




void varyingEpsilon() {
	barEvolutionKalnajs(0.005, 0.18, "../Plotting/Nbody_Sormani_Data/Varying_Ep/evolution_005_18_0.csv", false, 0.35);
	barEvolutionKalnajs(0.005, 0.18, "../Plotting/Nbody_Sormani_Data/Varying_Ep/evolution_005_18_1.csv", false, 0.35);
	barEvolutionKalnajs(0.005, 0.18, "../Plotting/Nbody_Sormani_Data/Varying_Ep/evolution_005_18_2.csv", false, 0.35);

	barEvolutionKalnajs(0.005, 0.18, "../Plotting/Nbody_Sormani_Data/Varying_Ep/evolutionW_005_18_0.csv", false, 0.45);
	barEvolutionKalnajs(0.005, 0.18, "../Plotting/Nbody_Sormani_Data/Varying_Ep/evolutionW_005_18_1.csv", false, 0.45);
	barEvolutionKalnajs(0.005, 0.18, "../Plotting/Nbody_Sormani_Data/Varying_Ep/evolutionW_005_18_2.csv", false, 0.45);
	

	
}

