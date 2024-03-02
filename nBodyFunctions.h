#ifndef NBODYFUNCTIONS
#define NBODYFUNCTIONS
#include <string>

void barEvolutionGaussian(); 
//void barEvolutionKalnajs(const std::string & stem, const bool isSelfConsistent, const double littleSigma);

void somaniTrappedOrbits(const double littleSigma, const std::string & outFilename);
void sormaniBoxOrbit(int particleIndex); 

void checkingConservedJacobi();
void orbitSection(const std::string & bodiesFile, const std::string & outFile);  

void kalanajTest();

void spiralTesting();
void calculateDiscAM();


void testingPhaseSpace();


void varyingEpsilon(); 
#endif