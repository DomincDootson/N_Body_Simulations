#include "nBodyFunctions.h"

#include <iostream>
#include <cmath>

int main(){
	//testingPhaseSpace();
	//barEvolutionKalnajs("Evolution_Test_Cold", false, 0.35); 
	//barEvolutionKalnajs("Evolution_Test_Warm", false, 0.45);

	// barEvolutionKalnajs("Evolution_Consistent_Cold", true, 0.35); 
	// barEvolutionKalnajs("Evolution_Consistent_Warm", true, 0.45);
	
	//somaniTrappedOrbits(0.35, "../Plotting/BO_Amp_realOmega_cold_large.csv");
	//somaniTrappedOrbits(0.45, "../Plotting/BO_Amp_realOmega_warm_large.csv");
	//sormaniBoxOrbit(41);
	//barEvolutionKalnajs("Large_Bar", true, 0.35);

	//orbitSection("Bodies/sormaniCR.out", "sormaniCRSections.csv");
	// orbitSection("Bodies/sormaniNR.out", "sormaniNRSections.csv"); 
	//orbitSection("Bodies/sormaniILR.out", "sormaniILRSections.csv");
	
	
	//barEvolutionGaussian();
	/* 
	barEvolutionKalnajs("Evolution_Self_Cold", true, 0.35); 
	barEvolutionKalnajs("Evolution_Self_Warm", true, 0.45); */
	
	//spiralTesting();

	varyingEpsilon();

	return 0;
} 
