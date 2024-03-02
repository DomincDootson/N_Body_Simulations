#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>

#include "Bodies.h"
#include "../../DF_Class/Mestel.h"


void cumulativeDensity(const Mestel & df){
	df.cumulativeDensity("Bodies/cumulativeDensity.out");
}


void massSample(Bodies & ptle){
	std::ifstream inFile; inFile.open("Bodies/cumulativeDensity.out");
	int n; double spacing, normalisation;
	inFile >> n >> spacing >> normalisation;
	inFile.close();
	for (int n=0; n < ptle.n; ++n){ptle.m[n] = normalisation / ((double) ptle.n);} 
}

void radiusSample(Bodies & ptle, const Mestel & df){
	std::vector cumulative = df.readInCumulativeDensity("Bodies/cumulativeDensity.out");
	std::random_device generator; std::uniform_real_distribution<double> uniform(0, 1);

	for (int n=0; n < ptle.n; ++n){
		ptle.xy[2*n] = df.radiusSampling(uniform(generator), cumulative);
	}
}

void phiSample(Bodies & ptle){
	std::random_device generator; std::uniform_real_distribution<double> uniform(0, 2 * M_PI);
	for (int n=0; n < ptle.n; ++n){
		ptle.xy[2*n+1] = uniform(generator);
	}
}

void vRSample(Bodies & ptle, const Mestel & df){
	for (int n = 0; n <ptle.n; ++n){
		ptle.vxvy[2*n] = df.vRSampling(); 
	}
}

void vPhiSample(Bodies & ptle, const Mestel & df){
	for (int n = 0; n <ptle.n; ++n){
		ptle.vxvy[2*n+1] = df.vPhiSampling(ptle.xy[2*n], ptle.vxvy[2*n]);
	}
}

void sampleParticles(Bodies & ptle, const Mestel & df){
	std::cout << "Sampling particles\n";
	cumulativeDensity(df);
	massSample(ptle); 
	radiusSample(ptle, df);
	phiSample(ptle); 
	vRSample(ptle, df);
	vPhiSample(ptle, df); 
	ptle.convert2Cartesian();
}

void writeSamples2file(Bodies & ptle){ 
	std::ofstream out("Bodies/particleSamples.out");
	out << ptle.n << '\n';
	for (int n =0; n < ptle.n; ++n){
		out << ptle.m[n] << " " << ptle.xy[2*n] << " " << ptle.xy[2*n+1] << " " << ptle.vxvy[2*n] << " " << ptle.vxvy[2*n+1] << '\n';
	}
	out.close();
}

int main(int argc, char *argv[]) {
	double nPtle{atof(argv[1])}, sigma{atof(argv[2])};
	std::cout << "Temp of disk: " << sigma << '\n';
	Mestel DF(1, 1, sigma);
	Bodies ptle(nPtle);
	sampleParticles(ptle, DF);
	writeSamples2file(ptle);

	return 0; 
} 