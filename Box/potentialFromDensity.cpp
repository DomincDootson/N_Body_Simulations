#include <iostream>
#include <fstream>


#include <fftw3.h>
#include "Box.h"
#include "../Bodies/Bodies.h"

std::vector<double> densities() {
	std::ifstream inFile("../DF_Class/Densities/scarredMestelDensity.csv");
	int nPoint; double spacing; inFile >> nPoint >> spacing;
	std::vector<double> den(nPoint);
	for (int i = 0; i < nPoint; ++i) {inFile >> den[i]; std::cout << i *spacing << " " << den[i] << '\n';}
	inFile.close();

	return den; 
}

double density(const std::vector<double> & den, const double spacing, const double radius) {
	double index{radius/spacing}; int indLower{(int) floor(index)};
	return den[indLower] + (index-indLower) * (den[indLower+1] - den[indLower]); 
}

double massDisk(const std::vector<double> & den, const double spacing) {
	double integral{0};
	for (int i = 0; i < den.size(); ++i) {integral += spacing * 2 * 3.14 * (i*spacing) * den[i];}
	return integral; 
}


Bodies density2Bodies(const std::vector<double> & den, const double spacing, const double maxRadius) {
	int nR{5000}; //nphi = 1 
	double rSpacing{(maxRadius/((double) nR))};
	Bodies ptle(nR);

	for (int i = 0; i < nR; ++i) {
		double radius{i * rSpacing}, denistyAtR{2 * 3.14 * radius * density(den, spacing, radius) * rSpacing}; //
		ptle.m[i] = denistyAtR; ptle.xy[2*i] = radius, ptle.xy[2*i+1] = 0; 
	}
	ptle.convert2Cartesian(); 

	return ptle; 
}



int main() {
	
	std::cout << "Saving background potential.\n"; 
	std::ifstream inFile("/Users/dominicdootson/Documents/PhD/phd/Linear_Stability_Clean/DF_Class/Densities/scarredMestelDensity.csv");
	if (!inFile.good()) {exit(1);}

	int nPoint; double spacing; inFile >> nPoint >> spacing;
	
	std::vector<double> den(nPoint); std::cout << nPoint << '\n';
	for (int i = 0; i < nPoint; ++i) {inFile >> den[i]; std::cout << i << " " << den[i] << '\n';}
	inFile.close(); 

	Bodies ptle = density2Bodies(den, spacing, spacing * nPoint); 

	Box diskBox(240, 26.0, 0.18, 1);
	diskBox.zero();
	diskBox.bodies2density_m2(ptle, 0, 2400, 720);
	diskBox.density2pot();

	std::ofstream out("/Users/dominicdootson/Documents/PhD/phd/Linear_Stability_Clean/DF_Class/Densities/scarredMestelPotential.csv");
	out << nPoint << " " << spacing << '\n';
	for (int i = 0; i < nPoint; ++i) {out << diskBox.potential(i*spacing,0) << '\n';}
	out.close();
	

	return 0; 
}