#ifndef ORBITSECTIONS
#define ORBITSECTIONS

#include <vector>
#include <iostream>
#include <valarray>
#include <algorithm>

#include "../Bodies/Bodies.h"
#include "../Perturbation_Grids/BarGrid.h"

class OrbitSections
{
public:
	OrbitSections(const int nPtle = 20, const std::string & filename = "Bodies/particleSamplesSections.out") : m_nParticles{nPtle}, m_nPoints{140}, m_oldBarAngle{0},
	v_xValues(m_nParticles * m_nPoints), v_pxValues(m_nParticles * m_nPoints), v_indicies(m_nParticles),
	m_ptle(filename, nPtle, 1),
	m_oldXY(2*m_nParticles), m_oldVXVY(2*m_nParticles)
	{}

	~OrbitSections() {}
	
	Bodies m_ptle;  
	
	void driftStep(const double timeStep) {m_ptle.xy += m_ptle.vxvy * timeStep * 0.5;}
	void angularMomentumSections(double barAngle);

	void setUpforCounting(); 
	void countingSections(double barAngle);
	double fractionOfOrbitsTrapped() const {return std::accumulate(v_indicies.begin(), v_indicies.end(), 0)/ ((double) m_nParticles);}
	double fractionOfOrbitsAroundBar(const double barPosition) const; 
	

	bool continueSections() const {if (minIndex() == m_nPoints) {return false;} else {return true;}}
	int minIndex() const { return *std::min_element(v_indicies.begin(), v_indicies.end()); }

	void outputSections(const std::string & filename); 
	void saveAngularMomentum(std::ofstream & out) const; 
	void whichAreTrapped() const {for (const auto& element : v_indicies) {std::cout << element << '\n';}}


private:
	const int m_nParticles, m_nPoints;
	double m_oldBarAngle; 
	std::vector<double> v_xValues, v_pxValues; std::vector<int> v_indicies; 
	std::valarray<double> m_oldXY, m_oldVXVY; 

	bool goneThroughApo(int i) const;
	double pr(const double vx, const double vy, const double ang) const {return vx*cos(ang) + vy*sin(ang);}
	double corotatingAngle(const double barAngle, const int ptleIndex) const; 
	double inertialAngle(const int i) const {return atan2(m_ptle.xy[2*i+1], m_ptle.xy[2*i]);}
	int index(int i) const{return i + v_indicies[i] * m_nParticles;}
};

bool OrbitSections::goneThroughApo(int i) const {
	if ((pr(m_ptle.vxvy[2*i], m_ptle.vxvy[2*i+1], inertialAngle(i)) < 0) && ((pr(m_oldVXVY[2*i], m_oldVXVY[2*i+1], atan2(m_oldXY[2*i+1], m_oldXY[2*i])) > 0))) {return true;}
	else {return false;}
}

double OrbitSections::corotatingAngle(const double barAngle, const int ptleIndex) const {
	double BA{barAngle};
	if (BA > M_PI) {BA -= 2*M_PI;}
	double ang{inertialAngle(ptleIndex)}, holding{ang - (barAngle)};
	if (abs(holding) > M_PI) {holding -= 2 * M_PI * (holding/abs(holding));}
	return holding; 
} 

void OrbitSections::angularMomentumSections(double barAngle) {
	std::valarray<double> angles{m_ptle.angle()}, angMom{m_ptle.angularMomentum()}; barAngle -= 2*M_PI * floor(barAngle/(2*M_PI));
	for (int i = 0; i < m_nParticles; ++i) {
		if (goneThroughApo(i)) {
			v_xValues[index(i)] = corotatingAngle(barAngle, i); v_pxValues[index(i)] = angMom[i]; v_indicies[i] += 1; 
		}
	}
	m_oldXY = m_ptle.xy; m_oldVXVY = m_ptle.vxvy;
}

void OrbitSections::outputSections(const std::string & filename) {
	std::ofstream out(filename);
	for (int i = 0; i < v_xValues.size(); i++){
		 
		if ((i+1) % m_nParticles == 0) {out << v_xValues[i] << ',' << v_pxValues[i] << '\n';}
		else {out << v_xValues[i] << ',' << v_pxValues[i] << ',' ;}
	}
	std::cout << "Sections saved to: " << filename << '\n'; 
	out.close();
}


void OrbitSections::saveAngularMomentum(std::ofstream & out) const {
	std::valarray<double> angMom{m_ptle.angularMomentum()};
	for (int i = 0; i < m_nParticles; ++i) {out << angMom[i] << ',';}
	out << 0 <<'\n';
}


// Counting Sections //
// ----------------- // 

int whichRegion(const double phiApo) { // 0 if abs(\phi) < 0.5 * pi 
	if (abs(phiApo) < 0.5 * M_PI) {return 0;}
	else {return 1;}
}

void OrbitSections::setUpforCounting() {
	v_indicies.resize(m_nParticles);    for (int i = 0; i < v_indicies.size(); ++i) {v_indicies[i] = 1;}
	v_xValues.resize(m_nParticles, 10); for (int i = 0; i <  v_xValues.size(); ++i) { v_xValues[i] = 10;}
}


void OrbitSections::countingSections(double barAngle) {
	std::valarray<double> angles{m_ptle.angle()}; barAngle -= 2*M_PI * floor(barAngle/(2*M_PI));

	for (int i = 0; i < m_nParticles; ++i) {
		if (goneThroughApo(i) && v_indicies[i] != 0) {int areaLabel{whichRegion(corotatingAngle(barAngle, i))};
	
			if (v_xValues[i] == 10) {v_xValues[i] = areaLabel;} // This is for the first time that the ptle reach apo centre. 
			if (v_xValues[i] != areaLabel) {v_indicies[i] = 0;} 
			v_xValues[i] = areaLabel; 
		}
	}
	m_oldXY = m_ptle.xy; m_oldVXVY = m_ptle.vxvy;
	//whichAreTrapped(); std::cout << '\n';
}

double OrbitSections::fractionOfOrbitsAroundBar(const double barPosition) const {
	double count{0}; std::valarray<double> apo{m_ptle.apocentre()}, per{m_ptle.pericentre()};
	auto aroundBar = [barPosition] (double pericentre, double apocentre) {if ((pericentre < barPosition) && (barPosition < apocentre)) {return true;} else {return false;}};

	for (int i = 0; i < m_nParticles; ++i) {if (aroundBar(per[i], apo[i])) {count +=1;}}
	std::cout << count << '\n';
	return count /((double) m_nParticles); 
}

#endif