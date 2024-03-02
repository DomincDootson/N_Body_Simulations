#ifndef BODIES
#define BODIES

#include <string>
#include <valarray>
#include <Eigen/Dense>
#include <complex>

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "../../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../../DF_Class/Mestel.h"
#include "../../DF_Class/DFClass.h"
 
class Bodies {
public:
  Bodies(const int N) : m(N), xy(2*N), vxvy(2*N), n(N) {}
  Bodies(const int N, const double xUpper) : m(N), xy(2*N), vxvy(2*N), n(N) // Can we pass a df?? 
  {
    double energy = log(xUpper); 
    for (int i = 0; i < n; ++i) {
      xy[2*i]   = 1 + (xUpper - 1) * (i/((double) n)) ; xy[2*i+1]   = 0;
      vxvy[2*i] = 0; vxvy[2*i+1] = sqrt(2 * (energy - log(xy[2*i])));
      std::cout << xy[2*i] << " " << xy[2*i+1] << " " << vxvy[2*i] << " " << vxvy[2*i+1] << '\n';
    }

  }

  Bodies(const double x, const double y, const double vx, const double vy) : m(1), xy(2), vxvy(2), n(1) {xy[0] = x; xy[1] = y; vxvy[0] = vx; vxvy[1] = vy; m[0] = 1;}


  Bodies(std::string fileName, int numbParticles, const double xi, const double sigma = 0.35);

  ~Bodies() {}
 
  std::valarray<double> m, xy, vxvy;
  int n;
  
  void samplingDF();
  void copy(const Bodies & ptle);



  // some public functions that do the outputting of the coefficents
  std::valarray<double> angularMomentum() const;
  std::valarray<double> angularMomentumCorotating(double omegaP) const {return angularMomentum() - omegaP * radius();} 
  
  std::valarray<double> energy() const;
  std::valarray<double> energy(const std::valarray<double> & pot) const {return energy() + pot;} 

  std::valarray<double> jacobiIntegral(const std::valarray<double> & pot, const double omegaP) const {return (energy() + pot - angularMomentum()* omegaP);} 
  std::valarray<double> jacobiIntegral(const double omegaP) const {return (energy() - angularMomentum()* omegaP);} 

  double angularMomentum(int i) const {return xy[2*i] * vxvy[2*i +1] - xy[2*i +1] * vxvy[2*i];}

  std::valarray<double> radius(double softening2 = 0) const;
  std::valarray<double>  angle() const ;

  template <class Tbf>
  Eigen::VectorXcd responseCoefficents(const Tbf & bf) const;

  void particlePosition(std::ofstream & out, int n = 0) const {out << xy[2*n] <<',' << xy[2*n+1] << ',' << vxvy[2*n] << ',' << vxvy[2*n+1] << '\n';}
  void convert2Cartesian();
  void print(int n) const {std::cout << xy[2*n] << " " << xy[2*n+1] << " " << vxvy[2*n] << " " << vxvy[2*n+1] << '\n';}
  
  std::valarray<double>  apocentre() const;
  std::valarray<double> pericentre() const;

  void savePhaseSpace(const std::string & filename, const std::valarray<double> & pot, const double omegaP) const;

private:
  double  apocentre(const double energy, const double angularMomentum) const;
  double pericentre(const double energy, const double angularMomentum) const;
};


Bodies::Bodies(std::string fileName, int numbParticles, const double xi, const double sigma) : m(numbParticles), xy(2*numbParticles), vxvy(2*numbParticles), n(numbParticles)
{
  std::string argument = "./particleSampling " + std::to_string(numbParticles) + " " + std::to_string(sigma); 
  std::system(argument.c_str());  
  std::cout << "Reading particles in from: " << fileName << '\n';
  std::ifstream inFile; inFile.open(fileName);
  int nFile; 
  inFile >> nFile;
  if (nFile < n) {std::cout << "Particle file doesn't contain enough particles.\n"; std::cout << "It contains: " << nFile <<'\n';exit(1);}

  for (int i = 0; i < n; ++i){
    inFile >> m[i] >> xy[2*i] >> xy[2*i+1] >> vxvy[2*i] >> vxvy[2*i+1]; m[i] *= xi;
  }
  inFile.close(); 
}

void Bodies::copy(const Bodies & ptle) {
  assert(ptle.n >= n && "Trying to copy more values that exist"); // Check that this is the right way around
  for (int i = 0; i < n; i++) {
    m[i] = ptle.m[i];
    xy[2*i] = ptle.xy[2*i]; xy[2*i+1] = ptle.xy[2*i+1];
    vxvy[2*i] = ptle.vxvy[2*i]; vxvy[2*i+1] = ptle.vxvy[2*i+1];
  }
}

template <class Tbf>
Eigen::VectorXcd Bodies::responseCoefficents(const Tbf & bf) const {
  
  Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(bf.maxRadialIndex()+1);
  std::valarray<double> rad{radius(0)}, ang{angle()};
  std::complex<double> unitComplex(0,1);

  for (int i = 0; i < n; ++i){
    for (int j = 0; j < coef.size(); ++j){
      coef[j] += 2*m[i] * exp(-unitComplex*(ang[i]*(bf.fourierHarmonic())))*bf.potential(rad[i], j);
    }
  }
  return -coef;
}



void Bodies::convert2Cartesian() {

  for (int i =0; i < 2*n; i += 2){
    double r{xy[i]}, theta{xy[i+1]}, vr{vxvy[i]}, vPhi{vxvy[i+1]};
    xy[i]     = r * cos(theta); 
    xy[i+1]   = r * sin(theta);
    vxvy[i]   = cos(theta)*vr - sin(theta)*vPhi;
    vxvy[i+1] = sin(theta)*vr + cos(theta)*vPhi;
  }
}

std::valarray<double> Bodies::radius(double softening2) const 
{
  std::valarray<double> radius(n);
  
  for (int i = 0; i<n; i++)
  {
    radius[i] = sqrt(xy[2*i]*xy[2*i] +xy[2*i+1]*xy[2*i+1] +softening2);
  }

  return radius;
}

std::valarray<double>  Bodies::angle() const {
  std::valarray<double> angles(n);

  for (int i = 0; i < n; ++i) {angles[i] = atan2(xy[2*i+1], xy[2*i]);}
  return angles;
}

std::valarray<double>  angle(const std::valarray<double> & xy) {
  int n{((int) xy.size())/2};
  std::valarray<double> angles(n);

  for (int i = 0; i < n; ++i){angles[i] = atan2(xy[2*i+1], xy[2*i]);}
  return angles;
}

std::valarray<double> Bodies::angularMomentum() const
{
  std::valarray<double> angMom(n);
  for (int i = 0; i < n; ++i){
    angMom[i] = xy[2*i]*vxvy[2*i+1] - xy[2*i+1]*vxvy[2*i];
  }
  return angMom;
}


std::valarray<double> Bodies::energy() const
{
  std::valarray<double> energy(n), rad{radius(0)};
  for (int i = 0; i < n; ++i){
    energy[i] = 0.5 * (vxvy[2*i]*vxvy[2*i] + vxvy[2*i+1]*vxvy[2*i+1]) + log(rad[i]);
  }
  return energy;
}

std::valarray<double> Bodies::apocentre() const {
  std::valarray<double> apo(n), en{energy()}, am{angularMomentum()};
  for (int i = 0; i < n; ++i) {apo[n] = apocentre(en[n], am[n]);}
  return apo;
}

std::valarray<double> Bodies::pericentre() const {
  std::valarray<double> peri(n), en{energy()}, am{angularMomentum()};
  for (int i = 0; i < n; ++i) {peri[n] = pericentre(en[n], am[n]);}
  return peri;
}

/* Saving Phase Space */
/* ------------------ */

void Bodies::savePhaseSpace(const std::string & filename, const std::valarray<double> & pot, const double omegaP) const {
  std::ofstream out(filename);
  std::valarray<double> en{energy(pot)}, angMom{angularMomentum()}, jac{jacobiIntegral(pot, omegaP)}, angMomCo{angularMomentumCorotating(omegaP)}; 
  for (int i = 0; i < n; ++i) {out << en[i] << ',' << angMom[i] << ',' << jac[i] << ',' << angMomCo[i] << '\n';}
  out.close();
}



double Bodies::apocentre(const double energy, const double angularMomentum) const {
 
  double x1{pow(M_E, energy-0.5)}, x2{pow(M_E, energy/(1*1))}, holding{}; // Assume vc =1 
  auto radialMometum = [holding, energy, angularMomentum] () {return 2 * holding*holding * (energy - 1*1*log(holding)) - angularMomentum*angularMomentum;};

  for (int i = 0; i < 20; ++i) {
      holding = x1 + 0.5 * (x2 -x1);
      if (radialMometum() > 0) {x1 = holding;}
      else {x2 = holding;}
    }
  return x2;
}

double Bodies::pericentre(const double energy, const double angularMomentum) const {
 
  double x1{0}, x2{pow(M_E, energy-0.5)}, holding{}; 
  auto radialMometum = [holding, energy, angularMomentum] () {return 2 * holding*holding * (energy - 1*1*log(holding)) - angularMomentum*angularMomentum;};

  for (int i = 0; i < 20; ++i) {
      holding = x1 + 0.5 * (x2 -x1);
      if (radialMometum() < 0) {x1 = holding;}
      else {x2 = holding;}
    }
  return x2;
}

#endif 