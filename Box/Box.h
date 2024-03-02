#ifndef BOX
#define BOX 

#include <fftw3.h>
#include <vector>

#include "../Bodies/Bodies.h"

using namespace std; 

class Box {
 public:
  Box(int n, double xmax, double eps, double G);
  ~Box();
  void zero();
  void bodies2density(const Bodies &ptle);

  void bodies2density_m2(const Bodies &ptle, int fourierIndex,  int nring, int nphi);
  void density2pot();
  valarray<double> pot2accels(const Bodies &ptle);

  double potential(double x = 0, double y = 0);  

 //private:
  const double G; // Big G
  const int nx, ny; // Number of points in x and y direction
  const double xmin, ymin; // Minium x and y values, i.e. xmin = -xmax
  const double eps; // Softening
  const double dx, dy, odx, ody; // Grid spacings, 1/gridspacing
  const double xmin1, ymin1; // xmin1 = xmin+dx
  const double xmin2, ymin2; // xmin1 = xmin+2*dx
  const double xmax2, ymax2; // Max values
  
  double *fftw_mesh;
  fftw_complex *fftw_ftmesh, *fftw_ftkernel;
  fftw_plan pk, pkinv;

};

Box::Box(int n, double xmax, double eps, double G)
  : G(G), nx(n), ny(n), eps(eps),
    xmin(-xmax), ymin(-xmax),
    dx(2*xmax/n), dy(2*xmax/n),
    odx(1.0/dx), ody(1.0/dy),
    xmin1(xmin+dx), ymin1(ymin+dy),
    xmin2(xmin+2*dx), ymin2(ymin+2*dy),
    xmax2(-xmin2), ymax2(-ymin2)
{
  fftw_mesh = (double *)
    fftw_malloc((size_t) nx*4*ny*sizeof(double));
  fftw_ftkernel = (fftw_complex*)
    fftw_malloc(sizeof(fftw_complex)*nx*2*(ny+1));
  fftw_ftmesh = (fftw_complex*)
    fftw_malloc(sizeof(fftw_complex)*nx*2*(ny+1));
  
  pk = fftw_plan_dft_r2c_2d(nx*2, ny*2,
                            fftw_mesh, fftw_ftmesh, FFTW_MEASURE);
  pkinv = fftw_plan_dft_c2r_2d(nx*2, ny*2,
                               fftw_ftmesh, fftw_mesh, FFTW_MEASURE);

  double scale = G/(nx*ny*4);
  for(int i=0;i<nx*ny*4;i++) fftw_mesh[i]=0.0;
  for(int i=-nx+1;i<nx;i++)
    for(int j=-ny+1;j<ny;j++) {
      double r=sqrt(eps*eps+i*i*dx*dx+j*j*dy*dy); // softened radius
      int i0=(i+2*nx)%(2*nx), j0=(j+2*ny)%(2*ny);
      int ij=i0*2*ny+j0;
      fftw_mesh[ij]=-scale/r;
    }

  fftw_execute(pk);

  for(int i=0;i<nx*2*(ny+1);i++) {
    fftw_ftkernel[i][0]=fftw_ftmesh[i][0];
    fftw_ftkernel[i][1]=fftw_ftmesh[i][1];
  }
}

Box::~Box() {
  fftw_free(fftw_mesh);
  fftw_free(fftw_ftmesh);
  fftw_free(fftw_ftkernel);
}

void Box::zero() {
  for(int i=0;i<nx*4*ny;i++) fftw_mesh[i] = 0.0;
}

void Box::density2pot() {
  // transform density -> FT(density)
  fftw_execute(pk);

  // Multiply by FT poisson kernel
  fftw_complex *msh = fftw_ftmesh, *krnl = fftw_ftkernel;
  for(int i=0;i<2*nx;i++)
    for(int j=0;j<ny+1;j++) {
      int ij=i*(ny+1)+j;
      double tmp=msh[ij][0];
      msh[ij][0]=krnl[ij][0]*msh[ij][0]-krnl[ij][1]*msh[ij][1];
      msh[ij][1]=krnl[ij][0]*msh[ij][1]+krnl[ij][1]*tmp;
    }

  // Transform back
  fftw_execute(pkinv);
}

void Box::bodies2density(const Bodies &ptle) {
  for(int n=0;n<ptle.n;n++) {
    double x = ptle.xy[2*n];
    double y = ptle.xy[2*n+1];
    double m = ptle.m[n];
    const int DX=2*ny, DY=1;
    if(fabs(x)<xmax2 && fabs(y)<ymax2) {
      int ix = int(odx*(x-xmin));
      int iy = int(ody*(y-ymin));
      int ixy = ix*2*ny+iy;
      double fx = odx*(x-xmin)-ix;
      double fy = ody*(y-ymin)-iy;

      fftw_mesh[ixy]    += (1-fx)*(1-fy)*m;
      fftw_mesh[ixy+DX] += fx*(1-fy)*m;
      fftw_mesh[ixy+DY]    += (1-fx)*fy*m;
      fftw_mesh[ixy+DX+DY] += fx*fy*m;
    }
  }
}



valarray<double> Box::pot2accels(const Bodies &ptle) { // Gives the accels of a load of particles from their accel on a grid
  const int DX=2*ny, DY=1;
  valarray<double> ans(2*ptle.n);
  ans = 0.0;
  for(int n=0;n<ptle.n;n++) {
    double x = ptle.xy[2*n];
    double y = ptle.xy[2*n+1];
    double m = ptle.m[n];
    if(fabs(x)<xmax2 && fabs(y)<ymax2) {
      int ix = int(odx*(x-xmin)); // SHOULD THIS REALLY BE ROUND INSTEAD OF INT? 
      int iy = int(ody*(y-ymin));
      int ixy = ix*2*ny+iy;
      double fx = odx*(x-xmin)-ix;
      double fy = ody*(y-ymin)-iy;
      double *mxy = fftw_mesh+ixy;
      ans[2*n+0] = -0.5*odx*
  ( (1-fx)*(1-fy)*(mxy[+DX]-mxy[-DX])
    +(1-fx)*fy*(mxy[+DX+DY]-mxy[-DX+DY])
    +fx*(1-fy)*(mxy[+2*DX]-mxy[0*DX])
    +fx*fy*(mxy[+2*DX+DY]-mxy[+DY]) );
      ans[2*n+1] = -0.5*ody*
  ( (1-fx)*(1-fy)*(mxy[+DY]-mxy[-DY])
    +(1-fx)*fy*(mxy[+2*DY]-mxy[0*DY])
    +fx*(1-fy)*(mxy[+DX+DY]-mxy[+DX-DY])
    +fx*fy*(mxy[+DX+2*DY]-mxy[+DX]) );
    }
  }
  return ans;
}



void Box::bodies2density_m2(const Bodies &ptle, int fourierIndex, int nring, int nphi) {
  // Fourier analysis: calculate rho_k in shells whose
  // inner boundaries are given by Rmin+i*dR;

  double Rmin = 0.0, Rmax = -xmin;
  double dR = (Rmax-Rmin)/nring;
  double soft2 = 1e-6;
  valarray<double> rhok_cos(nring), rhok_sin(nring);
  rhok_cos = 0.0; rhok_sin = 0.0;

  valarray<double> phis = ptle.angle();
  valarray<double> radii = ptle.radius(soft2);
  
  for(int n=0;n<ptle.n;n++) {
    double cosMphi = cos(fourierIndex * phis[n]);
    double sinMphi = sin(fourierIndex * phis[n]);
    int ndx = ((int) ((radii[n]-Rmin)/dR));
    if(ndx>=0 && ndx<nring) {
      rhok_cos[ndx] += ptle.m[n]*cosMphi;
      rhok_sin[ndx] += ptle.m[n]*sinMphi;
    }
  }


  // pre-calculate trig constants
  valarray<double> cosphi(nphi), sinphi(nphi), cosMphi(nphi), sinMphi(nphi);
  for(int iphi=0;iphi<nphi;iphi++) {
    double phi = iphi*2*M_PI/nphi;
    cosphi[iphi] = cos(phi);
    sinphi[iphi] = sin(phi);
    cosMphi[iphi] = cos(fourierIndex * phi);
    sinMphi[iphi] = sin(fourierIndex * phi);
  }
  double normalisation{1};
  if (fourierIndex ==0) {normalisation = .5;} // Due to the factor of a half in the a_{0} term in defn of Fourier sine and cosine series.

  // now assign Fourier-reconstructed mass to the mesh, ring by ring
  for(int iring=0;iring<nring;iring++){
    double R = Rmin+iring*dR;
    double rhok_cos_here = 2*M_PI*rhok_cos[iring]/(nphi*M_PI);
    double rhok_sin_here = 2*M_PI*rhok_sin[iring]/(nphi*M_PI);
    Bodies ring(nphi);
    for(int iphi=0;iphi<nphi;iphi++) {
      ring.xy[2*iphi+0] = R*cosphi[iphi];
      ring.xy[2*iphi+1] = R*sinphi[iphi];
      ring.m[iphi] = normalisation * rhok_cos_here * cosMphi[iphi] + rhok_sin_here*sinMphi[iphi];
    }
    bodies2density(ring);
  }
}


double Box::potential(double x, double y) {

      double dx{odx*(x-xmin)}, dy{ody*(y-ymin)}; 
      int ix0{(int) floor(dx)}, ix1{(int) floor(dx+1)}, iy0{(int) floor(dy)}, iy1{(int) floor(dy+1)}, 
      ix0y0{ix0*2*ny+iy0}, ix1y0{ix1*2*ny+iy0}, ix0y1{ix0*2*ny+iy1}, ix1y1{ix1*2*ny+iy1}; 

      double pot = fftw_mesh[ix0y0] * (ix1-dx) * (iy1-dy) + fftw_mesh[ix1y0] * (dx-ix0) * (iy1-dy)
      +fftw_mesh[ix0y1] * (ix1-dx) * (dy-iy0) + fftw_mesh[ix1y1] * (dx-ix0) * (dy-iy0);

      return pot; 


      // int ix = round(odx*(x-xmin));
      // int iy = round(ody*(y-ymin));
      // int ixy = ix*2*ny+iy;
      // return fftw_mesh[ixy];
}
#endif 