// ROOT
#include <TMath.h>
using namespace TMath;
// UNIC
#include <UNIC/Units.h>
#include <UNIC/Constants.h>
using namespace UNIC;
// PLASMA
#include "D1.h"
using namespace PLASMA;
// infinite flat plasma sheet
int main(int argc, char** argv)
{
   int nSteps = atoi(argv[1]);
   double dt=0.001*ns;
   double dx=0.1*nm;
   double Ee=1000*volt/cm;
   double mean=0,sigma=10*nm, height=10/nm;
   bool norm;
   double epsilon=epsilon0;
   double mu_e = 40000*cm2/(volt*s);
   double mu_h = mu_e;
   double dp_dt, dn_dt;

   // initialize arrays
   const int N = 1000;
   D1 x(2*N+1), p(2*N+1), n(2*N+1), V(2*N+1), E(2*N+1);
   for (int i=-N; i<=N; i++) {
      x[i]=i*dx;
      p[i]=n[i]=height*Gaus(x[i],mean=0,sigma,norm=kTRUE);
      V[i]=0;
      for (int j=-N; j<=N; j++) V[i]+=(p[j]-n[j])/Abs(i-j);
      V[i]*=e_SI/4/Pi()/epsilon/dx;
      E[i]=V.Slope(i)/dx+Ee;
   }

   // evolve
   int iStep=0;
   while (iStep<nSteps) {
      for (int i=-N; i<=N; i++) {
         // update electron distribution
         dn_dt = mu_e*(n.Slope(i)/dx*E[i]+n[i]*E.Slope(i)/dx);
         n[i]+=dn_dt*dt;

         // update hole distribution
         dp_dt = -mu_h*(p.Slope(i)/dx*E[i]+p[i]*E.Slope(i)/dx);
         p[i]+=dp_dt*dt;

         // update electric field distribution
         V[i]=0;
         for (int j=-N; j<=N; j++) V[i]+=(p[j]-n[j])/Abs(i-j);
         V[i]*=e_SI/4/Pi()/epsilon/dx;
         E[i]=V.Slope(i)/dx+Ee;
      }
      iStep++;
   }

   return 0;
}

