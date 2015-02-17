#include "Units.h"
#include "Constants.h"
using namespace UNIC;

#include "D1.h"
using namespace PLASMA;

void main(int argc, char** argv)
{
   double dt=0.01*ns;
   double dx=0.1*nm;
   double Ee=1000*volt/cm;

   // initialize arrays
   const int N = 1000;
   D1 p(2*N+1), n(2*N+1), E(2*N+1);
   for (int i=-N; i<=N; i++) {
      p[i]=n[i]=gaus(i);
      E[i]=0;
      for (int j=-N; j<=i; j++) {
         dE_dx[j]=q*(p[j]-n[j])/epsilon;
         E[i]+=dE_dx[j]*dx;
      }
      E[i]+=Ee;
   }

   // evolve
   double t=0;
   while (t<100*ns) {
      for (int i=-N; i<=N; i++) {
         // update electron distribution
         dn_dx = (n[i+1]-n[i])/dx;
         dn_dt = mu_e*(dn_dx*E[i]+n[i]*dE_dx[i]);
         n[i]+=dn_dt*dt;

         // update hole distribution
         dp_dx = (p[i+1]-p[i])/dx;
         dp_dt = -mu_h*(dp_dx*E[i]+p[i]*dE_dx[i]);
         p[i]+=dp_dt*dt;

         // update electric field distribution
         E[i]=0;
         for (int j=-N; j<=i; j++) {
            dE_dx[j]=q*(p[j]-n[j])/epsilon;
            E[i]+=dE_dx[j]*dx;
         }
         E[i]+=Ee;
      }
      t+=dt;
   }

   return;
}

double gaus(int)
{
   return 1;
}
