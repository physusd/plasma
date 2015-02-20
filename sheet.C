#include <iostream>
using namespace std;
// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
using namespace TMath;
// UNIC
#include <UNIC/Units.h>
#include <UNIC/Constants.h>
using namespace UNIC;
// infinite thin flat plasma sheet
int main(int argc, char** argv)
{
   int nSteps = atoi(argv[1]);
   double dt=0.1*ns;
   double dx=0.1*nm;
   double Ee=1000*volt/cm;
   double mean=0,sigma=1*nm, height=1/nm;
   bool norm;
   double epsilon=epsilon0;
   double mu_e = 40000*cm2/(volt*s);
   double mu_h = mu_e;
   double dp_dt, dn_dt;

   // initialize arrays
   const int N = 100;
   double x[2*N+1], p[2*N+1], n[2*N+1], V[2*N+1], E[2*N+1];
   double dp,dn,dV,dE;
   for (int i=0; i<2*N+1; i++) {
      x[i]=(i-N)*dx;
      p[i]=n[i]=height*Gaus(x[i],mean=0,sigma,norm=kTRUE);
   }
   for (int i=0; i<2*N+1; i++) {
      V[i]=0;
      for (int j=0; j<2*N+1; j++) {
         if (j==i) continue; // what to do?
         V[i]+=(p[j]-n[j])/Abs(i-j);
      }
      V[i]*=e_SI/4/Pi()/epsilon/dx;
   }
   for (int i=0; i<2*N+1; i++) {
      if (i==2*N) dV = V[i]-V[i-1];
      else dV = V[i+1]-V[i];
      E[i]=dV/dx+Ee;
   }

   TFile *output = new TFile("sheet.root","recreate");
   TTree *t = new TTree("t","time slices");
   t->Branch("x",x,"x[201]/D");
   t->Branch("p",p,"p[201]/D");
   t->Branch("n",n,"n[201]/D");
   t->Branch("E",E,"E[201]/D");

   // evolve
   int iStep=0;
   while (iStep<nSteps) {
      t->Fill();
      // update electron and hole distributions
      for (int i=0; i<2*N+1; i++) {
         if (i==2*N) {
            dp = p[i]-p[i-1];
            dn = n[i]-n[i-1];
            dE = E[i]-E[i-1];
         }
         else {
            dp = p[i+1]-p[i];
            dn = n[i+1]-n[i];
            dE = E[i+1]-E[i];
         }
         dn_dt = mu_e*(dn/dx*E[i]+n[i]*dE/dx);
         n[i]+=dn_dt*dt;
         dp_dt = -mu_h*(dp/dx*E[i]+p[i]*dE/dx);
         p[i]+=dp_dt*dt;
      }

      // update electric potential and field distributions
      for (int i=0; i<2*N+1; i++) {
         V[i]=0;
         for (int j=0; j<2*N+1; j++) {
            if (j==i) continue; // what to do?
            V[i]+=(p[j]-n[j])/Abs(i-j);
         }
         V[i]*=e_SI/4/Pi()/epsilon/dx;
      }
      for (int i=0; i<2*N+1; i++) {
         if (i==2*N) dV = V[i]-V[i-1];
         else dV = V[i+1]-V[i];
         E[i]=dV/dx+Ee;
      }
      iStep++;
   }
   t->Write("t",6);
   output->Close();

   return 0;
}

