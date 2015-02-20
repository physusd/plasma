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
#include <MAD/GeCrystal.h>
using namespace MAD;
// infinite thin flat plasma sheet
int main(int argc, char** argv)
{
   double dt=0.1*ns;
   double dx=0.1*nm;
   double Ee=1000*volt/cm;

   double mean=0, sigma=1*nm, height=1/nm;
   bool norm;

   GeCrystal Ge;
   double epsilon=Ge.Epsilon()*epsilon0;
   double mu_e = Ge.MuE(100);
   double mu_h = Ge.MuH(100);

   double dp_dt, dn_dt;

   // initialize arrays
   const int N = 100;
   // p and n are number densities in unit area
   double x[2*N+1], p[2*N+1], n[2*N+1], E[2*N+1];
   double dp,dn,dE;
   for (int i=0; i<2*N+1; i++) {
      x[i]=(i-N)*dx;
      p[i]=n[i]=height*Gaus(x[i],mean=0,sigma,norm=kTRUE);
   }
   for (int i=0; i<2*N+1; i++) {
      E[i]=0;
      for (int j=0; j<i; j++) E[i]+=p[j]-n[j];
      for (int j=i+1; j<2*N+1; j++) E[i]+=n[j]-p[j];
      E[i]/=2*epsilon;
      E[i]+=Ee;
   }

   TFile *output = new TFile("sheet.root","recreate");
   TTree *t = new TTree("t","time slices");
   t->Branch("x",x,"x[201]/D");
   t->Branch("p",p,"p[201]/D");
   t->Branch("n",n,"n[201]/D");
   t->Branch("E",E,"E[201]/D");

   // evolve
   int iStep=0, nSteps = 10;
   if (argc>1) nSteps = atoi(argv[1]);
   while (iStep<nSteps) {
      t->Fill();
      // update electron and hole distributions
      for (int i=0; i<2*N+1; i++) {
         if (i==2*N) {
            dp = p[i]-p[i-1];
            dn = n[i]-n[i-1];
            dE = E[i]-E[i-1];
         } else {
            dp = p[i+1]-p[i];
            dn = n[i+1]-n[i];
            dE = E[i+1]-E[i];
         }
         dn_dt = mu_e*(dn/dx*E[i]+n[i]*dE/dx);
         n[i]+=dn_dt*dt;
         dp_dt = -mu_h*(dp/dx*E[i]+p[i]*dE/dx);
         p[i]+=dp_dt*dt;
      }

      // update electric field distribution
      for (int i=0; i<2*N+1; i++) {
         for (int j=0; j<i; j++) E[i]+=p[j]-n[j];
         for (int j=i+1; j<2*N+1; j++) E[i]+=n[j]-p[j];
         E[i]/=2*epsilon;
         E[i]+=Ee;
      }
      iStep++;
   }
   t->Write("t",6);
   output->Close();

   return 0;
}

