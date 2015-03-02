#include <iostream>
using namespace std;
// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
using namespace TMath;
// UNIC
#include <UNIC/Units.h>
using namespace UNIC;
#include <MAD/GeCrystal.h>
using namespace MAD;
// infinite thin flat plasma sheet
int main(int argc, char** argv)
{
   double dt=1e-5;//ns
   double dx=10;//nm
   double Ee=100;//volt/cm

   double mean=0, sigma=100,/*nm*/ height=1e13;//1/nm2
   bool norm;

   double epsilon=16.2*8.854187817620e-14;//C/V/cm
   double Q=1.60217657e-19;//C
   GeCrystal Ge;
   double mu_e = Ge.MuE(100)/cm2*volt*s;
   double mu_h = Ge.MuH(100)/cm2*volt*s;

   // initialize arrays
   const int N = 100;
   // p and n are number densities
   double x[2*N+1], p[2*N+1], n[2*N+1], E[2*N+1];
   for (int i=0; i<2*N+1; i++) {
      x[i]=(i-N)*dx;
      p[i]=n[i]=height*Gaus(x[i],mean=0,sigma,norm=kTRUE);
      E[i]=Ee;
   }
   double dp[2*N+1], dn[2*N+1], dE[2*N+1];
   for (int i=0; i<2*N; i++) {
      dp[i] = p[i+1]-p[i];
      dn[i] = n[i+1]-n[i];
      dE[i] = E[i+1]-E[i];
   }
   dp[2*N] = p[2*N]-p[2*N-1];
   dn[2*N] = n[2*N]-n[2*N-1];
   dE[2*N] = E[2*N]-E[2*N-1];

   // output
   TFile *output = new TFile("sheet.root","recreate");
   TTree *t = new TTree("t","time slices");
   t->Branch("x",x,Form("x[%d]/D",2*N+1));
   t->Branch("p",p,Form("p[%d]/D",2*N+1));
   t->Branch("n",n,Form("n[%d]/D",2*N+1));
   t->Branch("E",E,Form("E[%d]/D",2*N+1));
   t->Branch("dp",dp,Form("dp[%d]/D",2*N+1));
   t->Branch("dn",dn,Form("dn[%d]/D",2*N+1));
   t->Branch("dE",dE,Form("dE[%d]/D",2*N+1));

   // evolve
   int iStep=0, nSteps = 10;
   if (argc>1) nSteps = atoi(argv[1]);
   while (iStep<nSteps) {
      t->Fill();
      // update electron and hole distributions
      for (int i=0; i<2*N+1; i++) {
         double dn_dt = mu_e*(dn[i]/dx*E[i]+n[i]*dE[i]/dx);
         n[i]+=dn_dt*dt;
         double dp_dt = -mu_h*(dp[i]/dx*E[i]+p[i]*dE[i]/dx);
         p[i]+=dp_dt*dt;
      }
      // update electric field distribution
      for (int i=0; i<2*N+1; i++) {
         E[i]=0;
         for (int j=0; j<i; j++) E[i]+=p[j]-n[j];
         for (int j=i+1; j<2*N+1; j++) E[i]+=n[j]-p[j];
         E[i]/=(2*epsilon/Q/dx/1e-7);
         E[i]+=Ee;
      }
      // update slopes
      for (int i=0; i<2*N; i++) {
         dp[i] = p[i+1]-p[i];
         dn[i] = n[i+1]-n[i];
         dE[i] = E[i+1]-E[i];
      }
      dp[2*N] = p[2*N]-p[2*N-1];
      dn[2*N] = n[2*N]-n[2*N-1];
      dE[2*N] = E[2*N]-E[2*N-1];
      iStep++;
   }
   t->Write("t",6);
   output->Close();

   return 0;
}

