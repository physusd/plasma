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
// MAD
#include <MAD/GeCrystal.h>
using namespace MAD;
// thin flat infinite plasma sheet
int main(int argc, char** argv)
{
   double dt=1e-1;//ns
   double dx=1;//nm
   double Ee=100;//volt/cm

   double mean=0, sigma=10,/*nm*/ height=25;//1/nm2
   bool norm;

   double epsilon=16.2*8.854187817620e-14;//C/V/cm
   double Q=1.60217657e-19;//C
   double De=100, Dh=50;//cm2/s

   // initialize arrays
   const int N = 100;
   // p and n are number densities
   double x[2*N+1], p[2*N+1], n[2*N+1], E[2*N+1];
   for (int i=0; i<2*N+1; i++) {
      x[i]=(i-N)*dx;
      p[i]=n[i]=height*Gaus(x[i],mean=0,sigma,norm=kTRUE); //1/nm3
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
   double d2p[2*N+1], d2n[2*N+1];
   for (int i=0; i<2*N; i++) {
      d2p[i] = dp[i+1]-dp[i];
      d2n[i] = dn[i+1]-dn[i];
   }
   d2p[2*N] = dp[2*N]-dp[2*N-1];
   d2n[2*N] = dn[2*N]-dn[2*N-1];

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
   t->Branch("d2p",d2p,Form("d2p[%d]/D",2*N+1));
   t->Branch("d2n",d2n,Form("d2n[%d]/D",2*N+1));

   // evolve
   int iStep=0, nSteps = 100;
   if (argc>1) nSteps = atoi(argv[1]);
   double mu_e=40000;
   double mu_h=40000;
   while (iStep<nSteps) {
      t->Fill();
      // update electron and hole distributions
      for (int i=0; i<2*N+1; i++) {
         if (n[i]>1e-3)
            mu_e = 2.14e36/Power(n[i]*1e21,1.82);//cm2/(Vs)
         if (p[i]>1e-3)
            mu_h = 2.14e36/Power(p[i]*1e21,1.82);//cm2/(Vs)
	      double dn_dt = mu_e*(dn[i]/dx*E[i]+n[i]*dE[i]/dx); // 1e7/nm3/s
	      double dp_dt = -mu_h*(dp[i]/dx*E[i]+p[i]*dE[i]/dx);
         dn_dt+=De*d2n[i]/dx/dx/1e7; // 1e7/nm3/s
         dp_dt+=Dh*d2p[i]/dx/dx/1e7;
	      n[i]+=dn_dt*dt;
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
      // update second derivative
      for (int i=0; i<2*N; i++) {
         d2p[i] = dp[i+1]-dp[i];
         d2n[i] = dn[i+1]-dn[i];
      }
      d2p[2*N] = dp[2*N]-dp[2*N-1];
      d2n[2*N] = dn[2*N]-dn[2*N-1];
      iStep++;
   }
   t->Write("t",6);
   output->Close();

   return 0;
}
