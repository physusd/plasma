// std
#include <iostream>
using namespace std;
// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>
using namespace TMath;
// UNIC
#include <UNIC/Units.h>
#include <UNIC/Constants.h>
using namespace UNIC;
// MAD
#include <MAD/GeCrystal.h>
using namespace MAD;
// trace movement of individual charge carriers
int main(int argc, const char *argv[])
{
   double dt=1e-7*ns;
   double E=1000*volt/cm; // external field
   int n = 300; // number of charge carriers
   GeCrystal ge;

   // distribute n e&h's along x with a Gaus PDF
   const int N=10000;
   double xe[N], xh[N], Ee[N], Eh[N], mean, sigma;
   TRandom3 rge(0), rgh(0);
   for (int i=0; i<n; i++) {
      xe[i]=rge.Gaus(mean=0, sigma=10*nm);
      xh[i]=xe[i];
      //xh[i]=rgh.Gaus(mean=0, sigma=10*nm);
   }

   TFile *output = new TFile("carrier.root","recreate");
   TTree *t = new TTree("t","time slices");
   t->Branch("n",&n,"n/I");
   t->Branch("xe",xe,"xe[n]/D");
   t->Branch("xh",xh,"xh[n]/D");
   t->Branch("Ee",Ee,"Ee[n]/D");
   t->Branch("Eh",Eh,"Eh[n]/D");
   double NM=nm, CM=cm, V=volt;
   t->Branch("nm",&NM,"nm/D");
   t->Branch("cm",&CM,"cm/D");
   t->Branch("V",&V,"V/D");

   int iStep=0, nSteps = 100;
   if (argc>1) nSteps = atoi(argv[1]);
   // evolve
   while (iStep<nSteps) {
      // calculate induced E fields on each carrier
      for (int i=0; i<n; i++) { // for the i-th carrier
         Ee[i]=0; Eh[i]=0; // initialize
         for (int j=0; j<n; j++) { // the effect from the j-th carrier
            if (j==i) continue; // skip itself
            Eh[i]-=1/Abs(xh[i]-xe[j])/(xh[i]-xe[j]); // E from electrons
            Eh[i]+=1/Abs(xh[i]-xh[j])/(xh[i]-xh[j]); // E from other holes
            Ee[i]+=1/Abs(xe[i]-xe[j])/(xe[i]-xe[j]); // E from other electrons
            Ee[i]-=1/Abs(xe[i]-xh[j])/(xe[i]-xh[j]); // E from holes
         }
         Eh[i]*=Abs(electron_charge)/4/Pi()/ge.Epsilon()/epsilon0;
         Ee[i]*=Abs(electron_charge)/4/Pi()/ge.Epsilon()/epsilon0;
         Eh[i]+=E; Ee[i]+=E; // plus external E
      }
      t->Fill();
      for (int i=0; i<n; i++) {
         xe[i]-=ge.Mu('e',1e21/cm3)*Ee[i]*dt;
         xh[i]+=ge.Mu('h',1e21/cm3)*Eh[i]*dt;
      }
      iStep++;
   }
   t->Write("t",6);
   output->Close();

   return 0;
}
