#include <iostream>
using namespace std;
// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TGraph.h>
using namespace TMath;
// UNIC
#include <UNIC/Units.h>
#include <UNIC/Constants.h>
using namespace UNIC;
// MAD
#include <MAD/GeCrystal.h>
using namespace MAD;
// thin plasma tube
//GeCrystal *ge = new GeCrystal;
int main(int argc, char** argv)
{
   double dt=1e-7*ns; // large dt causes E to decrease too fast
   double dx=0.5*nm; // large dx may cause asymmetry 
   double Ee=1000*volt/cm;
   double R=1e-6*cm, mean=0, sigma=R/3, height=300./(Pi()*R*R);
   bool norm;

   GeCrystal ge;
   double epsilon=ge.Epsilon()*epsilon0;
   double Q=Abs(electron_charge);
   double De=0.046*cm2/s, Dh=0.046*cm2/s;

   // initialize arrays
   const int N = 500;
   // p and n are number densities
   double x[2*N+1], p[2*N+1], n[2*N+1], E[2*N+1], pE[2*N+1],nE[2*N+1];
   for (int i=0; i<2*N+1; i++) {
      x[i]=(i-N)*dx;
      p[i]=n[i]=height*Gaus(x[i],mean=0,sigma,norm=kTRUE);
      E[i]=Ee;
      pE[i]=p[i]*E[i];
      nE[i]=n[i]*E[i];
   }
   // slopes
   double dp[2*N+1], dn[2*N+1], dE[2*N+1], dpE[2*N+1], dnE[2*N+1];
   for (int i=1; i<2*N; i++) {
      dp[i] = (p[i+1]-p[i-1])/2;
      dn[i] = (n[i+1]-n[i-1])/2;
      dE[i] = (E[i+1]-E[i-1])/2;
      dpE[i] = (pE[i+1]-pE[i-1])/2;
      dnE[i] = (nE[i+1]-nE[i-1])/2;
   }
   dp[2*N] = dp[0] = 0;
   dn[2*N] = dn[0] = 0;
   dE[2*N] = dE[0] = 0;
   dpE[2*N] = dpE[0] = 0;
   dnE[2*N] = dnE[0] = 0;

   double d2p[2*N+1], d2n[2*N+1];
   for (int i=1; i<2*N; i++) {
      d2p[i] = (dp[i+1]-dp[i-1])/2;
      d2n[i] = (dn[i+1]-dn[i-1])/2;
   }
   d2p[2*N] = d2p[0] = 0;
   d2n[2*N] = d2n[0] = 0;

   double phi[2*N+1]; // weighting potential
   for (int i=0; i<2*N+1; i++) phi[i]=float(i)/2/N;

   //Transverse diffusion 
   double Re[2*N+1], Rh[2*N+1];
   for (int i=0; i<2*N+1; i++){
      Re[i]=18*n[i]*De/R/R;
      Rh[i]=18*p[i]*Dh/R/R;
   }

   // output
   TFile *output = new TFile("tube.root","recreate");
   TTree *t = new TTree("t","time slices");
   t->Branch("x",x,Form("x[%d]/D",2*N+1));
   t->Branch("p",p,Form("p[%d]/D",2*N+1));
   t->Branch("n",n,Form("n[%d]/D",2*N+1));
   t->Branch("E",E,Form("E[%d]/D",2*N+1));
//   t->Branch("pE",pE,Form("pE[%d]/D",2*N+1));
//   t->Branch("nE",nE,Form("nE[%d]/D",2*N+1));
//   t->Branch("dp",dp,Form("dp[%d]/D",2*N+1));
//   t->Branch("dn",dn,Form("dn[%d]/D",2*N+1));
//   t->Branch("dE",dE,Form("dE[%d]/D",2*N+1));
//   t->Branch("dpE",dpE,Form("dpE[%d]/D",2*N+1));
//   t->Branch("dnE",dnE,Form("dnE[%d]/D",2*N+1));
//   t->Branch("d2p",d2p,Form("d2p[%d]/D",2*N+1));
//   t->Branch("d2n",d2n,Form("d2n[%d]/D",2*N+1));
//   t->Branch("Re",Re,Form("Re[%d]/D",2*N+1));
//   t->Branch("Rh",Rh,Form("Rh[%d]/D",2*N+1));
//   t->Branch("dx",&dx,"dx/D");
//   t->Branch("dt",&dt,"dt/D");
   // save units into tree
   double CM3=cm3, CM2=cm2, CM=cm, UM=um, V=volt, SEC=s, NS=ns, NM=nm;
   t->Branch("cm3",&CM3,"cm3/D");
   t->Branch("cm2",&CM2,"cm2/D");
   t->Branch("V",&V,"V/D");
   t->Branch("cm",&CM,"cm/D");
   //   t->Branch("um",&UM,"um/D");
   t->Branch("nm",&NM,"nm/D");
   //   t->Branch("ns",&NS,"ns/D");
   t->Branch("s",&SEC,"s/D");
   //   t->Branch("eps",&epsilon,"eps/D");
   //   t->Branch("q",&Q,"q/D");

   int iStep=0, nSteps = 100;
   double time[1000]={0}, charge[1000]={0};

   // evolve
   if (argc>1) nSteps = atoi(argv[1]);
//  t->Branch("mue",&mu_e,"mue/D");
//  t->Branch("muh",&mu_h,"muh/D");
 
   while (iStep<nSteps) {
      t->Fill();
      // update charge
      for (int i=0; i<2*N+1; i++) {
         charge[iStep]+=(p[i]-n[i])*phi[i];
      }
      time[iStep]=iStep*dt;
      // update electron and hole distributions
      for (int i=0; i<2*N+1; i++) {
         double dn_dt = ge.Mu('e', n[i])*(dn[i]/dx*E[i]+n[i]*dE[i]/dx);
         double dp_dt = -ge.Mu('h', p[i])*(dp[i]/dx*E[i]+p[i]*dE[i]/dx);
         n[i]+=dn_dt*dt;
         p[i]+=dp_dt*dt;
      }
      // update electric field distribution
      for (int i=0; i<2*N+1; i++) {
         E[i]=0;
         for (int j=0; j<i; j++)
            E[i]+=(p[j]-n[j])*(1-(i-j)*dx/sqrt((i-j)*(i-j)*dx*dx+R*R));
         for (int j=i+1; j<2*N+1; j++) 
            E[i]+=(n[j]-p[j])*(1-(j-i)*dx/sqrt((j-i)*(j-i)*dx*dx+R*R));
         E[i]*=Q*dx/2/epsilon;
         if (iStep==(nSteps-1))
            cout<<"E["<<i<<"] = "<<E[i]<<endl;
         E[i]+=Ee;
         if (E[i]<0) E[i]=0; // large dt may over evolve things
      }
      
      // update pE[i]and nE[i]
      for (int i=0; i<2*N+1; i++){
         pE[i] = p[i]*E[i];
         nE[i] = n[i]*E[i];

      }
      // update slopes
      for (int i=1; i<2*N; i++) {
         dp[i] = (p[i+1]-p[i-1])/2;
         dn[i] = (n[i+1]-n[i-1])/2;
         dE[i] = (E[i+1]-E[i-1])/2;
         dpE[i] =(pE[i+1]-pE[i-1])/2;
         dnE[i] =(nE[i+1]-nE[i-1])/2;
      }
      dp[2*N] = dp[0] = 0;
      dn[2*N] = dn[0] = 0;
      dE[2*N] = dE[0] = 0;
      dpE[2*N] = dpE[0] = 0;
      dnE[2*N] = dnE[0] = 0;

      // update second derivative
      for (int i=1; i<2*N; i++) {
         d2p[i] = (dp[i+1]-dp[i-1])/2;
         d2n[i] = (dn[i+1]-dn[i-1])/2;
      }
      d2p[2*N] = d2p[0] = 0;
      d2n[2*N] = d2n[0] = 0;

      //update transverse diffusion
      for (int i=0; i<2*N+1; i++){
         Re[i]=(18*n[i]*De/R/R)*((18*De*iStep*dt/R/R)+1);
         Rh[i]=(18*p[i]*Dh/R/R)*((18*Dh*iStep*dt/R/R)+1);
      }

      iStep++;
   }
   t->Write("t",6);

  // TGraph *g = new TGraph(nSteps,time,charge);
  // g->Write("g");

   output->Close();

   return 0;
}
