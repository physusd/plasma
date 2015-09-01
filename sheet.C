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

// parameters
double dx;
double dt=5e-7*ns; // large dt causes E to decrease too fast
double Ee=1000*volt/cm;
double R=1.97e-7*cm, mean=0, sigma=R/3, height=84./(Pi()*R*R);

const int Nmidd = 50; // # of points in (0,5*sigma)
const int Nwing = 95; // # of points in (5*sigma, 100*sigma)
const int Ntotal = 2*(Nmidd + Nwing) + 1;

GeCrystal ge;
double epsilon=ge.Epsilon()*epsilon0;
double Q=Abs(electron_charge);

// df at i
double d(double *f, int i)
{
   if (i<1 || i>Ntotal-1) return 0;
   return (f[i+1]-f[i-1])/2; // 2 points
   //if (i<2 || i>Ntotal-2) return 0;
   //return (-f[i+2]+8*f[i+1]-8*f[i-1]+f[i-2])/12; // 4
   //if (i<4 || i>Ntotal-4) return 0;
   //return (-f[i+4]-8*f[i+2]+128*f[i+1]-128*f[i-1]+8*f[i-2]+f[i-4])/180; // 6
}

// thin flat infinite plasma sheet
int main(int argc, char** argv)
{
   // initialize arrays
   double x[Ntotal], E[Ntotal], p[Ntotal], n[Ntotal], pE[Ntotal], nE[Ntotal];
   for (int i=0; i<=Nwing; i++) { // wings
      dx = sigma;
      x[i]=-100*sigma + i*dx;
      x[Ntotal-1-i]= 100*sigma - i*dx;
      p[i]=n[i]=height*Gaus(x[i],mean=0,sigma,kTRUE);
      p[Ntotal-1-i]=height*Gaus(x[Ntotal-1-i],mean=0,sigma,kTRUE);
      n[Ntotal-1-i]=height*Gaus(x[Ntotal-1-i],mean=0,sigma,kTRUE);
      E[i]=E[Ntotal-1-i]=Ee;
      pE[i]=p[i]*E[i];
      nE[i]=n[i]*E[i];
      pE[Ntotal-1-i]=p[Ntotal-1-i]*E[Ntotal-1-i];
      nE[Ntotal-1-i]=n[Ntotal-1-i]*E[Ntotal-1-i];
   }
   for (int i=Nwing+1; i<Ntotal-Nwing; i++) { // middle
      dx = sigma/10;
      x[i]=x[Nwing] + (i-Nwing)*dx;
      p[i]=n[i]=height*Gaus(x[i],mean=0,sigma,kTRUE);
      E[i]=Ee;
      pE[i]=p[i]*E[i];
      nE[i]=n[i]*E[i];
   }

   // output
   TFile *output = new TFile("sheet.root","recreate");
   TTree *t = new TTree("t","time slices");
   t->Branch("x",x,Form("x[%d]/D",Ntotal));
   t->Branch("p",p,Form("p[%d]/D",Ntotal));
   t->Branch("n",n,Form("n[%d]/D",Ntotal));
   t->Branch("E",E,Form("E[%d]/D",Ntotal));
   t->Branch("dt",&dt,"dt/D");
   // save units into tree
   double V=volt, nanosec=ns, nanometer=nm;
   t->Branch("V",&V,"V/D");
   t->Branch("nm",&nanometer,"nm/D");
   t->Branch("ns",&nanosec,"ns/D");

   // evolve
   int iStep=0, nSteps = 100;
   if (argc>1) nSteps = atoi(argv[1]);
   while (iStep<nSteps) {
      t->Fill();
      for (int i=0; i<Ntotal; i++) {
         if (i==Ntotal-1) dx = x[i]-x[i-1];
         else if (i==Nwing || i==Ntotal-1-Nwing) dx = (x[i+1]-x[i-1])/2;
         else dx = (x[i+1]-x[i]);
         double dn_dt = ge.Mu('e', n[i])*d(nE,i)/dx;
         double dp_dt =-ge.Mu('h', p[i])*d(pE,i)/dx;
         double dE_dt = E[i]*(-ge.Mu('e',n[i])*n[i]-ge.Mu('h',p[i])*p[i])/epsilon/Q;
         n[i]+=dn_dt*dt;
         p[i]+=dp_dt*dt;
         E[i]+=dE_dt*dt;
      }
      iStep++;
   }
   t->Write("t",6);

   output->Close();

   return 0;
}
