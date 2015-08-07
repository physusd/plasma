#include <iostream>
using namespace std;
// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TChain.h>
#include <TGraph.h>
using namespace TMath;
// UNIC
#include <UNIC/Units.h>
#include <UNIC/Constants.h>
using namespace UNIC;
// MAD
#include <MAD/GeCrystal.h>
using namespace MAD;
const int N = 2000;
// df at i
double d(double *f, int i)
{
   if (i<1 || i>2*N-1) return 0;
   return (f[i+1]-f[i-1])/2; // 2 points
   //if (i<2 || i>2*N-2) return 0;
   //return (-f[i+2]+8*f[i+1]-8*f[i-1]+f[i-2])/12; // 4
   //if (i<4 || i>2*N-4) return 0;
   //return (-f[i+4]-8*f[i+2]+128*f[i+1]-128*f[i-1]+8*f[i-2]+f[i-4])/180; // 6
}
// thin plasma tube
int main(int argc, char** argv)
{
   double dt=5e-7*ns; // large dt causes E to decrease too fast
   double dx=0.2*nm; // small dx may cause artificial oscillation
   double Ee=500*volt/cm;
   double R=1.97e-7*cm, mean=0, sigma=R/3, height=84./(Pi()*R*R);
   bool norm;

   GeCrystal ge;
   double epsilon=ge.Epsilon()*epsilon0;
   double Q=Abs(electron_charge);
   //double De=0.046*cm2/s, Dh=0.046*cm2/s;

   // initialize arrays
   // p and n are number densities
   double x[2*N+1], p[2*N+1], n[2*N+1], E[2*N+1], pE[2*N+1],nE[2*N+1], T[10000], J[10000];
   if (argc==3) {
      TChain *ti = new TChain("t");
      ti->Add(argv[2]);
      ti->SetBranchAddress("x",x);
      ti->SetBranchAddress("p",p);
      ti->SetBranchAddress("n",n);
      ti->SetBranchAddress("E",E);
      ti->GetEntry(ti->GetEntries()-1);

   } else {
      for (int i=0; i<2*N+1; i++) {
         x[i]=(i-N)*dx;
         p[i]=n[i]=height*Gaus(x[i],mean=0,sigma,norm=kTRUE);
         if (p[i]<1e-5/cm3) p[i]=0;
         if (n[i]<1e-5/cm3) n[i]=0;
         E[i]=Ee;
      }
   }
   for (int i=0; i<2*N+1; i++) {
      pE[i]=p[i]*E[i];
      nE[i]=n[i]*E[i];
   }

   for (int i=0; i<1000; i++){
      J[i] = 0;
      T[i] = 0;
   }

   // output
   TFile *output = new TFile("tube.root","recreate");
   TTree *t = new TTree("t","time slices");
   t->Branch("x",x,Form("x[%d]/D",2*N+1));
   t->Branch("p",p,Form("p[%d]/D",2*N+1));
   t->Branch("n",n,Form("n[%d]/D",2*N+1));
   t->Branch("E",E,Form("E[%d]/D",2*N+1));
   // save units into tree
   double CM3=cm3, CM2=cm2, CM=cm, UM=um, V=volt, SEC=s, NS=ns, NM=nm;
   t->Branch("cm3",&CM3,"cm3/D");
   t->Branch("cm2",&CM2,"cm2/D");
   t->Branch("V",&V,"V/D");
   t->Branch("cm",&CM,"cm/D");
   t->Branch("um",&UM,"um/D");
   t->Branch("nm",&NM,"nm/D");
   t->Branch("ns",&NS,"ns/D");
   t->Branch("s",&SEC,"s/D");

   int iStep=0, nSteps = 100;

   // evolve
   if (argc>1) nSteps = atoi(argv[1]);
   while (iStep<nSteps) {
      if ((iStep+1)%100==0) cout<<iStep+1<<" steps passed"<<endl;
      t->Fill();
      // update electron and hole distributions
      for (int i=0; i<2*N+1; i++) {
         //double dn_dt = ge.Mu('e', n[i])*(dn[i]/dx*E[i]+n[i]*dE[i]/dx);
         //double dp_dt = -ge.Mu('h', p[i])*(dp[i]/dx*E[i]+p[i]*dE[i]/dx);
         double dn_dt = ge.Mu('e', n[i])*d(nE,i)/dx;
         double dp_dt = -ge.Mu('h', p[i])*d(pE,i)/dx;
         if (Abs(dn_dt*dt)<1e-10/cm3) dn_dt=0;
         if (Abs(dp_dt*dt)<1e-10/cm3) dp_dt=0;
         n[i]+=dn_dt*dt;
         p[i]+=dp_dt*dt;
         if (n[i]<0) n[i]=0;
         if (p[i]<0) p[i]=0;
      }

      if (iStep<=1000&&(iStep+1)%20==0){
         double aa = 0; //initialization 
         for (int i=N; i<2*N+1; i++)
            aa += (p[i]-n[i]);
         J[(iStep+1)/20-1] = aa*Q*dx/(iStep+1)/dt/(coulomb/s/cm2);
         T[(iStep+1)/20-1] = (iStep+1);
         cout<<"test: "<<(iStep+1)/20-1<<",J: "<<J[(iStep+1)/20-1]<<", T: "<<T[(iStep+1)/20-1]<<endl;
      }

 if (iStep>1000 && iStep<=2000 && (iStep+1)%100==0){
          double aa = 0; //initialization 
          for (int i=N; i<2*N+1; i++)
             aa += (p[i]-n[i]);
          J[(iStep+1)/100+39] = aa*Q*dx/(iStep+1)/dt/(coulomb/s/cm2);
          T[(iStep+1)/100+39] = (iStep+1);
          cout<<"test: "<<(iStep+1)/100+39<<",J: "<<J[(iStep+1)/100+39]<<", T: "<<T[(iStep+1)/100+39]<<endl;
                }
          //

  if (iStep>2000&&(iStep+1)%1000==0){
            double aa = 0; //initialization 
            for (int i=N; i<2*N+1; i++)
               aa += (p[i]-n[i]);
            J[(iStep+1)/1000+57] = aa*Q*dx/(iStep+1)/dt/(coulomb/s/cm2);
            T[(iStep+1)/1000+57] = (iStep+1);
            cout<<"test: "<<(iStep+1)/1000+57<<",J: "<<J[(iStep+1)/1000+57]<<", T: "<<T[(iStep+1)/1000+57]<<endl;
                  }

      // update electric field distribution
      for (int i=0; i<2*N+1; i++) {
         E[i]=0;
         for (int j=0; j<i; j++)
            E[i]+=(p[j]-n[j])*(1-(i-j)*dx/sqrt((i-j)*(i-j)*dx*dx+R*R));
         for (int j=i+1; j<2*N+1; j++) 
            E[i]+=(n[j]-p[j])*(1-(j-i)*dx/sqrt((j-i)*(j-i)*dx*dx+R*R));
         E[i]*=Q*dx/2/epsilon;
         E[i]+=Ee;
         if (E[i]<1e-10*volt/cm) E[i]=0; // large dt may over evolve things
      }

      // update pE[i]and nE[i]
      for (int i=0; i<2*N+1; i++){
         pE[i] = p[i]*E[i];
         nE[i] = n[i]*E[i];
      }

      iStep++;
   }
   t->Write("t",6);

   TGraph *g = new TGraph (nSteps, T, J);
   g->Write("g");
   output->Close();

   return 0;
}
