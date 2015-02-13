void main(int argc, int** argv)
{
   double dt=0.1*ns;
   double dx=0.1*nm;
   double Ee=1000*V/cm;

   double t=0;
   for (int i=-imin; i<imax; i++) {
      p[i]=-n[i]=gaus(i);
   }
   while (t<100*ns) {
      for (int i=-imin; i<imax; i++) {
         E[i]=0;
         for (int j=-imin; j<i; j++) {
            E[i]+=q/dielec*(p[j]-n[j])*dx;
         }
         E[i]+=Ee;

         pSlope_x[i]=(p[i+1]-p[i])/dx;
         nSlope_x[i]=(n[i+1]-n[i])/dx;

         pSlope_t[i] - nSlope_t[i] = - mu*E[i]*(pSlope_x[i]-nSlope_x[i]) - q*mu/dielec*(p[i]-n[i])*(p[i]-n[i]);
      }

      t=t+dt;
   }
}
double gaus(int)
{
   return 1;
}
