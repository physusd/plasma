#include "D1.h"

using namespace PLASMA;

D1::Slope(Int_t i)
{
   if (i==N()) i--;
   return (At[i+1]-At[i])/dx;
}

//______________________________________________________________________________
//

D1:At(Int_t i)
{
   if (i<-N() || i>N()) {
      Warning("i=%d is out of [-%d, %d]! Return 0", i, N(), N());
      return 0;
   }
   return vec.at(i+N());
}
