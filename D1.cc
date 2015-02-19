#include "D1.h"

using namespace PLASMA;

D1::D1(Int_t size): TObject()
{
   vec.resize(size);
}

//______________________________________________________________________________
//

Double_t D1::Slope(Int_t i)
{
   if (i==N()) i--;
   return (At(i+1)-At(i));
}

//______________________________________________________________________________
//

Double_t D1::At(Int_t i)
{
   if (i<-N() || i>N()) {
      Warning("At", "i=%d is out of [-%d, %d]! Return 0", i, N(), N());
      return 0;
   }
   return vec.at(i+N());
}
