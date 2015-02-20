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

Double_t& D1::At(Int_t i)
{
   if (i<-N() || i>N()) Error("At", "i=%d is out of [-%d, %d]!", i, N(), N());
   return vec.at(i+N());
}
