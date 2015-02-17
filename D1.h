#ifndef D1_H
#define D1_H

#include <vector>

#include <TObject.h>

namespace PLASMA {
  class D1; // 1 dimensional distribution
}

class D1: public TObject
{
   public:
      const Int_t Nmax = 1000; // maximal number of grids in one direction

   protected:
      std::vector<Double_t> vec;

   public:
      D1(Int_t size=2*Nmax+1);
      virtual ~D1() {};

      Int_t N() { return (vec.size()-1)/2; }
      void SetN(Int_t n=1000) { vec.resize(2*n+1); }

      /**
       * Get slope around i.
       */
      Double_t Slope(Int_t i);

      /**
       * Get the i-th value;
       */
      Double_t At(Int_t i);
      Double_t operator[](Int_t i) { At(i); }

      ClassDef(D1,1);
};

#endif
