#ifndef D1_H
#define D1_H

#include <vector>

#include <TObject.h>

namespace PLASMA {
  class D1; // 1 dimensional distribution
  const Int_t Nmax = 1000; // maximal number of grids in one direction
}

class D1: public TObject
{
   protected:
      std::vector<Double_t> vec;
      Double_t dx; // interval

   public:
      D1(Int_t size=2*Nmax+1, Double_t interval=1): dx(interval) { vec.resize(size); }
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

      void SetInterval(Double_t interval) { dx=interval; }
      Double_t Interval() { return dx; }

      ClassDef(D1,1);
};

#endif
