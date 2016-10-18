//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <functional>
using std::function;

namespace fiatlux
{
  using partons = function<double(int,double,double)>;

  class FiatLux
  {
  public:
    FiatLux(partons const&f);

  private:
    partons _xfxQ;
  };
}
