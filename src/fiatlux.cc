//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/fiatlux.h"
#include "fiatlux/tools.h"
#include "fiatlux/settings.h"
#include <exception>

namespace fiatlux
{
  //_________________________________________________________________________
  FiatLux::FiatLux(xfxq const& f):
    _xfxq(f),   
    _elastic{}
  {
    string me = "FiatLux::FiatLux";

    if (s().qed_running)
      {
        info(me, "QED running activated");
        throw std::runtime_error("not implemented");
      }
  }

  //_________________________________________________________________________
  luxqed FiatLux::evaluate(const double &x, const double &q2) const
  {
    luxqed e;
    e.elastic = _elastic.evaluate(x,q2);
    return e;
  }
}
