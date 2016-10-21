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
  FiatLux::FiatLux(const string &filename, xfxq const& f):
    _xfxq(f)
  {
    input().load(filename);
    _elastic = unique_ptr<ElasticPhoton>(new ElasticPhoton{});
    _inelastic = unique_ptr<InelasticPhoton>(new InelasticPhoton{});

    if (input().get<bool>("qed_running"))
      throw runtime_exception("FiatLux:FiatLux", "QED running activated");
  }

  //_________________________________________________________________________
  luxqed FiatLux::evaluatephoton(const double &x, const double &q2) const
  {
    luxqed e;
    e.elastic = _elastic->evaluatephoton(x, q2);
    e.inelastic_pf = _inelastic->evaluatephoton(x, q2);
    e.msbar_pf = 0;
    e.total = e.elastic + e.inelastic_pf + e.msbar_pf;
    return e;
  }
}
