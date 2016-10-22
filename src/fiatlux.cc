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
    _inelastic->insert_inel_split(-2);

    const double lhtrq2 = input().get<double>("lhapdf_transition_q2");
    if (lhtrq2 > 0)
      _inelastic->insert_inel_split(lhtrq2);
    _inelastic->insert_inel_split(1E200);

    if (input().get<bool>("verbose"))
      {
        info("InelasticPhoton::insert_inel_split","q2_inel_split:");
        stringstream ss("");
        for (const auto &i: _inelastic->get_q2_inel_split())
          ss << i <<  " ";
        info("InelasticPhoton::insert_inel_split", ss.str());
      }

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
