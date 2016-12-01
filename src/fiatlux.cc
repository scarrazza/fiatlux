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
  FiatLux::FiatLux(const string &filename)
  {
    input().load(filename);
    _elastic = unique_ptr<ElasticPhoton>(new ElasticPhoton{});
    _inelastic = unique_ptr<InelasticPhoton>(new InelasticPhoton{});
    _msbar = unique_ptr<MSbarPhoton>(new MSbarPhoton{});
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
  }

  //_________________________________________________________________________
  void FiatLux::plug_alphaqed(alpha_running const& a) const
  {
    _elastic->set_alpha_running(a);
    _inelastic->set_alpha_running(a);
    _msbar->set_alpha_running(a);
  }

  //_________________________________________________________________________
  luxqed FiatLux::evaluatephoton(const double &x, const double &mu2) const
  {
    luxqed e;
    e.elastic = _elastic->evaluatephoton(x, mu2);
    e.inelastic_pf = _inelastic->evaluatephoton(x, mu2);
    e.msbar_pf = _msbar->evaluatephoton(x, mu2);
    e.total = e.elastic + e.inelastic_pf + e.msbar_pf;
    return e;
  }
}
