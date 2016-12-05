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
  }

  //_________________________________________________________________________
  void FiatLux::plug_alphaqed(alpha_running const& a) const
  {
    _elastic->set_alpha_running(a);
    _inelastic->set_alpha_running(a);
    _msbar->set_alpha_running(a);
  }

  //_________________________________________________________________________
  void FiatLux::plug_f2_fl(ext_sf const& f2, ext_sf const& fl) const
  {
    _inelastic->set_sf(f2, fl);
    _msbar->set_sf(f2, fl);
  }

  //_________________________________________________________________________
  void FiatLux::insert_inel_split(vector<double> const& qthresholds) const
  {
    for (auto const& v: qthresholds)
      _inelastic->insert_inel_split(pow(v, 2));
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
