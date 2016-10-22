//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/proton.h"
#include "fiatlux/settings.h"

namespace fiatlux
{
  //_________________________________________________________________________
  ProtonStructure::ProtonStructure():
    _mproton2(pow(input().get<double>("mproton"), 2)),
    _mum_proton(input().get<double>("mum_proton")),
    _eps_base(input().get<double>("eps_base")),
    _eps_rel(input().get<double>("eps_rel")),
    _log_q2_max(log(input().get<double>("q2_max"))),
    _alpha_ref(input().get<double>("alpha_ref"))
  {
  }

  //_________________________________________________________________________
  StrucFunc ProtonStructure::compute_proton_structure(const double &x, const double &q2) const
  {
    StrucFunc sf;



    return sf;
  }
}
