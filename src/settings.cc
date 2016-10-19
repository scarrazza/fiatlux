//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/settings.h"

namespace fiatlux
{
  //_________________________________________________________________________
  Settings& s()
  {
    static Settings settings{};
    return settings;
  }

  //_________________________________________________________________________
  Settings::Settings():
    qed_running(false),
    q2_max(1E9),
    eps_base(1E-5),
    eps_rel(1E-1),
    mproton(0.938272046),
    alpha_ref(1.0/137.035999074)
  {
  }
}
