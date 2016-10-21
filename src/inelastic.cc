//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/inelastic.h"
#include "fiatlux/tools.h"
#include "fiatlux/integrator.h"
#include "fiatlux/settings.h"

using std::min;
using std::max;

namespace fiatlux
{
  //_________________________________________________________________________
  InelasticPhoton::InelasticPhoton():
    Integrator{},
    _eps_base(input().get<double>("eps_base")),
    _eps_rel(input().get<double>("eps_rel")),
    _alpha_ref(input().get<double>("alpha_ref")),
    _x(0),
    _q2(0),
    _q2min_inel_override(input().get<double>("q2min_inel_override")),
    _q2max_inel_override(input().get<double>("q2max_inel_override")),
    _use_mu2_as_upper_limit(input().get<bool>("use_mu2_as_upper_limit")),
    _inelastic_q2(InelasticQ2{})
  {
  }

  //_________________________________________________________________________
  double InelasticPhoton::evaluatephoton(const double &x, const double& q2)
  {
    _x = x; _q2 = q2;
    const double eps_local = _eps_base * _eps_rel * pow(1.0-x, 4);
    double res = integrate(0, log(1.0/x), eps_local) * _alpha_ref / M_PI / 2.0;
    return res;
  }

  //_________________________________________________________________________
  double InelasticPhoton::integrand(double const& ln1oz) const
  {
    double z = exp(-ln1oz);
    double q2min = _x*_x * _mproton2 / (1 - z);

    double q2max = _q2;
    if (!_use_mu2_as_upper_limit)
      q2max /= 1.0 - z;

    q2min = min(max(q2min,_q2min_inel_override),_q2max_inel_override);
    q2max = min(max(q2max,_q2min_inel_override),_q2max_inel_override);

    double mult = 1;
    if (q2min > q2max)
      {
        swap(&q2min, &q2max);
        mult = -1;
      }

    double res = 0;

    return res*mult;
  }

  //_________________________________________________________________________
  double InelasticQ2::integrand(const double &lnq2) const
  {
    return 0;
  }
}
