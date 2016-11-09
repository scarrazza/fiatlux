//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/msbar.h"
#include "fiatlux/tools.h"
#include "fiatlux/integrator.h"
#include "fiatlux/settings.h"

namespace fiatlux
{
  //_________________________________________________________________________
  MSbarPhoton::MSbarPhoton():
    Integrator{}, ProtonStructure{},
    _use_mu2_as_upper_limit(input().get<bool>("use_mu2_as_upper_limit"))
  {
  }

  //_________________________________________________________________________
  double MSbarPhoton::evaluatephoton(const double &x, const double &q2) const
  {
    double res = 0;
    const double eps_local = _eps_base * _eps_rel * pow(1.0-x, 4);
    res = integrate(0, log(1.0/x), eps_local, {x,q2}) * _alpha_ref / M_PI / 2.0;
    return res;
  }

  //_________________________________________________________________________
  double MSbarPhoton::integrand(double const& ln1oz, const vector<double> &e) const
  {
    const double x = e[0];
    const double q2 = e[1];
    const double z = exp(-ln1oz);
    const auto kin = compute_proton_structure(x/z, q2);

    double res = -pow(z,2)*kin.F2;
    if (_use_mu2_as_upper_limit)
      res += log(1.0/(1.0-z))*(1+pow(1.0-z,2))*kin.F2;

    return res;
  }
}
