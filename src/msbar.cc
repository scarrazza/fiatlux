//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/msbar.h"
#include "fiatlux/tools.h"
#include "fiatlux/integrator.h"
#include "fiatlux/settings.h"
#include <cmath>

namespace fiatlux
{
  //_________________________________________________________________________
  MSbarPhoton::MSbarPhoton(const unique_ptr<ProtonStructure> &proton):
    Integrator{},
    _proton(proton),
    _use_mu2_as_upper_limit(input().get<bool>("use_mu2_as_upper_limit"))
  {
  }

  //_________________________________________________________________________
  double MSbarPhoton::evaluatephoton(const double &x, const double &mu2) const
  {
    double res = 0;
    const double eps_local = _proton->_eps_base * _proton->_eps_rel * pow(1.0-x, 4);
    res = integrate(0, log(1.0/x), eps_local, {x,mu2}) * _proton->_alpha_ref / M_PI / 2.0;

    if (_proton->_qed_running)
      res *= _proton->_alpha_running(sqrt(mu2))/_proton->_alpha_ref;

    return res;
  }

  //_________________________________________________________________________
  double MSbarPhoton::integrand(double const& ln1oz, const vector<double> &e) const
  {
    const double x = e[0];
    const double q2 = e[1];
    const double z = exp(-ln1oz);
    const auto kin = _proton->compute_proton_structure(x/z, q2, true);

    // following LUX17 we use F2 at LO instead of NLO/NNLO.
    double res = -pow(z,2)*kin.F2LO;
    if (_use_mu2_as_upper_limit)
      res += log(1.0/(1.0-z))*(1+pow(1.0-z,2))*kin.F2LO;

    return res;
  }
}
