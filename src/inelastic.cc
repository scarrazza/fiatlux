//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/inelastic.h"
#include "fiatlux/tools.h"
#include "fiatlux/integrator.h"
#include "fiatlux/settings.h"
#include <algorithm>
#include <cmath>

using std::min;
using std::max;
using std::cout;
using std::endl;

namespace fiatlux
{
  //_________________________________________________________________________
  InelasticPhoton::InelasticPhoton(const unique_ptr<ProtonStructure> &proton):
    Integrator{},
    _proton(proton),
    _q2min_inel_override(input().get<double>("q2min_inel_override")),
    _q2max_inel_override(input().get<double>("q2max_inel_override")),
    _use_mu2_as_upper_limit(input().get<bool>("use_mu2_as_upper_limit")),
    _inelastic_q2(InelasticQ2{proton})
  {
  }

  //_________________________________________________________________________
  double InelasticPhoton::evaluatephoton(const double &x, const double& mu2) const
  {
    const double eps = _proton->_eps_base * pow(1.0-x, 4);
    double res = integrate(0, log(1.0/x), eps, {x, mu2, eps}) * _proton->_alpha_ref / M_PI / 2.0;

    if (_proton->_qed_running)
      res *= _proton->_alpha_ref/_proton->_alpha_running(sqrt(mu2));

    return res;    
  }

  //_________________________________________________________________________
  double InelasticPhoton::integrand(double const& ln1oz, const vector<double> &e) const
  {
    const double z = exp(-ln1oz);
    const double x = e[0];
    const double q2 = e[1];
    const double eps = e[2];

    double q2min = x*x * _proton->_mproton2 / (1.0 - z);
    double q2max = q2;
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

    if (_q2_inel_split.size() < 2)
      throw runtime_exception("InelasticPhoton::integrand", "inconsistent n_inel_split");

    double res = 0;
    for (size_t i = 0; i < _q2_inel_split.size()-1; i++)
      {
        const double q2lo = max(q2min, _q2_inel_split[i]);
        const double q2hi = min(q2max, _q2_inel_split[i+1]);
        if (q2lo < q2hi)
          res += _inelastic_q2.integrate(log(q2lo), log(q2hi), eps*_proton->_eps_rel, {x, z});
      }

    return res*mult;
  }

  //_________________________________________________________________________
  void InelasticPhoton::insert_inel_split(const double &q2)
  {
    if (_q2_inel_split.size() > 20)
      throw runtime_exception("InelasticPhoton::insert_inel_split", "no room to insert further inel split");
    _q2_inel_split.push_back(q2);
    std::sort(_q2_inel_split.begin(), _q2_inel_split.end());

    if (input().get<bool>("verbose"))
      {
        stringstream ss("");
        for (const auto &i: _q2_inel_split)
          ss << i <<  " ";
        info("InelasticPhoton::insert_inel_split", ss.str());
      }
  }

  //_________________________________________________________________________
  double InelasticQ2::integrand(const double &lnq2, vector<double> const& e) const
  {
    const double x = e[0];
    const double z = e[1];
    const double q = exp(0.5*lnq2);
    const double q2 = q*q;
    const auto kin = _proton->compute_proton_structure(x/z, q2);

    double res = -z*z * kin.FL + (2 - 2*z + z*z + (2 * x*x * _proton->_mproton2)/q2)* kin.F2;

    if (_proton->_qed_running)
      res *= pow(_proton->_alpha_running(q)/_proton->_alpha_ref, 2);

    return res;
  }
}
