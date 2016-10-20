//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/elastic.h"
#include "fiatlux/tools.h"
#include "fiatlux/integrator.h"
#include "fiatlux/settings.h"

#include <array>
using namespace std::placeholders;

namespace fiatlux
{
  //_________________________________________________________________________
  ElasticPhoton::ElasticPhoton():
    Integrator{},
    _mproton2(pow(input().get<double>("mproton"), 2)),
    _eps_base(input().get<double>("eps_base")),
    _eps_rel(input().get<double>("eps_rel")),
    _log_q2_max(log(input().get<double>("q2_max"))),
    _alpha_ref(input().get<double>("alpha_ref")),
    _elastic_param(input().get<int>("elastic_param"))
  {
  }

  //_________________________________________________________________________
  double ElasticPhoton::evaluatephoton(const double &x, const double &q2) const
  {
    double res = 0;
    if (x != 1.0)
      {
        const double eps_local = _eps_base * _eps_rel * pow(1.0-x, 4);
        const double Q2min = x*x*_mproton2/(1.0-x);
        res = integrate(log(Q2min), _log_q2_max, eps_local, x) * _alpha_ref / M_PI / 2.0;
      }
    return res;
  }

  //_________________________________________________________________________
  double ElasticPhoton::integrand(double const& lnQ2, const double &x) const
  {
    const double Q2 = exp(lnQ2);
    const auto   ge_gm = elastic_ge_gm(Q2);
    const double tau = Q2/(4*_mproton2);
    const double ge2 = pow(ge_gm[0], 2);
    const double gm2 = pow(ge_gm[1], 2);
    const double F2 = (ge2 + gm2*tau)/(1.0 + tau);
    const double FL = ge2/tau;

    double res = (2 - 2*x + x*x*(1.0 + 0.5/tau))*F2 - x*x * FL;

    return res;
  }

  //_________________________________________________________________________
  array<double,2> ElasticPhoton::elastic_ge_gm(double const& q2) const
  {
    array<double, 2> res;

    switch (_elastic_param) {
      case elastic_dipole:
        break;

      case elastic_A1_world_spline:
        break;

      default:
        info("ElasticPhoton::elastic_ge_gm", "Unrecognised elastic_param");
        throw std::runtime_error("Unrecognised elastic_param");
      }      

    return res;
  }
}
