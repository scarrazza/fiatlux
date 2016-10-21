//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/elastic.h"
#include "fiatlux/tools.h"
#include "fiatlux/integrator.h"
#include "fiatlux/settings.h"

#include <array>
#include <fstream>
using namespace std::placeholders;
using std::fstream;
using std::ios;
using std::getline;

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
    _mum_proton(input().get<double>("mum_proton")),
    _elastic_electric_rescale(input().get<double>("elastic_electric_rescale")),
    _elastic_magnetic_rescale(input().get<double>("elastic_magnetic_rescale")),
    _x(0),
    _elastic_param(input().get<int>("elastic_param"))
  {
    if (_elastic_param != elastic_dipole)
      {

        // read data from 1307.6227
        // interpolate latter.
        fstream f;
        switch (_elastic_param) {
          case elastic_A1_world_spline:
            f.open("data/elastic/World/CrossSectionsOnly/SplinesWithVariableKnots.dat", ios::in);
            break;

          case elastic_A1_world_pol_spline:
            f.open("data/elastic/World/CrossSectionsAndPolarized/SplinesWithVariableKnots.dat", ios::in);
            break;
          }

        if (f.fail())
          throw runtime_exception("ElasticPhoton::ElasticPhoton", "data file not found");

        // initial nodes for interpolation
        _fit.push_back({0,1,1});
        _fit_uperr.push_back({0,1,1});
        _fit_downerr.push_back({0,1,1});

        string tmp;
        double a, b, c, d, e, g;
        for (;;)
          {
            array<double, 3> cv, up, dn;

            // order: q2, ge, ge stat. error, upper model error, lower model error, similar for gm/mup and mup ge/gm
            f >> cv[0] >> cv[1] >> a >> b >> c >> cv[2] >> d >> e >> g; getline(f, tmp);
            if (f.eof()) break;
            up[0] = dn[0] = cv[0];
            up[1] = cv[1] + sqrt(a*a+b*b);
            up[2] = cv[2] + sqrt(d*d+e*e);
            dn[1] = cv[1] - sqrt(a*a+c*c);
            dn[2] = cv[1] - sqrt(d*d+g*g);
            _fit.push_back(cv);
            _fit_uperr.push_back(up);
            _fit_downerr.push_back(dn);
          }
        f.close();

        if (input().get<bool>("verbose"))
          info("ElasticPhoton::ElasticPhothon", "elastic fit loaded npoints " + std::to_string(_fit.size()));
      }
  }

  //_________________________________________________________________________
  double ElasticPhoton::evaluatephoton(const double &x, const double &q2)
  {
    double res = 0;
    if (x != 1.0)
      {
        _x = x;
        const double eps_local = _eps_base * _eps_rel * pow(1.0-x, 4);
        const double q2min = x*x*_mproton2/(1.0-x);
        res = integrate(log(q2min), _log_q2_max, eps_local) * _alpha_ref / M_PI / 2.0;
      }
    return res;
  }

  //_________________________________________________________________________
  double ElasticPhoton::integrand(double const& lnq2) const
  {
    const double q2 = exp(lnq2);
    const array<double,2> ge_gm = elastic_ge_gm(q2);
    const double tau = q2/(4*_mproton2);
    const double ge2 = pow(ge_gm[0], 2);
    const double gm2 = pow(ge_gm[1], 2);
    const double F2 = (ge2 + gm2*tau)/(1.0 + tau); // eq. (7a) 1607.04266
    const double FL = ge2/tau; // eq. (7b) 1607.04266

    // eq. (6) 1607.04266
    double res = (2 - 2*_x + _x*_x*(1.0 + 0.5/tau))*F2 - _x*_x * FL;

    return res;
  }

  //_________________________________________________________________________
  array<double,2> ElasticPhoton::elastic_ge_gm(double const& q2) const
  {
    array<double, 2> ge_gm;

    if (_elastic_param == elastic_dipole)
      {
        // use dipole model
        ge_gm[0] = elastic_dipole_factor(q2);
        ge_gm[1] = _mum_proton * ge_gm[0];
      }
    else if (_elastic_param == elastic_A1_world_spline ||
             _elastic_param == elastic_A1_world_pol_spline)
      {
        // perform linear interpolation for q2 = [0.005, 10] GeV2
        // returns error if q2 < 0.05 GeV2
        // returns dipole model in extrapolation region q2 > 10. GeV2
        int iQ2 = 1000;
        if (q2 <= 10) iQ2 = floor((q2+_fit[1][0])/(_fit[2][0]-_fit[1][0]));

        if (iQ2 < 0)
          throw runtime_exception("ElasticPhoton::elastic_ge_gm", "iQ2 < 0");
        else if (iQ2 >= 1000)
          {
            const double q2max = _fit.back()[0];
            const double frac = elastic_dipole_factor(q2)/elastic_dipole_factor(q2max);
            ge_gm[0] = _fit.back()[1] * frac;
            ge_gm[1] = _fit.back()[2] * frac;
          }
        else
          {
            const double weight = (_fit[iQ2+1][0] - q2)/(_fit[iQ2+1][0]-_fit[iQ2][0]);
            ge_gm[0] = _fit[iQ2][1] * weight + _fit[iQ2+1][1] * (1.0-weight);
            ge_gm[1] = _fit[iQ2][2] * weight + _fit[iQ2+1][2] * (1.0-weight);
          }

        ge_gm[1] *= _mum_proton;
      }
    else
      throw runtime_exception("ElasticPhoton::elastic_ge_gm", "Unrecognised elastic_param");

    ge_gm[0] *= _elastic_electric_rescale;
    ge_gm[1] *= _elastic_magnetic_rescale;
    return ge_gm;
  }

  //_________________________________________________________________________
  double ElasticPhoton::elastic_dipole_factor(const double &q2) const
  {
    // eq. (30) 1307.6227
    return 1.0/pow(1.0+q2/0.71, 2);
  }
}
