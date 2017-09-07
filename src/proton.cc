//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/proton.h"
#include "fiatlux/settings.h"
#include "fiatlux/tools.h"
#include <cmath>

extern "C" {
  void gd_fit_11_(double const&, double const&,
                  double const&, const char*, double &, double &, double &, double&);
  void f2sf_(int const&, double const&, double const&, double &);
}

namespace fiatlux
{
  //_________________________________________________________________________
  ProtonStructure::ProtonStructure():
    _qed_running(input().get<bool>("qed_running")),
    _alpha_ref(1/137.035999074), // @ electron mass = 0.000510998946
    _mproton2(pow(input().get<double>("mproton"), 2)),
    _mum_proton(input().get<double>("mum_proton")),
    _eps_base(input().get<double>("eps_base")),
    _eps_rel(input().get<double>("eps_rel")),
    _log_q2_max(log(input().get<double>("q2_max"))),
    _rescale_r(input().get<double>("rescale_r")),
    _rescale_r_twist4(input().get<double>("rescale_r_twist4")),
    _lhapdf_transition_q2(input().get<double>("lhapdf_transition_q2")),
    _allm_limits(input().get<double>("allm_limits")),
    _rescale_non_resonance(input().get<double>("rescale_non_resonance")),
    _rescale_resonance(input().get<double>("rescale_resonance")),
    _HAC_loW2(3.0),
    _HAC_hiW2(4.0),
    _inelastic_param(input().get_inelastic_param())
  {
  }

  //_________________________________________________________________________
  double ProtonStructure::w2_from_xq2(double const& x, double const& q2) const
  {
    return (1.0-x)/x*q2+_mproton2;
  }

  //_________________________________________________________________________
  double ProtonStructure::x_from_w2q2(double const& w2, double const& q2) const
  {
    return 1.0/(1.0+(w2-_mproton2)/q2);
  }

  //_________________________________________________________________________
  void ProtonStructure::rescaler(double &rvalue, double const& q2) const
  {
    double rescale;
    if ((q2 > _lhapdf_transition_q2 && _lhapdf_transition_q2 >= 0) ||
        _inelastic_param == inelastic_LHAPDF_Hermes_ALLM_CLAS)
       rescale = (1 + _rescale_r_twist4/q2);
    else
       rescale = _rescale_r;
    rvalue *= rescale;
    // R must be positive (because sigmaL >= 0)
    if (rvalue < 0) rvalue = 0;
  }

  //_________________________________________________________________________
  StrucFunc ProtonStructure::compute_proton_structure(const double &x, const double &q2, bool f2lo_only) const
  {
    StrucFunc sf;
    sf.x = x;
    sf.q2 = q2;
    sf.w2 = w2_from_xq2(x, q2);

    switch (_inelastic_param) {
      case inelastic_Hermes_ALLM_CLAS:
        Hermes_ALLM_CLAS(sf);
        break;
      case inelastic_LHAPDF_Hermes_ALLM_CLAS:
        LHAPDF_Hermes_ALLM_CLAS(sf, f2lo_only);
        break;
      default:
        throw runtime_exception("ProtonStructure::compute_proton_structure", "Unrecognised inelastic_param");
      }
    return sf;
  }

  //_________________________________________________________________________
  void ProtonStructure::Hermes_ALLM_CLAS(StrucFunc & sf) const
  {
    double rvalue = R1998_protected_lowQ2(sf.x, sf.q2, sf.w2);
    rescaler(rvalue, sf.q2);

    double sigmaTL = 0;
    if (sf.w2 > _HAC_hiW2)
      {
        // get Hermes ALLM
        sigmaTL = allm_sigma(sf.x, sf.q2, sf.w2);
        sigmaTL *= _rescale_non_resonance;
       }
    else if (sf.w2 < _HAC_loW2)
      {
        // get CLAS
        double F2;
        f2sf_(1,sf.q2,sf.x,F2);
        F2 *= _rescale_resonance;
        sigmaTL = sigmaTL_from_F2(sf.x, sf.q2, F2);
      }
    else
      {        
        double F2_hiW, F2_loW;
        // get ALLM                

        sigmaTL = allm_sigma(sf.x, sf.q2, sf.w2);
        F2_hiW = F2_from_sigmaTL(sf.x, sf.q2, sigmaTL);
        F2_hiW *= _rescale_non_resonance;

        // get CLAS
        f2sf_(1, sf.q2, sf.x, F2_loW);
        F2_loW *= _rescale_resonance;

        // then interpolate between them
        const double u = (sf.w2 - _HAC_loW2)/(_HAC_hiW2 - _HAC_loW2);

        // a function that has zero derivative at u = 0, 1
        // and is respectively 0 and 1 there
        const double wgt_hiW = 2*pow(u,2) - pow(u,4);
        const double wgt_loW = 1.0 - wgt_hiW;
        const double F2 = F2_hiW * wgt_hiW + F2_loW * wgt_loW;

        sigmaTL = sigmaTL_from_F2(sf.x, sf.q2, F2);
      }

    // fill F2 and FL
    sf.F2 = sf.F2LO = F2_from_sigmaTL(sf.x, sf.q2, sigmaTL);
    sf.FL = FL_from_F2R(sf.x, sf.q2, sf.F2, rvalue);
  }

  //_________________________________________________________________________
  void ProtonStructure::LHAPDF_Hermes_ALLM_CLAS(StrucFunc & sf, bool f2lo_only) const
  {
    if (sf.w2 > _HAC_hiW2 && sf.q2 > _lhapdf_transition_q2)
      {
        double rvalue = 1;
        rescaler(rvalue, sf.q2);
        const double q = sqrt(sf.q2);

        // call external functions
        if (f2lo_only)
          sf.F2LO = _f2lo(sf.x, q)*_rescale_non_resonance; // following LUX17
        else
          {
            sf.F2 = _f2(sf.x, q)*_rescale_non_resonance;
            sf.FL = rvalue*_fl(sf.x, q)*_rescale_non_resonance;
          }
      }
    else
       Hermes_ALLM_CLAS(sf);
  }

  //_________________________________________________________________________
  double ProtonStructure::R1998_protected_lowQ2(double const& x, double const& q2, double const& w2) const
  {
    double rvalue = 0;
    const double q2_boundary = 0.34;
    if (q2 < q2_boundary)
      {
        const double x_rescaled = x_from_w2q2(w2, q2);
        rvalue = R1998(x_rescaled, q2_boundary);
        const double u = q2 / q2_boundary;
        rvalue *= 0.5*(3*u-pow(u,3));
      }
    else
      rvalue = R1998(x, q2);
    return rvalue;
  }

  //_________________________________________________________________________
  double ProtonStructure::R1998(const double &x, const double &q2) const
  {
    const double Theta = 1.+12.*q2/(1.+q2)*pow(0.125,2)/(pow(0.125,2) + pow(x,2));
    const double lnQ2_004 = log(q2/0.04);

    // table II from hep-ph/9808028v1
    const double a[6] = { 0.0485, 0.5470, 2.0621,-0.3804, 0.5090,-0.0285 };
    const double b[6] = { 0.0481, 0.6114,-0.3509,-0.4611, 0.7172,-0.0317 };
    const double c[6] = { 0.0577, 0.4644, 1.8288,12.3708,-43.1043,41.7415 };

    const double Ra = a[0]/lnQ2_004*Theta + a[1]/pow(pow(q2,4) + pow(a[2],4), 0.25) * (1 + a[3]*x + a[4]*pow(x,2))*pow(x,a[5]);
    const double Rb = b[0]/lnQ2_004*Theta +(b[1]/q2 + b[2]/(pow(q2,2) + pow(0.3,2)))*(1 + b[3]*x + b[4]*pow(x,2))*pow(x,b[5]);
    const double q2thr = c[3]*x + c[4]*pow(x,2) + c[5]*pow(x,3);
    const double Rc = c[0]/lnQ2_004 * Theta + c[1] * pow( pow(q2 - q2thr, 2) + pow(c[2],2), -0.5);
    return (Ra + Rb + Rc)/3.0;
  }

  //_________________________________________________________________________
  double ProtonStructure::allm_sigma(double const& x, double const& q2, double const& w2) const
  {
    const double hbarc2 = 389.379323;
    double sigmatot, F2, dsigmatot, dF2;
    gd_fit_11_(x,q2,w2,"p",sigmatot,F2,dsigmatot,dF2);
    double sigma = sigmatot + _allm_limits*dsigmatot;
    return sigma/hbarc2;
  }

  //_________________________________________________________________________
  double ProtonStructure::F2_from_sigmaTL(double const& x, double const& q2, double const& sigmaTL) const
  {
    return 1.0/(4*pow(M_PI,2)*_alpha_ref)*q2*(1.0-x)/(1.0+4*pow(x,2)*_mproton2/q2)*sigmaTL;
  }

  //_________________________________________________________________________
  double ProtonStructure::FL_from_F2R(double const& x, double const& q2, double const& F2, double const& rvalue) const
  {
    return F2*(1+4*_mproton2*pow(x,2)/q2)*rvalue/(1.0+rvalue);
  }

  //_________________________________________________________________________
  double ProtonStructure::sigmaTL_from_F2(double const& x, double const& q2, double const& F2) const
  {
    double sigmaTL = 0;
    if (F2 != 0)
      sigmaTL = (4*pow(M_PI,2)*_alpha_ref)*(1+4*pow(x,2)*_mproton2/q2)/(q2*(1.0-x))*F2;
    return sigmaTL;
  }

  //_________________________________________________________________________
  double ProtonStructure::R_from_F2FL(double const& x, double const& q2, double const& F2, double const& FL) const
  {
    const double scaledf2 = F2*(1+4*_mproton2*pow(x,2)/q2);
    return FL/(scaledf2-FL);
  }


}
