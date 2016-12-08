//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fiatlux/fiatlux.h>
#include <fiatlux/settings.h>
#include <APFEL/APFEL.h>
#include <LHAPDF/LHAPDF.h>
using namespace std;
using namespace fiatlux;

double y_of_zeta(double const& y, double const& a)
{
  return y + a*(1.0-exp(-y));
}

// hoppet related stuff
extern "C" double __hoppet_MOD_initialize(double const&,double const&, double const&);
extern "C" double __hoppet_MOD_alphaqed(double const&);
extern "C" double __hoppet_MOD_f2(double const&, double const&);
extern "C" double __hoppet_MOD_fl(double const&, double const&);
extern "C" double __hoppet_MOD_masses(int const&);

LHAPDF::PDF* _lhapdf;
extern "C" void externalsetapfel_(double const& x, double const& Q, double *xf)
{
  for (int i = 0; i < 13; i++)
    xf[i] = _lhapdf->xfxQ(i-6, x, Q);
  xf[13] = 0.0;
}

double APFELF2(double const& x, double const& Q)
{
  APFEL::ComputeStructureFunctionsAPFEL(Q,Q);
  return APFEL::F2total(x);
}

double APFELFL(double const& x, double const& Q)
{
  APFEL::ComputeStructureFunctionsAPFEL(Q,Q);
  return APFEL::FLtotal(x);
}

int main()
{
  FiatLux lux{"examples/runcard.yml"};

  bool apfel = input().get<bool>("apfel");

  const double mcharm = 1.275, mbottom = 4.18, mtop = 173.07;
  if (apfel)
    {
      cout << "Using APFEL" << endl;
      APFEL::SetPerturbativeOrder(0);
      APFEL::SetTheory("QUniD");
      APFEL::EnableNLOQEDCorrections(true);
      APFEL::SetMSbarMasses(mcharm, mbottom, mtop);
      APFEL::SetAlphaQEDRef(1/137.035999074, 0.000510998946);
      _lhapdf = LHAPDF::mkPDF("PDF4LHC15_nnlo_100");
      APFEL::SetPDFSet("external");
      //APFEL::SetMassScheme("ZMFV");
      //APFEL::SetProcessDIS("NC");
      APFEL::SetQLimits(1, 1e100);
      APFEL::InitializeAPFEL_DIS();
      lux.PlugAlphaQED(APFEL::AlphaQED);
      lux.PlugStructureFunctions(APFELF2, APFELFL);
      lux.InsertInelasticSplitQ({mbottom, 1e100});
    }
  else
    {
      cout << "Using HOPPET" << endl;
      __hoppet_MOD_initialize(mcharm, mbottom, mtop);
      lux.PlugAlphaQED(__hoppet_MOD_alphaqed);
      lux.PlugStructureFunctions(__hoppet_MOD_f2, __hoppet_MOD_fl);
      lux.InsertInelasticSplitQ({__hoppet_MOD_masses(5),__hoppet_MOD_masses(6)});
    }

  /*
  __hoppet_MOD_initialize(mcharm, mbottom, mtop);
  cout << APFEL::AlphaQED(0.000510998946) << " " << __hoppet_MOD_alphaqed(0.000510998946) << endl;
  cout << APFEL::AlphaQED(1.777) << " " << __hoppet_MOD_alphaqed(1.777) << endl;
  cout << APFEL::AlphaQED(81.9) << " " << __hoppet_MOD_alphaqed(81.9) << endl;  
  */

  // print results
  double q2 = 10000;
  cout << setprecision(15) << scientific;
  for (auto i = 1; i <= 100; i++)
    {
      const auto y = y_of_zeta(i*0.1, 0);
      const auto x = exp(-y);

      const auto pht = lux.EvaluatePhoton(x,q2);
      cout << x << "\t"
           << q2 << "\t"
           << pht.elastic << "\t"
           << pht.inelastic_pf << "\t"
           << pht.msbar_pf << "\t"
           << pht.total << "\t"
           << endl;
    }

  return 0;
}
