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

double APFELF2(double const& x, double const& Q)
{
  return APFEL::StructureFunctionxQ("EM", "F2", "total", x, Q);
}

double APFELFL(double const& x, double const& Q)
{
  return APFEL::StructureFunctionxQ("EM", "FL", "total", x, Q);
}

int main()
{
  FiatLux lux{"examples/runcard.yml"};

  bool apfel = input().get<bool>("apfel");

  const double mcharm = 1.275, mbottom = 4.18, mtop = 173.07;
  if (apfel)
    {
      cout << "Using APFEL" << endl;
      APFEL::SetPerturbativeOrder(2);
      APFEL::SetTheory("QUniD");
      APFEL::EnableNLOQEDCorrections(true);
      APFEL::SetPoleMasses(mcharm, mbottom, mtop);
      APFEL::SetAlphaQEDRef(1/137.035999074, 0.000510998946);
      APFEL::SetPDFSet("PDF4LHC15_nnlo_100.LHgrid");
      APFEL::SetQLimits(1, 1e6);
      APFEL::SetQGridParameters(50, 3);
      APFEL::InitializeAPFEL_DIS();
      APFEL::CacheStructureFunctionsAPFEL(1);
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
