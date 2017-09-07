//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <APFEL/APFEL.h>
using namespace std;

// hoppet related stuff
extern "C" double __hoppet_MOD_initialize(double const&,double const&, double const&);
extern "C" double __hoppet_MOD_alphaqed(double const&);
extern "C" double __hoppet_MOD_f2(double const&, double const&);
extern "C" double __hoppet_MOD_fl(double const&, double const&);
extern "C" double __hoppet_MOD_f2lo(double const&, double const&);
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
  const double mcharm = 1.275, mbottom = 4.18, mtop = 173.07;

  // loading apfel
  APFEL::SetPerturbativeOrder(2);
  APFEL::SetTheory("QUniD");
  APFEL::EnableNLOQEDCorrections(true);
  APFEL::EnableSFNLOQEDCorrections(false);
  APFEL::SetAlphaQCDRef(0.118, 91.2);
  APFEL::SetPoleMasses(mcharm, mbottom, mtop);
  APFEL::SetAlphaQEDRef(1/137.035999074, 0.000510998946);
  APFEL::SetPDFSet("PDF4LHC15_nnlo_100.LHgrid");
  APFEL::SetQLimits(1, 1e6);
  APFEL::SetQGridParameters(50, 3);
  APFEL::SetProcessDIS("NC");
  APFEL::InitializeAPFEL_DIS();
  APFEL::CacheStructureFunctionsAPFEL(-1);
  APFEL::CachePDFsAPFEL(-1);

  // loading hoppet
  cout << "Using HOPPET" << endl;
  __hoppet_MOD_initialize(mcharm, mbottom, mtop);

  // Testing SF
  vector<double> Q = {1, 5, 10, 15, 100};
  vector<double> x = {1e-3, 1e-1, 0.9};

  cout << scientific;
  cout << "\nTesting F2" << endl;
  for (const auto& iq: Q)
    for (const auto& ix: x)
    {
      const double apfel = APFELF2(ix, iq);
      const double hoppet = __hoppet_MOD_f2(ix, iq);
      cout << "Q=" << iq << "\tx=" << ix << "\tF2 apfel=" << apfel << "\thoppet=" << hoppet << "\tratio=" << apfel/hoppet << endl;
    }

  cout << "\nTesting FL" << endl;
  for (const auto& iq: Q)
    for (const auto& ix: x)
    {
      const double apfel = APFELFL(ix, iq);
      const double hoppet = __hoppet_MOD_fl(ix, iq);
      cout << "Q=" << iq << "\tx=" << ix << "\tFL apfel=" << apfel << "\thoppet=" << hoppet << "\tratio=" << apfel/hoppet << endl;
    }

  cout << "\nTesting F2LO" << endl;
  for (const auto& iq: Q)
    for (const auto& ix: x)
    {
      const double apfel = APFEL::F2LO(ix, iq);
      const double hoppet = __hoppet_MOD_f2lo(ix, iq);
      cout << "Q=" << iq << "\tx=" << ix << "\tF2LO apfel=" << apfel << "\thoppet=" << hoppet << "\tratio=" << apfel/hoppet << endl;
    }

  return 0;
}
