//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include <iostream>
#include <vector>
#include <cmath>
#include <APFEL/APFEL.h>
using namespace std;

// hoppet related stuff
extern "C" double __hoppet_MOD_light_initialize(double const&,double const&, double const&);
extern "C" double __hoppet_MOD_alphaqed(double const&);
extern "C" double __hoppet_MOD_f2(double const&, double const&);
extern "C" double __hoppet_MOD_fl(double const&, double const&);
extern "C" double __hoppet_MOD_masses(int const&);

int main()
{
  const double mcharm = 1.275, mbottom = 4.18, mtop = 173.07;

  // APFEL loading
  APFEL::SetPerturbativeOrder(2);
  APFEL::SetTheory("QUniD");
  APFEL::EnableNLOQEDCorrections(true);
  APFEL::SetPoleMasses(mcharm, mbottom, mtop);
  APFEL::SetAlphaQEDRef(1/137.035999074, 0.000510998946);
  APFEL::SetPDFSet("PDF4LHC15_nnlo_100.LHgrid");
  APFEL::SetQLimits(1, 1e6);
  APFEL::SetQGridParameters(50, 3);
  APFEL::InitializeAPFEL();

  // hoppet loading
  __hoppet_MOD_light_initialize(mcharm, mbottom, mtop);

  // printing alpha
  cout << "\n- Testing alpha evolution:" << endl;
  vector<double> Q = {0.000510998946, 1, 1.777, 10, 100, 1000};
  for (const auto &i: Q)
  {
    const double apfel = APFEL::AlphaQED(i);
    const double hoppet = __hoppet_MOD_alphaqed(i);
    cout << scientific;
    cout << "Q=" << i << "\t apfel=" << apfel << "\t hoppet=" << hoppet << "\t ratio=" << apfel/hoppet << endl;
  }

  return 0;
}
