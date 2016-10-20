//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include <iostream>
#include <vector>
using namespace std;

#include <fiatlux/fiatlux.h>
#include <fiatlux/settings.h>
using namespace fiatlux;
using namespace std::placeholders;

#include <LHAPDF/LHAPDF.h>
using namespace LHAPDF;
double (PDF::*fn)(int, double, double) const = &PDF::xfxQ; // extracting LHAPDF methods

double y_of_zeta(double const& y, double const& a)
{
  return y + a*(1.0-exp(-y));
}

int main()
{
  //PDF *pdf = mkPDF("PDF4LHC15_nnlo_100", 0);
  FiatLux lux{"examples/runcard.yml",NULL};

  vector<double> q2= {input().get<double>("q2")};

  cout << setprecision(15) << scientific;

  //cout << "Columns: x Q2 elastic inelastic_PF MSbar-PF total" << endl;
  for (auto const& iq2: q2)
    for (auto i = 1; i <= 100; i++)
      {
        const auto y = y_of_zeta(i*0.1, 0);
        const auto x = exp(-y);

        const auto pht = lux.evaluatephoton(x,iq2);
        cout << x << "\t"
             << iq2 << "\t"
             << pht.elastic << "\t"
             << pht.inelastic_pf << "\t"
             << pht.msbar_pf << "\t"
             << pht.total << "\t"
             << endl;
      }

  return 0;
}
