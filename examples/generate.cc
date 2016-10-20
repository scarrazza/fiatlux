//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include <iostream>
#include <vector>
using namespace std;

#include <fiatlux/fiatlux.h>
using namespace fiatlux;
using namespace std::placeholders;

#include <LHAPDF/LHAPDF.h>
using namespace LHAPDF;
double (PDF::*fn)(int, double, double) const = &PDF::xfxQ; // extracting LHAPDF methods

int main()
{
  PDF *pdf = mkPDF("PDF4LHC15_nnlo_100", 0);
  FiatLux lux{"examples/runcard.yml", std::bind(fn, pdf, _1, _2, _3)};

  vector<double> x = {0.90483741668764672};
  vector<double> q2= {2.0};

  cout << "Columns: x Q2 elastic" << endl;
  for (auto const& iq2: q2)
    for (auto const& ix: x)
      {
        const auto pht = lux.evaluatephoton(ix,iq2);
        cout << scientific << ix << "\t" << iq2 << "\t" << pht.elastic << endl;
      }

  return 0;
}
