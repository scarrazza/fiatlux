//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include <iostream>
using namespace std;

#include <fiatlux/fiatlux.h>
using namespace fiatlux;
using namespace std::placeholders;

#include <LHAPDF/LHAPDF.h>
using namespace LHAPDF;
double (PDF::*xfxq)(int, double, double) const = &PDF::xfxQ; // extracting LHAPDF methods

int main()
{
  PDF *pdf = mkPDF("PDF4LHC15_nnlo_100", 0);
  FiatLux lux(std::bind(xfxq, pdf, _1, _2, _3));

  return 0;
}
