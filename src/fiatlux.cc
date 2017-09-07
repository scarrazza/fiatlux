//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/fiatlux.h"
#include "fiatlux/tools.h"
#include "fiatlux/settings.h"
#include "fiatlux/config.h"
#include <exception>
#include <cmath>

namespace fiatlux
{
  //_________________________________________________________________________
  FiatLux::FiatLux(const string &filename)
  {

    cout << "   _ _ _     ______ _       _   _                  " << endl;
    cout << "  | (_) |   |  ____(_)     | | | |                " << endl;
    cout << "  | |_| |__ | |__   _  __ _| |_| |    _   ___  __ " << endl;
    cout << "  | | | '_ \\|  __| | |/ _` | __| |   | | | \\ \\/ / " << endl;
    cout << "  | | | |_) | |    | | (_| | |_| |___| |_| |>  <  " << endl;
    cout << "  |_|_|_.__/|_|    |_|\\__,_|\\__|______\\__,_/_/\\_\\ " << endl;
    cout << "        __version__ " << VERSION << " : S. Carrazza" << endl;
    cout << endl;

    input().load(filename);
    _proton= unique_ptr<ProtonStructure>(new ProtonStructure{});
    _elastic = unique_ptr<ElasticPhoton>(new ElasticPhoton{_proton});
    _inelastic = unique_ptr<InelasticPhoton>(new InelasticPhoton{_proton});
    _msbar = unique_ptr<MSbarPhoton>(new MSbarPhoton{_proton});

    _inelastic->insert_inel_split(-2);
    const double lhtrq2 = input().get<double>("lhapdf_transition_q2");
    if (lhtrq2 > 0)
      _inelastic->insert_inel_split(lhtrq2);
    _inelastic->insert_inel_split(1E200);  
  }

  //_________________________________________________________________________
  void FiatLux::PlugAlphaQED(alpha_running const& a, double qref) const
  {
    _proton->set_alpha_running(a);
    _proton->set_alpha_ref(a(qref));
  }

  //_________________________________________________________________________
  void FiatLux::PlugStructureFunctions(ext_sf const& f2, ext_sf const& fl, ext_sf const& f2lo) const
  {
    _proton->set_sf(f2,fl,f2lo);
  }

  //_________________________________________________________________________
  void FiatLux::InsertInelasticSplitQ(vector<double> const& qthresholds) const
  {
    for (auto const& v: qthresholds)
      _inelastic->insert_inel_split(pow(v, 2));
  }

  //_________________________________________________________________________
  luxqed FiatLux::EvaluatePhoton(const double &x, const double &mu2) const
  {
    luxqed e;
    e.elastic = _elastic->evaluatephoton(x, mu2);
    e.inelastic_pf = _inelastic->evaluatephoton(x, mu2);
    e.msbar_pf = _msbar->evaluatephoton(x, mu2);
    e.total = e.elastic + e.inelastic_pf + e.msbar_pf;
    return e;
  }
}
