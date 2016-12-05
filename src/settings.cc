//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/settings.h"
#include "fiatlux/tools.h"

namespace fiatlux
{
  //_________________________________________________________________________
  Settings& input()
  {
    static Settings settings{};
    return settings;
  }

  //_________________________________________________________________________
  void Settings::load(string const& filename)
  {
    try
    {
      _config = YAML::LoadFile(filename);
    }
    catch (YAML::BadFile &)
    {
      throw runtime_exception("Settings::load", "cannot find file " + filename);
    }
  }

  //_________________________________________________________________________
  int Settings::get_elastic_param()
  {
    const string eopt = get<string>("elastic_param");
    if ( eopt.compare("dipole") == 0)
      return elastic_dipole;
    else if (eopt.compare("A1_world_spline") == 0)
      return elastic_A1_world_spline;
    else if (eopt.compare("A1_world_pol_spline") == 0)
      return elastic_A1_world_pol_spline;
    else
      throw runtime_exception("Settings::get_elastic_param", "option not recognised");
  }

  //_________________________________________________________________________
  int Settings::get_inelastic_param()
  {
    const string eopt = get<string>("inelastic_param");
    if ( eopt.compare("Hermes_ALLM_CLAS") == 0)
      return inelastic_Hermes_ALLM_CLAS;
    else if (eopt.compare("LHAPDF_Hermes_ALLM_CLAS") == 0)
      return inelastic_LHAPDF_Hermes_ALLM_CLAS;
    else
      throw runtime_exception("Settings::get_inelastic_param", "option not recognised");
  }
}
