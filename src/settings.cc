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
}
