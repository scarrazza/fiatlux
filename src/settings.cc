//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include "fiatlux/settings.h"
#include "fiatlux/tools.h"
using std::fstream;
using std::ios;
using std::vector;
using std::istream_iterator;

namespace fiatlux
{
  //_________________________________________________________________________
  vector<string> split(string const &input)
  {
        stringstream strstr(input);
        istream_iterator<string> it(strstr);
        istream_iterator<string> end;
        vector<string> results(it, end);
        return results;
  }

  //_________________________________________________________________________
  Settings& input()
  {
    static Settings settings{};
    return settings;
  }

  //_________________________________________________________________________
  void Settings::load(string const& filename)
  {
    fstream f;
    f.open(filename.c_str(), ios::in);
    if (!f.good())
      throw runtime_exception("Settings::load", "cannot find file " + filename);

    string line;
    while (f.good())
    {
      getline(f, line);
      vector<string> v = split(line);
      vector<string> entry;
      for (const auto& i: v)
      {
        if (i == "#") break;
        entry.push_back(i);
      }
      if (entry.size() != 2 && entry.size() != 0)
        throw runtime_exception("Settings::load", "syntax error for entry " + entry[0]);
      if (entry.size() == 2)
        _config[entry[0].substr(0, entry[0].size() - 1)] = entry[1];
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
