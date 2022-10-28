//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <string>
#include <sstream>
#include <map>
using std::string;
using std::stringstream;
using std::map;

namespace fiatlux
{
  /**
   * @brief The eparam enum
   */
  enum {elastic_dipole, elastic_A1_world_spline, elastic_A1_world_pol_spline};
  enum {inelastic_Hermes_ALLM_CLAS, inelastic_LHAPDF_Hermes_ALLM_CLAS};

  /**
   * @brief The Settings class.
   *
   * This class stores and manipulates the configuration
   * options of the library.
   */
  class Settings
  {
  public:
    Settings() {}

    /**
     * @brief load
     * @param filename
     */
    void load(string const& filename);

    /**
     * @brief get
     * @param key
     * @return
     */
    template<class T>
    T get(string const& key)
    {
      stringstream ss;
      T result;
      string val = _config.at(key);
      if (val == "true" || val == "on" || val == "True") val = "1";
      if (val == "false" || val == "off" || val == "False") val = "0";
      ss << val;
      ss >> result;
      return result;
    }

    /**
     * @brief get_elastic_param
     * @return
     */
    int get_elastic_param();

    /**
     * @brief get_inelastic_param
     * @return
     */
    int get_inelastic_param();

    /**
     * @brief print
     */
    void print() const;

  private:
    map<string, string> _config;
  };

  /**
   * @brief A singleton alternative to access Settings.
   *
   * By calling s(), this function returns an unique reference
   * to the Settings class and its attributes.
   *
   * @return a single instance of the Settings class, available from
   * all parts of the code.
   */
  Settings& input();
}
