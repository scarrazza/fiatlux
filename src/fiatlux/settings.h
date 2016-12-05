//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <string>
using std::string;

#include <yaml-cpp/yaml.h>

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
      T result = _config[key].as<T>();
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
    YAML::Node _config;
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
