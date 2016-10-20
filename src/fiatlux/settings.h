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
  enum {elastic_dipole, elastic_A1_world_spline};

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
     *
     */
    template<class T> T get(string const& key)
    {
      T result = _config[key].as<T>();
      return result;
    }

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
