//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

namespace fiatlux
{
  /**
   * @brief The Settings class.
   *
   * This class stores and manipulates the configuration
   * options of the library.
   */
  class Settings
  {
  public:
    Settings();
    bool qed_running; //!< determines the running of alpha.
    double q2_max;   //!< the maximum allowed Q2.
    double eps_base; //!< precision on final integration of double integral.
    double eps_rel;  //!< extra precision on any single integration.
    double mproton;  //!< the proton mass.
    double alpha_ref;//!< the reference alpha constant.
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
  Settings& s();
}
