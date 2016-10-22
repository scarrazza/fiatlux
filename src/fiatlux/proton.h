//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

namespace fiatlux
{
  /**
   * @brief A simple container for the output results
   */
  struct StrucFunc
  {
    double x, q2, y, s, nu_esw, nu_exp, kinw2, ein, eout, theta;
    double FL, F2;
  };

  /**
   * @brief The ProtonStructure class.
   *
   * This class provides two functionalities to the library:
   * - a set of common attributes from the configuration file related
   * to the proton structure and general setup of the code.
   * - the methods to compute the structure function based on the inelastic
   * parametrization selected by the user.
   */
  class ProtonStructure
  {
  public:
    ProtonStructure();

    /**
     * @brief Provides access to the proton structure functions.
     *
     * @param x the momentum fraction.
     * @param q2 the energy scale
     * @return a StrucFunc which contains the F2 and FL structure functions.S
     */
    StrucFunc compute_proton_structure(double const& x, double const& q2) const;

  protected:
    const double _mproton2;   //!< the proton mass^2
    const double _mum_proton; //!< proton magnetic momentum.
    const double _eps_base;   //!< precision on final integration of double integral.
    const double _eps_rel;    //!< extra precision on any single integration.
    const double _log_q2_max; //!< the maximum allowed Q2.
    const double _alpha_ref;  //!< the reference alpha.
  };
}
