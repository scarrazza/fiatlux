//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <fiatlux/integrator.h>

namespace fiatlux
{
  /**
   * @brief The InelasticQ2 class
   */
  class InelasticQ2: Integrator
  {
  public:
    InelasticQ2(): Integrator{} {}
    double integrand(double const& lnq2) const;
  };

  /**
   * @brief The specialized class for the construction of
   * the inelastic content of the photon.
   *
   * This class provides the inelastic piece of the photon
   * structre following the LUXqed description.
   */
  class InelasticPhoton: Integrator
  {
  public:
    /**
     * @brief The InelasticPhoton class constructor.
     *
     * Associates the integrand function to the integrator
     * object.
     */
    InelasticPhoton();

    /**
     * @brief Evaluates the inelastic piece of the photon PDF.
     *
     * @param x the momentum fraction.
     * @param q2 the energy scale.
     * @return the inelastic integral for the photon PDF.
     */
    double evaluatephoton(double const&x, double const& q2);

    /**
     * @brief The inelastic integrand.
     *
     * @param q2 the Q^2.
     * @return the elastic integrand.
     */
    double integrand(double const& ln1oz) const;

  private:
    double _mproton2;   //!< the square of the proton mass.
    double _eps_base;   //!< precision on final integration of double integral.
    double _eps_rel;    //!< extra precision on any single integration.
    double _alpha_ref;  //!< the reference alpha.
    double _x; //!< tmp storage.
    double _q2; //!< tmp storage.
    double _q2min_inel_override;
    double _q2max_inel_override;
    bool _use_mu2_as_upper_limit; //!< upper bound switcher
    InelasticQ2 _inelastic_q2;
  };
}
