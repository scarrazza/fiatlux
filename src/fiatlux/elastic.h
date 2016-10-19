//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <array>
using std::array;

#include <fiatlux/integrator.h>

namespace fiatlux
{
  /**
   * @brief The specialized class for the construction of
   * the elastic content of the photon.
   *
   * This class provides the elastic piece of the photon
   * structre following the LUXqed description.
   */
  class ElasticPhoton
  {
  public:
    /**
     * @brief The ElasticPhoton class constructor.
     *
     * Associates the integrand function to the integrator
     * object.
     */
    ElasticPhoton();

    /**
     * @brief Evaluates the elastic piece of the photon PDF.
     *
     * @param x the momentum fraction.
     * @param q2 the energy scale.
     * @return the elastic integral for the photon PDF.
     */
    double evaluate(double const&x, double const& q2) const;

    /**
     * @brief The elastic integrand.
     *
     * @param q2 the Q^2.
     * @return the elastic integrand.
     */
    double integrand(double const& lnQ2, double const& x) const;

    /**
     * @brief Computes the electric and magnetic form factors.
     *
     * For a given Q^2, return the electric and magnetic form factors
     * in the output array. GM includes the mum_proton factor.
     * @param q2 the input energy.
     * @return a 2D stack array with GE (0) and GM (1)
     */
    array<double,2> elastic_ge_gm(double const& q2) const;

  private:
    Integrator _integrator; //!< the integrator for the elastic component.
  };
}
