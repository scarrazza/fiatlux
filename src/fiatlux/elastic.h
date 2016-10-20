//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <array>
#include <string>
using std::array;
using std::string;

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
  class ElasticPhoton: Integrator
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
    double evaluatephoton(double const&x, double const& q2) const;

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
    double _mproton2;  //!< the square of the proton mass.
    double _eps_base;  //!< precision on final integration of double integral.
    double _eps_rel;   //!< extra precision on any single integration.
    double _log_q2_max;//!< the maximum allowed Q2
    double _alpha_ref; //!< the reference alpha
    int _elastic_param;//! the elastic parametrization.
  };
}
