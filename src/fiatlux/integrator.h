//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <functional>
#include <array>
using std::function;
using std::array;

namespace fiatlux
{
  /**
   * @brief The Integrator class which uses Adaptative Gaussian Quadrature.
   *
   * This class takes as input the integrand function and provides
   * the integrate method which performs the integration.
   */
  class Integrator
  {
  public:
    Integrator();  //!< The class constructor

    /**
     * @brief Integrates the integrand passed during initialization
     * between xmin and xmax with tolerance eps.
     *
     * @param xmin the lower bound integration value.
     * @param xmax the upper bound integration value.
     * @param eps the required relative error.
     * @return the integral value.
     */
    double integrate(double const& xmin, double const& xmax, double const& eps = 1E-4, double const& extra = 0) const;

  protected:
    /**
     * @brief Protected virtual integrand function.
     *
     * @param x the integration variable.
     * @param extra an optional extra double
     * @return the integrand evaluated at x.
     */
    virtual double integrand(double const& x, double const& extra) const = 0;

    /**
     * @brief The dgauss integrator from cernlib.
     *
     * @param a the lower integration bound.
     * @param b the upper integratio bound.
     * @param eps the required accuracy.
     * @param extra optional parameter for 1D integrand which require extra information.
     * @return the integral from a to b.
     */
    double dgauss(double const& a, double const& b, double const& eps, double const& extra) const;

  private:
    double _cst;         //!< param for dgauss
    array<double,12> _w; //!< param for dgauss
    array<double,12> _x; //!< param for dgauss
  };
}
