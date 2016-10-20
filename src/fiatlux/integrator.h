//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <functional>
using std::function;

#include <gsl/gsl_integration.h>

namespace fiatlux
{
  /**
   * @brief The Integrator class which uses GSL (QAG).
   *
   * This class takes as input the integrand function and provides
   * the integrate method which performs the integration.
   */
  class Integrator
  {
  public:
    Integrator();  //!< The class constructor
    ~Integrator(); //!< The class destructor

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
     * @brief Protected auxiliary function which provides
     * access to the GSL integration routines in a C++ style.
     *
     * @param x the integration variable.
     * @param extra an optional extra double
     * @return the integrand evaluated at x.
     */
    virtual double integrand(double const& x, double const& extra) const = 0;

  private:
    gsl_integration_workspace * _gslwork; //!< the gsl workpace for integration.
    friend double int_gsl(double x, void *p); //!< the auxiliary function for gsl integration.
  };
}
