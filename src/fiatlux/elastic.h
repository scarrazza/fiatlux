//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <array>
#include <string>
#include <vector>
#include <memory>
using std::array;
using std::string;
using std::vector;
using std::unique_ptr;

#include <fiatlux/integrator.h>
#include <fiatlux/proton.h>

namespace fiatlux
{
  /**
   * @brief The specialized class for the construction of
   * the elastic content of the photon.
   *
   * This class provides the elastic piece of the photon
   * structre following the LUXqed description.
   */
  class ElasticPhoton: public Integrator
  {
  public:
    /**
     * @brief The ElasticPhoton class constructor.
     *
     * Associates the integrand function to the integrator
     * object.
     */
    ElasticPhoton(unique_ptr<ProtonStructure> const& proton);

    /**
     * @brief Evaluates the elastic piece of the photon PDF.
     *
     * @param x the momentum fraction.
     * @return the elastic integral for the photon PDF.
     */
    double evaluatephoton(double const&x, const double &mu2) const;

    /**
     * @brief The elastic integrand.
     *
     * @param q2 the Q^2.
     * @return the elastic integrand.
     */
    double integrand(double const& lnq2, vector<double> const& e) const;

    /**
     * @brief Computes the electric and magnetic form factors.
     *
     * For a given Q^2, return the electric and magnetic form factors
     * in the output array. GM includes the mum_proton factor.
     * This implementation follows arXiv:1307.6227.
     *
     * @param q2 the input energy.
     * @return a 2D stack array with GE (0) and GM (1)
     */
    array<double,2> elastic_ge_gm(double const& q2) const;

    /**
     * @brief The standard elastic dipole factor
     * \f[
     *    G(Q^2) = \left( 1 + \frac{Q^2}{0.71 \text{GeV}^2} \right)^2
     * \f]
     *
     * @param q2 the input energy.
     * @return the elastic dipole factor
     */
    double elastic_dipole_factor(double const& q2) const;

  private:
    unique_ptr<ProtonStructure> const& _proton;
    double _elastic_electric_rescale; //!< the ge rescale.
    double _elastic_magnetic_rescale; //!< the gm rescale.
    int    _elastic_param;            //!< the elastic parametrization.

    vector<array<double, 3>> _fit;         //!< the A1 fit container.
    vector<array<double, 3>> _fit_uperr;   //!< the upper bound fit error.
    vector<array<double, 3>> _fit_downerr; //!< the down bound fit error.
  };
}
