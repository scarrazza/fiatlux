//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <vector>
#include <memory>
using std::vector;
using std::unique_ptr;

#include <fiatlux/integrator.h>
#include <fiatlux/proton.h>

namespace fiatlux
{
  /**
   * @brief The specialized class for the construction of
   * the msbar correction to the photon.
   *
   * This class provides the msbar piece of the photon
   * structre following the LUXqed description.
   */
  class MSbarPhoton: public Integrator
  {
  public:
    MSbarPhoton(unique_ptr<ProtonStructure> const& proton);

    /**
     * @brief Evaluates the msbar piece of the photon PDF.
     *
     * @param x the momentum fraction.
     * @return the elastic integral for the photon PDF.
     */
    double evaluatephoton(double const&x, double const& mu2) const;

    /**
     * @brief The msbar integrand.
     *
     * @param q2 the Q^2.
     * @return the elastic integrand.
     */
    double integrand(double const& ln1oz, vector<double> const& e) const;

  private:
    unique_ptr<ProtonStructure> const& _proton;
    bool _use_mu2_as_upper_limit; //!< upper bound switcher
  };
}
