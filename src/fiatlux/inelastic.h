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
   * @brief The InelasticQ2 class
   */
  class InelasticQ2: public Integrator
  {
  public:
    InelasticQ2(unique_ptr<ProtonStructure> const& p): Integrator{}, _proton(p) {}
    double integrand(double const& lnq2, const vector<double> &e) const;
  private:
    unique_ptr<ProtonStructure> const& _proton;
  };

  /**
   * @brief The specialized class for the construction of
   * the inelastic content of the photon.
   *
   * This class provides the inelastic piece of the photon
   * structre following the LUXqed description.
   */
  class InelasticPhoton: public Integrator
  {
  public:
    /**
     * @brief The InelasticPhoton class constructor.
     *
     * Associates the integrand function to the integrator
     * object.
     */
    InelasticPhoton(unique_ptr<ProtonStructure> const& proton);

    /**
     * @brief Evaluates the inelastic piece of the photon PDF.
     *
     * @param x the momentum fraction.
     * @param q2 the energy scale.
     * @return the inelastic integral for the photon PDF.
     */
    double evaluatephoton(double const&x, double const& mu2) const;

    /**
     * @brief The inelastic integrand.
     *
     * @param q2 the Q^2.
     * @return the elastic integrand.
     */
    double integrand(double const& ln1oz, vector<double> const& e) const;

    /**
     * @brief Inserts a inelastic split.
     *
     * This function splits the final integral in
     * subintegrals which are summed in the integrand method.
     * @param q2 the energy scale.
     */
    void insert_inel_split(double const& q2);

    /**
     * @brief Returns the array containing the inelastic splits.
     *
     * @return the constant referece to the vector containing the splits.
     */
    const vector<double>& get_q2_inel_split() const { return _q2_inel_split; }   

  private:
    unique_ptr<ProtonStructure> const& _proton; //!< the proton container
    double _q2min_inel_override;  //!< override value for the q2min
    double _q2max_inel_override;  //!< override value for the q2max
    bool _use_mu2_as_upper_limit; //!< upper bound switcher

    vector<double> _q2_inel_split;
    InelasticQ2 _inelastic_q2;
  };
}
