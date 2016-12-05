//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <vector>
using std::vector;

#include <fiatlux/integrator.h>
#include <fiatlux/proton.h>

namespace fiatlux
{
  /**
   * @brief The InelasticQ2 class
   */
  class InelasticQ2: public Integrator, public ProtonStructure
  {
  public:
    InelasticQ2(): Integrator{}, ProtonStructure{} {}
    double integrand(double const& lnq2, const vector<double> &e) const;
  };

  /**
   * @brief The specialized class for the construction of
   * the inelastic content of the photon.
   *
   * This class provides the inelastic piece of the photon
   * structre following the LUXqed description.
   */
  class InelasticPhoton: public Integrator, public ProtonStructure
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

    /**
     * @brief set_alpha_running
     * @param a
     */
    void set_alpha_running(const alpha_running &a) { _inelastic_q2.set_alpha_running(a); ProtonStructure::set_alpha_running(a); }

    /**
     * @brief set_sf
     * @param f2
     * @param fl
     */
    void set_sf(ext_sf const& f2, ext_sf const& fl) { _inelastic_q2.set_sf(f2,fl); ProtonStructure::set_sf(f2,fl); }

  private:
    double _q2min_inel_override;  //!< override value for the q2min
    double _q2max_inel_override;  //!< override value for the q2max
    bool _use_mu2_as_upper_limit; //!< upper bound switcher

    vector<double> _q2_inel_split;
    InelasticQ2 _inelastic_q2;
  };
}
