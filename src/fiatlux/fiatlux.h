//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <string>
#include <memory>
using std::string;
using std::unique_ptr;

#include <fiatlux/elastic.h>
#include <fiatlux/inelastic.h>
#include <fiatlux/msbar.h>
#include <fiatlux/proton.h>

namespace fiatlux
{  
  /**
   * @brief The luxqed struct which holds the 3 components.
   * elastic, inelastic, msbar-pf and the total.
   */
  struct luxqed
  {
    double elastic;
    double inelastic_pf;
    double msbar_pf;
    double total;
  };

  /**
   * @brief The main class of this library.
   *
   * Prepares and controls the procedure to generate the
   * photon PDF using the LUXqed methodology (arXiv:1607.04266).
   */
  class FiatLux
  {
  public:
    /**
     * @brief The fiatlux constructor
     *
     * Takes a PDF provider which returns x*f(fl,x,Q).
     * @param f a std::function object following the \c xfxq definition.
     */
    FiatLux(string const& filename);

    /**
     * @brief Set the external function which returns alphaQED running.
     * @param a the function for alpha running.
     */
    void plug_alphaqed(alpha_running const& a) const;

    /**
     * @brief Evaluates the photon PDF for a given x and Q2.
     *
     * This method computes the 3 components, elastic, inelastic-pf and msbar-pf
     * by performing integral and returns the overall sum in a \c luxqed structure.
     *
     * @param x the momentum fraction.
     * @param q2 the energy scale.
     * @return a luxqed structure with all integral pieces and the total sum.
     */
    luxqed evaluatephoton(double const&x, double const& mu2) const;

  protected:
    void load_settings(string const& filename) const;

  private:
    unique_ptr<ElasticPhoton> _elastic; //!< the integrator for the elastic component.
    unique_ptr<InelasticPhoton> _inelastic; //!< the integrator for the inelastic component.
    unique_ptr<MSbarPhoton> _msbar; //!< the integrator for the msbar component.
  };
}
