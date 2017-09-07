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
     * @param qref the reference energy scale (relevant when qed_running is false)
     */
    void PlugAlphaQED(alpha_running const& a, double qref = 0.000510998946) const;

    /**
     * @brief Set the external function which returns F2(x,Q) and FL(x,Q).
     * @param f2 the F2(x,Q) structure function.
     * @param fl the FL(x,Q) structure function.
     * @param f2lo the F2(x,Q) at leading order.
     */
    void PlugStructureFunctions(ext_sf const& f2, ext_sf const& fl, ext_sf const& f2lo) const;

    /**
     * @brief Splits inelastic integral into small pieces in order to avoid
     * discontinuity regions. Optional, useful when dealing with complicated sf.
     * @param thresholds a vector with the Q knots for the split.
     */
    void InsertInelasticSplitQ(vector<double> const& qthresholds) const;

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
    luxqed EvaluatePhoton(double const&x, double const& mu2) const;

  protected:
    void load_settings(string const& filename) const;

  private:
    unique_ptr<ProtonStructure> _proton; //!< the proton structure used by all components
    unique_ptr<ElasticPhoton> _elastic; //!< the integrator for the elastic component.
    unique_ptr<InelasticPhoton> _inelastic; //!< the integrator for the inelastic component.
    unique_ptr<MSbarPhoton> _msbar; //!< the integrator for the msbar component.
  };
}
