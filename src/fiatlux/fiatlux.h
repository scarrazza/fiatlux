//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <functional>
#include <string>
#include <memory>
using std::function;
using std::string;
using std::unique_ptr;

#include <fiatlux/elastic.h>
#include <fiatlux/inelastic.h>

namespace fiatlux
{
  /**
   * @brief Typename for parton input.
   *
   * The prototype function should return x*f(fl,x,Q)
   * where fl is PDG id code (integer), x and Q double precision.
   */
  using xfxq = function<double(int,double,double)>;

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
    FiatLux(string const& filename, xfxq const&f);

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
    luxqed evaluatephoton(double const&x, double const& q2) const;

  protected:
    void load_settings(string const& filename) const;

  private:
    xfxq _xfxq; //!< the function which holds the QCD parton information.
    unique_ptr<ElasticPhoton> _elastic; //!< the integrator for the elastic component.
    unique_ptr<InelasticPhoton> _inelastic; //!< the integrator for the inelastic component.
  };
}
