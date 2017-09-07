//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <iostream>
#include <functional>
using std::function;

namespace fiatlux
{
  /**
   * @brief Typename for alpha running alpha(Q)
   */
  using alpha_running = function<double(double const&)>;

  /**
   * @brief Structure function F(x,Q)
   */
  using ext_sf = function<double(double const&, double const&)>;

  /**
   * @brief A simple container for the output results
   */
  struct StrucFunc
  {
    double x, q2, w2;
    double FL, F2, F2LO;
  };

  /**
   * @brief The ProtonStructure class.
   *
   * This class provides two functionalities to the library:
   * - a set of common attributes from the configuration file related
   * to the proton structure and general setup of the code.
   * - the methods to compute the structure function based on the inelastic
   * parametrization selected by the user.
   */
  class ProtonStructure
  {
  public:
    ProtonStructure();

    /**
     * @brief Provides access to the proton structure functions.
     *
     * @param x the momentum fraction.
     * @param q2 the energy scale
     * @return a StrucFunc which contains the F2 and FL structure functions.S
     */
    StrucFunc compute_proton_structure(double const& x, double const& q2, bool f2lo_only = false) const;

    /**
     * @brief computes W^2 given x and Q^2
     * @param x the momentum fraction
     * @param q2 the energy scale
     * @return the W^2
     */
    double w2_from_xq2(double const& x, double const& q2) const;

    /**
     * @brief computes x given W^2 and Q^2
     */
    double x_from_w2q2(double const& w2, double const& q2) const;

    /**
     * @brief Takes Rvalue as input and modifies it, scaling it by a factor
     * involving two types of rescaling -- a constant (below
     * LHAPDF_transition_Q2) and a higher twist (above).
     * R_rescale_twist4 corresponds to the parameter A of 1605.08577
     */
    void rescaler(double &rvalue, double const& q2) const;

    /**
     * @brief Set the function for the alpha running
     * @param a the alpha function
     */
    void set_alpha_running(alpha_running const& a) { _alpha_running = a; }

    /**
     * @brief Set the alpha value at reference scale
     * @param aref the alpha value at reference scale
     */
    void set_alpha_ref(double aref) { _alpha_ref = aref; }

    /**
     * @brief Set the structure functions for the F2(x,Q) and FL(x,Q)
     * @param f2 the F2(x,Q) structure function
     * @param fl the FL(x,Q) structure function
     */
    void set_sf(ext_sf const& f2, ext_sf const& fl, ext_sf const& f2lo) { _f2 = f2; _fl = fl; _f2lo = f2lo; }

  private:
    /**
     * @brief Computes the structure functions using the Hermes_ALLM_CLAS model at low energy.
     * @param sf the structure function struct which will be filled by the code.
     */
    void Hermes_ALLM_CLAS(StrucFunc & sf) const;

    /**
     * @brief Computes the structure functions using the Hermes_ALLM_CLAS model + PDFs from LHAPDF for high scales.
     * @param sf the structure function struct which will be filled by the code.
     */
    void LHAPDF_Hermes_ALLM_CLAS(StrucFunc & sf, bool f2lo_only = false) const;

    /**
     * @brief Protected version of R1998 where Q2 is checked before computing.
     * If its value is smaller than a threshold interpolation is used.
     */
    double R1998_protected_lowQ2(double const& x, double const& q2, const double &w2) const;

    /**
     * @brief R1998 parametrisation of R = sigmaL/sigmaT as given
     * in hep-ph/9808028 from a fit to world data, their Eqs. 2 and 3.
     *
     * Q2 is assumed to be in GeV^2 and this appears not to be sensible
     * (crazy values) for small Q2, ie roughly below 0.3GeV^2 (where it is not
     * intended to be valid in the R1990 paper, Phys.Lett. B250 (1990) 193-198
     *
     * JLAB has a bunch of publications related to L/T
     * https://hallcweb.jlab.org/resdata/publications/
     * but the critical one, a Ph.D. thesis by Tvaskis,
     * Longitudinal-Transverse Separation of Deep-Inelastic Scattering at low Q2
     * on Nucleons and Nuclei seems not to be available...
     */
    double R1998(double const& x, double const& q2) const;

    /**
     * @brief Computes sigmaTL in GeV^-2
     */
    double allm_sigma(double const& x, double const& q2, const double &w2) const;

    /**
     * @brief Computes F2 from sigmaTL
     */
    double F2_from_sigmaTL(double const& x, double const& q2, double const& sigmaTL) const;

    /**
     * @brief Computes FL from F2 and R
     */
    double FL_from_F2R(double const& x, double const& q2, double const& F2, double const& rvalue) const;

    /**
     * @brief Computes sigmaTL from F2
     */
    double sigmaTL_from_F2(double const& x, double const& q2, double const& F2) const;

    /**
     * @brief Computes R for a given F2 and FL.
     */
    double R_from_F2FL(double const& x, double const& q2, double const& F2, double const& FL) const;

  public:
    alpha_running _alpha_running; //!< the alpha running function
    ext_sf _f2; //!< the external F2 function for Q2 > transition
    ext_sf _fl; //!< the external FL function for Q2 > transition
    ext_sf _f2lo; //!< the external F2 at LO function for the MSbar piece

    const bool _qed_running;  //!< switch qed alpha running
    double _alpha_ref;        //!< the reference alpha.
    const double _mproton2;   //!< the proton mass^2
    const double _mum_proton; //!< proton magnetic momentum.
    const double _eps_base;   //!< precision on final integration of double integral.
    const double _eps_rel;    //!< extra precision on any single integration.
    const double _log_q2_max; //!< the maximum allowed Q2.
    const double _rescale_r;  //!< the rescale coeff.
    const double _rescale_r_twist4; //!< the twist rescale
    const double _lhapdf_transition_q2; //!<    
    const double _allm_limits;
    const double _rescale_non_resonance;
    const double _rescale_resonance;
    const double _HAC_loW2; //!< low W2 CLAS F2
    const double _HAC_hiW2; //!< high W2 Hermes ALLM sigmaTL
    const int _inelastic_param;   //!< the inelastic parametrization
  };
}
