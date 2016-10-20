//
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#include "fiatlux/integrator.h"
#include "fiatlux/tools.h"

#include <gsl/gsl_errno.h>

namespace fiatlux
{
  struct int_param { const Integrator *it; double extra; };
  double int_gsl(double x, void *p)
  {
    struct int_param * params = (struct int_param *) p;
    return params->it->integrand(x, params->extra);
  }

  //_________________________________________________________________________
  Integrator::Integrator()
  {
    // allocate gsl integration workspace
    _gslwork = gsl_integration_workspace_alloc(10000);
  }

  //_________________________________________________________________________
  Integrator::~Integrator()
  {
    // free gsl workspace
    gsl_integration_workspace_free(_gslwork);
  }

  //_________________________________________________________________________
  double Integrator::integrate(const double &xmin, const double &xmax, const double &eps, const double &extra) const
  {
    double int_res = 0;
    double int_err = 0;

    // gsl parameters
    int_param gslparam = {this, extra};

    // gsl function
    gsl_function F;
    F.function = &int_gsl;
    F.params = &gslparam;

    // integration
    int status = gsl_integration_qags(&F, xmin, xmax, 0, eps, _gslwork->limit, _gslwork, &int_res, &int_err);

    if (status == GSL_EDIVERGE || status == GSL_ESING || status == GSL_EROUND)
      info("Integrator::integrate", "GSL error " + std::to_string(status));

    return int_res;
  }
}
