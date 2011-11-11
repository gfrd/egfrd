// Function findRoot iterates the GSL root finder until a root has been found
// with the requested precision.
//
// Author, amongst others: Laurens Bossen.
// FOM Insitute AMOLF


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdexcept>
#include <gsl/gsl_errno.h>

#include "Logger.hpp"
#include "findRoot.hpp"

// Iterates the solver until desired precision has been reached or a maximum
// number of iterations have been performed.
Real findRoot(gsl_function const& F, gsl_root_fsolver* solver, Real low,
              Real high, Real tol_abs, Real tol_rel, char const* funcName)
{
    // low and high should constitute an interval that straddles a root.
    Real l(low);
    Real h(high);

    gsl_root_fsolver_set(solver, const_cast<gsl_function*>(&F), l, h);

    const unsigned int maxIter(100);

    unsigned int i(0);
    for (;;)
    {
    
        // iterate        
        gsl_root_fsolver_iterate(solver);
               
        // get found bracketing interval
        l = gsl_root_fsolver_x_lower(solver);
        h = gsl_root_fsolver_x_upper(solver);     

        // see if this is acceptable
        const int status(gsl_root_test_interval(l, h, tol_abs,
                                                  tol_rel));

        // stop finder if convergence or too much iterations
        if (status == GSL_CONTINUE)
        {
            if (i >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                throw std::runtime_error(std::string(funcName) + ": failed to converge");
            }
        }
        else
        {
            break;
        }

        ++i;
    }
  

    const Real root(gsl_root_fsolver_root(solver));
        
    return root;    
}





