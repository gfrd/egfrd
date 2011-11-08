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

#include <iostream> // ONLY FOR DEBUG PURPOSES TODO

#include "Logger.hpp"
#include "findRoot.hpp"

Real findRoot(gsl_function const& F, gsl_root_fsolver* solver, Real low,
              Real high, Real tol_abs, Real tol_rel, char const* funcName)
{

    // TODO DEBUG REMOVE!!!    
    // std::cout << "Starting root finding. Initiated by: " << funcName << ";\n";
    // END TODO DEBUG REMOVE

    Real l(low);
    Real h(high);

    gsl_root_fsolver_set(solver, const_cast<gsl_function*>(&F), l, h);

    const unsigned int maxIter(100);

    unsigned int i(0);
    for (;;)
    {
    
        
        gsl_root_fsolver_iterate(solver);
               
        l = gsl_root_fsolver_x_lower(solver);
        h = gsl_root_fsolver_x_upper(solver);

        // TODO DEBUG REMOVE!!!
        // std::cout << "Current search: [" << l << ", " << h << "];\n";
        // END TODO DEBUG REMOVE
        

        const int status(gsl_root_test_interval(l, h, tol_abs,
                                                  tol_rel));

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

    // TODO DEBUG REMOVE!!!    
    // std::cout << "Final output: [" << root << "];\n";
    // END TODO DEBUG REMOVE
        
    return root;    
}





