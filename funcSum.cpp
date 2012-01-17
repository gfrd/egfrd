#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <cmath>
#include <boost/bind.hpp>
#include <gsl/gsl_sum.h>

#include "Logger.hpp"
#include "funcSum.hpp"

typedef std::vector<Real> RealVector;

static Logger& _log(Logger::get_logger("funcSum"));

Real 
funcSum_all(boost::function<Real(unsigned int i)> f, size_t max_i)
{
    Real sum(0.0);

    const Real p_0(f(0));
    if (p_0 == 0.0)
    {
        return 0.0;
    }

    sum = p_0;

    RealVector::size_type i(1); 
    while(i < max_i)
    {
        const Real p_i(f(i));
        sum += p_i;

        ++i;
    }

    return sum;
}


Real 
funcSum_all_accel(boost::function<Real(unsigned int i)> f,
                  size_t max_i, Real tolerance)
{
    RealVector pTable;
    pTable.reserve(max_i);

    const Real p_0(f(0));
    if (p_0 == 0.0)
    {
        return 0.0;
    }

    pTable.push_back(p_0);

    RealVector::size_type i(1);
    for(;  i < max_i; ++i)
    {
        const Real p_i(f(i));
        pTable.push_back(p_i);
    }

    Real sum;
    Real error;
    gsl_sum_levin_utrunc_workspace* 
        workspace(gsl_sum_levin_utrunc_alloc(i));
    gsl_sum_levin_utrunc_accel(&pTable[0], pTable.size(), workspace, 
                                &sum, &error);
    if (fabs(error) >= fabs(sum * tolerance))
    {
        _log.error("series acceleration error: %.16g"
                  " (rel error: %.16g), terms_used = %d (%d given)",
                  fabs(error), fabs(error / sum),
                  workspace->terms_used, pTable.size());
        // TODO look into this crashing behaviour
    }

    gsl_sum_levin_utrunc_free(workspace);

    return sum;
}


Real 
funcSum(boost::function<Real(unsigned int i)> f, size_t max_i, Real tolerance)
// funcSum
// ==
// Will simply calculate the sum over a certain function f, until it converges
// (i.e. the sum > tolerance*current_term for a CONVERGENCE_CHECK number of 
// terms), or a maximum number of terms is summed (usually 2000).
//
// Input:
// - f: A       function object
// - max_i:     maximum number of terms it will evaluate
// - tolerance: convergence condition, default value 1-e8 (see .hpp)
// 
// About Boost::function
// ===
// Boost::function doesn't do any type checking: It will take any object 
// and any signature you provide in its template parameter, and create an 
// object that's callable according to your signature and calls the object. 
// If that's impossible, it's a compile error.
// (From: http://stackoverflow.com/questions/527413/how-boostfunction-and-boostbind-work)
{
    // DEFAULT = 4
    const unsigned int CONVERGENCE_CHECK(4);    

    Real sum(0.0);
    RealVector pTable;

    const Real p_0(f(0));
    if (p_0 == 0.0)
    {
        return 0.0;
    }

    pTable.push_back(p_0);
    sum = p_0;

    bool extrapolationNeeded(true);

    unsigned int convergenceCounter(0);

    RealVector::size_type i(1); 
    while(i < max_i)
    {
        const Real p_i(f(i));
        pTable.push_back(p_i);
        sum += p_i;

        ++i;

        if (fabs(sum) * tolerance >= fabs(p_i)) // '=' is important
        {
            ++convergenceCounter;          
        }
        // this screws it up; why?
        else
        {
            convergenceCounter = 0;
        }


        if (convergenceCounter >= CONVERGENCE_CHECK)
        {
            extrapolationNeeded = false;
            break;
        }
        
    }

    if (extrapolationNeeded)
    {
        Real error;
        gsl_sum_levin_utrunc_workspace* 
            workspace(gsl_sum_levin_utrunc_alloc(i));
        gsl_sum_levin_utrunc_accel(&pTable[0], pTable.size(), workspace, 
        &sum, &error);
        if (fabs(error) >= fabs(sum * tolerance * 10))
        {
            _log.error("series acceleration error: %.16g"
                      " (rel error: %.16g), terms_used = %d (%d given)",
                      fabs(error), fabs(error / sum),
                      workspace->terms_used, pTable.size());
        }

        gsl_sum_levin_utrunc_free(workspace);
    }

    return sum;
}
