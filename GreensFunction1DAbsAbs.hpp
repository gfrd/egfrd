#if !defined( __GREENSFUNCTION1DABSABS_HPP )
#define __GREENSFUNCTION1DABSABS_HPP

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_roots.h>

#include <math.h>

#include "findRoot.hpp"
#include "Defs.hpp"
#include "OldDefs.hpp"			// TODO: this must be removed at some point!
#include "funcSum.hpp"
#include "freeFunctions.hpp"
#include "Logger.hpp"

#include "GreensFunction.hpp"
#include "PairGreensFunction.hpp"	// needed to declare EventType


class GreensFunction1DAbsAbs: public GreensFunction
{
public:
    typedef std::vector<Real> RealVector;
    typedef unsigned int uint;

private:
    // This is a typical length scale of the system, may not be true!
    static const Real L_TYPICAL = 1E-8;
    // The typical timescale of the system, may also not be true!!
    static const Real T_TYPICAL = 1E-6;
    // measure of 'sameness' when comparing floating points numbers
    static const Real EPSILON = 1E-10;
    //E3; Is 1E3 a good measure for the probability density?!
    static const Real PDENS_TYPICAL = 1;
    // The maximum number of terms in the sum
    static const uint MAX_TERMS = 500;
    // The minimum
    static const uint MIN_TERMS = 20;
    // Cutoff distance: When H * sqrt(2Dt) < 1/2*L, use free greensfunction instead of absorbing.
    static const Real CUTOFF_H = 6.0;


public:
    GreensFunction1DAbsAbs(Real D, Real r0, Real sigma, Real a)
	: GreensFunction(D), v(0.0), sigma(sigma), a(a), r0(r0), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
	;   // do nothing
    }

    // The constructor is overloaded and can be called with or without drift v
    GreensFunction1DAbsAbs(Real D, Real v, Real r0, Real sigma, Real a) // copy constructor including drift variable v
	: GreensFunction(D), v(v), sigma(sigma), a(a), r0(r0), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
	;   // do nothing
    }

    ~GreensFunction1DAbsAbs()
    { 
	;   // empty
    }

    // This also sets the scale
    void seta(Real a)
    {
	Real L( a - this->sigma );
	
	THROW_UNLESS( std::invalid_argument, L >= 0.0 && (this->r0 - sigma) <= L);

	// Use a typical domain size to determine if we are here
	// defining a domain of size 0.
	if ( L <= EPSILON * l_scale )
	{
	    // just some random value to show that the domain is zero
	    this->a = -INT_MAX;
	}
	else
	{   
	    // set the typical time scale (msd = sqrt(2*d*D*t) )
	    // this is needed by drawTime_f, do not get rid of it!
	    this->t_scale = (L*L)/this->getD();
	    // set a
	    this->a = a;
	}
    }

    Real getsigma() const
    {
	return this->sigma;
    }
 
    Real geta() const
    {
	return this->a;
    }

    Real getv() const
    {
	return this->v;
    }

    void setr0(Real r0)
    {
	if ( this->a - this->sigma < 0.0 )
	{
	    // if the domain had zero size    
	    THROW_UNLESS( std::invalid_argument,
	                  0.0 <= (r0-sigma) && (r0-sigma) <= EPSILON * l_scale );
	    this->r0 = 0.0;
	}
	else
	{
	    // The normal case
	    THROW_UNLESS( std::invalid_argument,
	                  0.0 <= (r0-sigma) && r0 <= this->a);
	    this->r0 = r0;
	}
    }

    Real getr0() const
    {
	return this->r0;
    }

    // Draws the first passage time from the propensity function
    Real drawTime (Real rnd) const;

    // Draws the position of the particle at a given time, assuming that 
    // the particle is still in the
    // domain
    Real drawR (Real rnd, Real t) const;

    // Calculates the amount of flux leaving the left boundary at time t
    Real leaves(Real t) const;

    // Calculates the amount of flux leaving the right boundary at time t
    Real leavea(Real t) const;

    // Determines based on the flux ratios if the particle left the left 
    // or right boundary
    EventKind drawEventType( Real rnd, Real t ) const;

    // Calculates the probability of finding the particle inside the 
    // domain at time t so, the survival probability
    Real p_survival (Real t) const;

    // Calculates the probability density of finding the particle at 
    // location z at timepoint t, given that the particle is still in the 
    // domain.
    Real calcpcum (Real r, Real t) const;
    
    // Calculates the probability density of finding the particle at 
    // location r at time t.
    Real prob_r (Real r, Real t) const;

    std::string dump() const;

    const char* getName() const
    {
        return "GreensFunction1DAbsAbs";
    }

private:
    struct drawT_params
    {
        GreensFunction1DAbsAbs const* gf;
        RealVector& psurvTable;
        Real rnd;
    };

   struct drawR_params
    {
        GreensFunction1DAbsAbs const* gf;
        const Real t;
        RealVector table;
        Real rnd;
        Real p_surv;
    };

    uint guess_maxi(Real const& t ) const;


    static Real drawT_f (Real t, void *p);

    Real p_survival_table( Real  t, RealVector& psurvTable ) const;

    Real p_survival_i(uint i, Real const& t, RealVector const& table ) const;

    Real p_survival_table_i_v( uint const& i ) const;

    Real p_survival_table_i_nov( uint const& i ) const;

    void createPsurvTable( uint const& maxi, RealVector& table) const;


    static Real drawR_f (Real z, void* p);

    Real p_int_r_table(Real const& r, Real const& t, Real const& p_surv, 
                       RealVector& table) const;

    void createP_int_r_Table( Real const& t, uint const& maxi, RealVector& table ) const;

private:
    // The diffusion constant and drift velocity
    Real v;
    // These are the dimensions of our domain; L is calculated as a-sigma
    Real sigma;
    Real a;
    Real r0;
    // This is the 'length scale' of your system (1e-14 or 1e6)
    // Although rescaling is discontinued, we use it to check whether a is well-chosen
    Real l_scale;
    // This is the time scale of the system, used by drawTime_f
    Real t_scale;

    static Logger& log_;
};
#endif // __GREENSFUNCTION1DABSABS_HPP
