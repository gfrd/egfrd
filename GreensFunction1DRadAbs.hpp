#if !defined( __GREENSFUNCTION1DRADABS_HPP )
#define __GREENSFUNCTION1DRADABS_HPP

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
#include "funcSum.hpp" 
#include "freeFunctions.hpp"
#include "Defs.hpp"
#include "GreensFunction.hpp"
#include "Logger.hpp"

class GreensFunction1DRadAbs: public GreensFunction
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
    // Is 1E3 a good measure for the probability density?!
    static const Real PDENS_TYPICAL = 1;
    // The maximum number of terms used in calculating the sum
    static const int MAX_TERMS = 500;
    // The minimum number of terms
    static const int MIN_TERMS = 20;
    /* Cutoff distance: When H * sqrt(2Dt) < a - r0 OR ro - sigma
       use free greensfunction instead of absorbing. */
    static const Real CUTOFF_H = 6.0;

public:
    GreensFunction1DRadAbs(Real D, Real k, Real r0, Real sigma, Real a)
	: GreensFunction(D), v(0.0), k(k), r0(r0), sigma(sigma), a(a), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
        //set first root.
        calculate_n_roots( 1 );
    }

    // The constructor is overloaded and can be called with or without drift v
    // copy constructor including drift variable v
    GreensFunction1DRadAbs(Real D, Real v, Real k, Real r0, Real sigma, Real a)
	: GreensFunction(D), v(v), k(k), r0(r0), sigma(sigma), a(a), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
        //set first root;
        calculate_n_roots( 1 );
    }

    ~GreensFunction1DRadAbs()
    {
	;   // empty
    }

    // This also sets the scale
    void seta(Real a)
    {
	THROW_UNLESS( std::invalid_argument, (a-this->sigma) >= 0.0 && this->r0 <= a);

	// Use a typical domain size to determine if we are here 
	// defining a domain of size 0.
	if ( (a-this->sigma) < EPSILON*this->l_scale )
	{
	    // just some random value to show that the domain is zero
	    this->a = -1.0;
	}
	else
	{
	    // set the l_scale to the given one
	    this->l_scale = a-sigma;
	    // set the typical time scale (MSD = sqrt(2*d*D*t) )
	    this->t_scale = (l_scale*l_scale)/this->getD();
	    this->a = a;
	}
    }

    Real geta() const
    {
        return this->a;
    }
    
    Real getsigma() const
    {
        return this->sigma;
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
        return r0;
    }

    Real getk() const
    {
        return this->k;
    }

    Real getv() const
    {
        return this->v;
    }

    // Calculates the probability density of finding the particle at 
    // location z at timepoint t, given that the particle is still in the 
    // domain.
    Real calcpcum (Real r, Real t) const;

    // Determine which event has occured, an escape or a reaction. Based 
    // on the fluxes through the boundaries at the given time. Beware: if 
    // t is not a first passage time you still get an answer!
    EventKind drawEventType( Real rnd, Real t ) const;

    // Draws the first passage time from the propensity function
    Real drawTime (Real rnd) const;

    // Draws the position of the particle at a given time, assuming that 
    // the particle is still in the domain
    Real drawR (Real rnd, Real t) const;


// These methods are both public and private, they are used by public methods 
// but can also be called from the 'outside'. This is mainly because of 
// debugging purposes.


    // Calculates the probability of finding the particle inside the 
    // domain at time t -> the survival probability
    Real p_survival (Real t) const;

    // Calculates the total probability flux leaving the domain at time t
    Real flux_tot (Real t) const;

    // Calculates the probability flux leaving the domain through the 
    // radiative boundary at time t
    Real flux_rad (Real t) const;

    // Calculates the flux leaving the domain through the radiative 
    // boundary as a fraction of the total flux. This is the probability 
    // that the particle left the domain through the radiative
    // boundary instead of the absorbing boundary.
    Real fluxRatioRadTot (Real t) const;

    // Calculates the probability density of finding the particle at 
    // location r at time t.
    Real prob_r (Real r, Real t) const;
    
// End of public/private mix methods

//private:	// method made public for testing

    std::string dump() const;

    const char* getName() const
    {
        return "GreensFunction1DRadAbs";
    }
    
private:

    Real An (Real a_n) const;

    Real Bn (Real a_n) const;

    Real Cn (Real a_n, Real t) const;

    struct tan_f_params
    {
        Real a;
        Real h;
    };

    struct drawT_params
    {
        GreensFunction1DRadAbs const* gf;
        RealVector& psurvTable;
        Real rnd;
    };

   struct drawR_params
    {
        GreensFunction1DRadAbs const* gf;
        const Real t;
        RealVector table;
        Real rnd;
    };


    /* Functions managing the rootList */
    
    /* return the rootList size */
    uint rootList_size() const
    {
        return rootList.size();
    }

    /* returns the last root. */
    Real get_last_root() const
    {
        return rootList.back();
    }

    /* ad a root to the rootList */
    void ad_to_rootList( Real const& root_i ) const
    {
        rootList.push_back( root_i );
    }

    /* remove n'th root from rootList */
    void remove_from_rootList(uint const& n) const
    {
        rootList.erase( rootList.begin() + n );
    }

    /* return the n + 1'th root */
    Real get_root( uint const& n ) const
    {
        if( n >= rootList.size() )
            calculate_n_roots( n+1 );

        return rootList[ n ];
    }

    /* Fills the rootList with the first n roots. */
    void calculate_n_roots(uint const& n) const;

    /* Guess the number of terms needed for convergence, given t. */
    uint guess_maxi( Real const& t ) const;

    /* this is the appropriate definition of the function in gsl. */
    static double tan_f (double x, void *p);


    /* functions for drawTime / p_survival */

    static double drawT_f (double t, void *p);

    Real p_survival_table( Real  t, RealVector& psurvTable ) const;

    Real p_survival_i(uint i, Real const& t, RealVector const& table ) const;

    Real p_survival_table_i_v( uint const& i ) const;

    Real p_survival_table_i_nov( uint const& i ) const;

    void createPsurvTable( RealVector& table) const;


    /* functions for drawR */

    static Real drawR_f (Real z, void* p);

    Real p_int_r_table(Real const& r, Real const& t, RealVector& table) const;

    Real p_int_r_i(uint i, Real const& r, Real const& t, RealVector& table) const;

    void create_p_int_r_Table( Real const& t, RealVector& table ) const;
    
    Real get_p_int_r_Table_i( uint& i, Real const& t, RealVector& table) const
    {
        if( i >= table.size() )
        {
            calculate_n_roots( i+1 );
            create_p_int_r_Table(t, table);
        }

        return table[i];
    }

    /* Member variables */

    // The diffusion constant and drift velocity
    Real v;
    // The reaction constant
    Real k;
    Real r0;
    // The left and right boundary of the domain (sets the l_scale, see below)
    Real sigma;
    Real a;
    // This is the length scale of the system
    Real l_scale;
    // This is the time scale of the system.
    Real t_scale;

    /* vector containing the roots 0f tan_f. */
    mutable RealVector rootList;

    static Logger& log_;
};
#endif // __GREENSFUNCTION1DRADABS_HPP
