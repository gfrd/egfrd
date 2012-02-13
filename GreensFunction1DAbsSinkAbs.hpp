#if !defined( __GREENSFUNCTION1DABSSINKABS_HPP )
#define __GREENSFUNCTION1DABSSINKABS_HPP

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
#include "GreensFunction.hpp"
#include "funcSum.hpp"
#include "freeFunctions.hpp"
#include "Logger.hpp"

class GreensFunction1DAbsSinkAbs: public GreensFunction
{
public:
    typedef std::pair<Real, Real> real_pair;
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
    GreensFunction1DAbsSinkAbs(Real D, Real k, Real r0, Real rsink, Real sigma, Real a)
	: GreensFunction(D), k(k), r0(r0), sigma(sigma), a(a), rsink(rsink), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
	    /* Set variables which define a domain with the sink at the origin. 
	       Futhermore r0 is assumed to be ringht from the sink. */
	    assert( a > sigma );
	    
	    L0 = fabs( r0 - rsink ); 
	    if( r0 >= rsink )
	    {
	        Lr = a - rsink; 
	        Ll = rsink - sigma;	        
	    }
	    else
	    {
	        Lr = rsink - sigma; 
	        Ll = a - rsink;
	    }

        calculate_n_roots( 1 );
    }

    ~GreensFunction1DAbsSinkAbs()
    {
	;   // empty
    }

    Real geta() const
    {
	    return a;
    }
    
    Real getLr() const
    {
	    return Lr;
    }
    
    Real getLl() const
    {
	    return Ll;
    }
    
    Real getL0() const
    {
	    return L0;
    }

    Real getrsink() const
    {
	    return rsink;
    }
    
    Real getsigma() const
    {
	    return sigma;
    }

    Real getr0() const
    {
	    return r0;
    }

    Real getk() const
    {
	    return k;
    }

    /* Calculates the probability density of finding the particle at 
       location z at timepoint t, given that the particle is still in the 
       domain. */
    Real calcpcum(Real r, Real t) const;

    /* Determine which event has occured at time t. Either an escape 
       (left or right boundary) or a reaction with the sink. 
       Based on the fluxes through the boundaries at the given time. */
    EventKind drawEventType( Real rnd, Real t ) const;

    /* Draws the first passage time from the propensity function */
    Real drawTime(Real rnd) const;

    /* Draws the position of the particle at a given time, assuming that 
       the particle is still in the domain. */
    Real drawR(Real rnd, Real t) const;

    /* Calculates the probability flux leaving the domain through the right 
       absorbing boundary at time t. */
    Real flux_leavea(Real t) const;

    /* Calculates the probability flux leaving the domain through the left 
       absorbing boundary at time t. */
    Real flux_leaves(Real t) const;

    /* Calculates the probability flux leaving the domain through the sink 
       at time t. */
    Real flux_sink(Real t) const;

    /* Calculates the probability of finding the particle inside the 
       domain at time t -> the survival probability */
    Real p_survival(Real t) const;

    /* c.d.f. of the greensfunction with respect to position. */
    Real p_int_r_table( Real const& r, Real const& t, RealVector& table ) const;

    Real p_int_r(Real const& r, Real const& t) const;

    /* TODO: Private methods which are public for now - debugging */

    /* Calculates the total probability flux leaving the domain at time t. */
    Real flux_tot(Real t) const;

    /* Calculates the probability flux leaving the domain through the 
       sub-domain containing r0 via the absorbing boundary and the flux 
       leaving the sub-domain not containing r0 via the absorbing boundary. */
    Real flux_abs_Lr(Real t, uint const& maxi) const;
    Real flux_abs_Ll(Real t, uint const& maxi) const;

    /* Calculates the probability density of finding the particle at 
       location r at time t. */
    Real prob_r(Real r, Real t) const;
    
    std::string dump() const;

    const char* getName() const
    {
        return "GreensFunction1DAbsSinkAbs";
    }

private:
    /* used structures */

    struct drawR_params
    {
        GreensFunction1DAbsSinkAbs const* gf;
        const Real t;
   	    RealVector& table;
	    const Real rnd;
    };

    struct drawT_params
    {
	    GreensFunction1DAbsSinkAbs const* gf;
	    RealVector& table;
	    const Real rnd;
    };

    struct root_f_params
    {
	    Real Lm_L;
	    Real h;
    };

    struct lower_upper_params
    {
        Real h;
        Real Lm_L;
        Real long_period;
        Real short_period;
        Real last_long_root;
        Real last_short_root;
        bool last_was_long;
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


    /* Functions concerned with finding the roots. */

    /* Fills the rootList with the first n roots */
    void calculate_n_roots( uint const& n ) const;
    
    /* Function returns two positions on the x-axis which straddle the next root. */
    real_pair get_lower_and_upper() const;

    /* Function of which we need the roots. */
    static Real root_f(Real x, void *p);

    /* Guess the number of terms needed for convergence, given t. */
    uint guess_maxi( Real const& t ) const;


    /* Function for calculating the survival probability. */
    Real p_survival_table(Real t, RealVector& psurvTable) const;

    Real p_survival_i( uint i, Real const& t, RealVector const& table ) const;

    Real p_survival_table_i( Real const& root_i ) const;
    
    void createPsurvTable( RealVector& table ) const;


    /* Functions for calculating the greensfunction. */
    Real prob_r_r0_i(uint i, Real const& rr, Real const& t) const;

    Real prob_r_nor0_i( uint i, Real const& rr, Real const& t) const;


    /* Functions for calculating the fluxes. */
    Real flux_tot_i(uint i, Real const& t) const;

    Real flux_abs_Lr_i(uint i, Real const& t) const;

    Real flux_abs_Ll_i(uint i, Real const& t) const;



    /* functions for calculating the c.d.f. */

    /* i'th term of p_int_r(r') for r' in left domain */
    Real p_int_r_leftdomain(uint i, Real const& rr, Real const& t, RealVector& table) const;

    /* i'th term of p_int_r(r') for r' in right domain, left of r0 */
    Real p_int_r_rightdomainA(uint i, Real const& rr, Real const& t, RealVector& table) const;

    /* i'th term of p_int_r(r') for r' in right domain, right of r0 */
    Real p_int_r_rightdomainB(uint i, Real const& rr, Real const& t, RealVector& table) const;

    /* Fills table with r-independent part of p_int_r_i. */
    void create_p_int_r_Table( Real const& t, RealVector& table ) const;

    /* Returns i'th r-independent term of p_int_r_i.
       Term is created if not in table. */
    Real get_p_int_r_Table_i( uint& i, Real const& t, RealVector& table) const
    {
        if( i >= table.size() )
        {
            calculate_n_roots( i+1 );
            create_p_int_r_Table( t, table );
        }

        return table[i];
    }

    /* Denominator of the Greens function */
    inline Real p_denominator_i( Real const& root_n ) const;
    
    /* Standard form of Greens Function: exp( -Dt root_n ** 2 ) / denominator */
    inline Real p_exp_den_i(Real const& t, Real const& root_n, Real const& root_n2) const;

    /* Function for drawR */
    static Real drawT_f (Real t, void *p);

    /* Function for drawTime */
    static Real drawR_f (Real r, void *p);

    
    /* Class variables */

    // The reaction constant
    const Real k;
    //starting position
    const Real r0;
    // The left and right boundary of the domain (sets the l_scale, see below)
    const Real sigma;
    const Real a;
    //Position of the sink in the domain.
    const Real rsink;
    // This is the length scale of the system
    Real l_scale;
    // This is the time scale of the system.
    Real t_scale;
    
    /* Greensfunction assumes that the sink is at the origin, and
       consists of two sub-domains: one between a boundary and the sink including
       r0, and one between boundary and sink not inlcuding r0. */

    // Length of sub-domain which does not include r0.
    Real Lr;
    // Length of sub-domain which does include r0.
    Real Ll;
    // Distance between the sink and r0.
    Real L0;

    // Stores all the roots.
    mutable RealVector rootList;

    // Stores params for rootfiner.
    mutable struct lower_upper_params lo_up_params;

    static Logger& log_;
};
#endif // __GREENSFUNCTION1DRADABS_HPP
