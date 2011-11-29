#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>

#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_roots.h>

#include "findRoot.hpp"
#include "funcSum.hpp"
#include "GreensFunction1DAbsSinkAbs.hpp"

/* This is the appropriate definition of the function defining
   the roots of our Green's functions in GSL.
   Later needed by the rootfinder. */
Real GreensFunction1DAbsSinkAbs::root_f (Real x, void *p)
{
    struct root_f_params *params = (struct root_f_params *)p;
    const Real Lm_L = (params->Lm_L);
    const Real h = (params->h);

    // L   = Lr + Ll
    // h    = L * k / (2 * D)
    // L_Lm = Lr + Ll / Lr - Ll
    // x    = q * L

    return x * sin(x) + h * ( cos(x * Lm_L) - cos(x) );

}

/* Calculates the first n roots of root_f */
void GreensFunction1DAbsSinkAbs::calculate_n_roots(uint const& n) const
{
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L( Lr + Ll );
    const Real Lm( Lr - Ll );
    const Real Lm_L( Lm / L );
    const Real h( getk() * L / ( 2 * getD() ) );
    
    Real root_i;
    uint i = 0;

    real_pair lower_upper_pair;

    lo_up_params.h = h;
    lo_up_params.Lm_L = Lm_L;
    lo_up_params.long_period = std::max( L/Lr * M_PI, L/Ll * M_PI );
    lo_up_params.short_period = std::min( L/Lr * M_PI, L/Ll * M_PI );

    /* Define the root function. */
    gsl_function F;
    struct root_f_params params = { Lm_L, h };
     
    F.function = &GreensFunction1DAbsSinkAbs::root_f;
    F.params = &params;

    /* define and ad a new solver type brent */
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    i = rootList_size();

    /* Stores the last root with long or short period. */
    if(i == 0)
    {
        lo_up_params.last_long_root = 0.0;
        lo_up_params.last_short_root = 0.0;
    }

    /* Find all the roots up to the nth */
    while(i++ < n)
    {
        lower_upper_pair = get_lower_and_upper();

        root_i = findRoot( F, solver, lower_upper_pair.first, lower_upper_pair.second, 
                           1.0*EPSILON, EPSILON, "GreensFunction1DAbsSinkAbs::root_f" );

        assert( root_i > std::max(lo_up_params.last_long_root, 
                                  lo_up_params.last_short_root) 
                - EPSILON );

        ad_to_rootList( root_i / L );

        if(lo_up_params.last_was_long)
            lo_up_params.last_long_root = root_i;
        else
            lo_up_params.last_short_root = root_i;

        std::cerr << "root#: " << rootList_size() << " - " << get_last_root() << std::endl;
    }

    gsl_root_fsolver_free( solver );
}


/* returns two points on the x-axis which straddle the next root. */
std::pair<Real, Real> GreensFunction1DAbsSinkAbs::get_lower_and_upper() const
{
    const Real root_n( std::max(lo_up_params.last_long_root, 
                                lo_up_params.last_short_root) );
    const Real safety( .75 );

    Real lower, upper, next_root_est, left_offset, right_offset;

    const Real last_root( root_n ==  0.0 ? M_PI : root_n );

    if( lo_up_params.h / last_root < 1 )
    {
        right_offset = M_PI;
        next_root_est = root_n + M_PI;
    }
    else
    {
        const Real next_root_long( lo_up_params.last_long_root 
                                   + lo_up_params.long_period );
        const Real next_root_short( lo_up_params.last_short_root
                                    + lo_up_params.short_period );

        if( next_root_long < next_root_short )
        {
            next_root_est = next_root_long;

            right_offset = std::min( next_root_short - next_root_est, 
                              lo_up_params.long_period );

            lo_up_params.last_was_long = true;
        }
        else
        {
            next_root_est = next_root_short;
            
            right_offset = std::min( next_root_long - next_root_est, 
                              lo_up_params.short_period );

            lo_up_params.last_was_long = false;                        
        }
    }
    
    left_offset = next_root_est - root_n - 1000 * EPSILON;      

    lower = next_root_est - left_offset;
    upper = next_root_est + safety * right_offset;

    struct root_f_params p = { lo_up_params.Lm_L, lo_up_params.h };

    Real f_lower( root_f( lower, &p ) );
    Real f_upper( root_f( upper, &p ) );

    /* set the parity operator for the next root: 
       +1 for an even root.
       -1 for an odd root. */
    const int parity_op( 2 * ( rootList_size()%2 ) - 1 );

    /* f_lower must have correct sign. */
    if( f_lower * parity_op > 0 )
    {
        std::cerr << "Parity error in lower of root# " << rootList_size() + 1 
                      << " in GF::AbsSinkAbs " << std::endl;
        std::cerr << "f_low(" << lower << ") = " << f_lower << ", f_high(" 
                  << upper << ") = " << f_upper << std::endl;
    }

    /* If the parity is incorrect -> correct it */
    if( f_upper * parity_op < 0)
    {
        int cntr = 0;

        const Real delta( .1 * std::min(left_offset, right_offset) );
        
        cntr = 0;
        //TODO: find out is upper is left or right from next_root.
        while(f_upper * parity_op < 0 && cntr++ < 10 )
        {
            //Assume overshoot for now.
            upper -= delta;
            f_upper = root_f( upper, &p );
        }

        if(cntr >= 10)
        {
            std::cerr << "Failed to straddle root # " << rootList_size() + 1 
                      << " in GF::AbsSinkAbs " << std::endl;
            std::cerr << "f_low(" << lower << ") = " << f_lower << ", f_high(" 
                  << upper << ") = " << f_upper << std::endl;
        }

    }

    return real_pair( lower, upper );
}

/* returns a guess for the number of terms needed for 
   the greensfunction to converge at time t */
uint GreensFunction1DAbsSinkAbs::guess_maxi(Real const& t) const
{
    const uint safety(2);

    if (t >= INFINITY)
    {
        return safety;
    }

    const Real D( getD() );

    const Real root0( get_root( 0 ) );
    const Real Dt(D * t);

    const Real thr(exp(- Dt * root0 * root0) * EPSILON * 1e-1);

    if (thr <= 0.0)
    {
        return MAX_TERMS;
    }

    const Real max_root( sqrt(root0 * root0 - log(thr) / Dt) );

    const uint maxi(safety + 
          static_cast<unsigned int>(max_root * ( getLr() + getLl() )  / M_PI));

    return std::min(maxi, MAX_TERMS);
}


/* Standart form of the greensfunction without numerator */
inline Real GreensFunction1DAbsSinkAbs::p_exp_den_i(Real const& t, 
                                                  Real const& root_i, 
                                                  Real const& root_i2) const
{
    return exp( - getD() * root_i2 * t ) / p_denominator_i( root_i );
}


/* Denominator of the greensfunction. */
inline Real GreensFunction1DAbsSinkAbs::p_denominator_i(Real const& root_n) const
{
    const Real Lm( getLr() - getLl() );
    const Real L( getLr() + getLl() );
    
    const Real term1( root_n * L * cos( root_n * L ) + sin( root_n * L ) );
    const Real term2( L * sin( root_n * L ) - Lm * sin( root_n * Lm ) );   

    return getD() * term1 + getk() / 2. * term2;
}


Real GreensFunction1DAbsSinkAbs::p_survival(Real t) const
{
    RealVector table;
    return p_survival_table(t, table);
}

/* Calculates survival probability using a table. 
   Switchbox for which greensfunction to use. */
Real GreensFunction1DAbsSinkAbs::p_survival_table(Real t, RealVector& psurvTable) const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    Real p;

    const Real D( getD() );
    //const Real sigma(getsigma());
    //const Real a( geta() );

    //    const Real L0( getL0() );
    //    const Real distToa(a - r0);
    //    const Real distTos(r0 - sigma);

    //    const Real H(6.0); // a fairly strict criterion for safety.
    //    const Real maxDist(H * sqrt(2.0 * D * t));

    if (t == 0.0 || (D == 0.0 && v == 0.0) )
    {
	    //particle can't escape.
	    return 1.0;
    }

    /*
    if (distToa > maxDist)
    {
        if (distTos > maxDist) // far from anything; it'll survive.
        {
            if( L0 > maxDist )
                p = 1.0;
        }
        else // close only to s, ignore a
        {
            const Real sigma(this->getSigma());
            const Real kf(this->getkf());
            p = p_survival_irr(t, r0, kf, D, sigma);
        }
    }
    else
    {
        if (distTos > maxDist)  // close only to a.
        {
            p = p_survival_nocollision(t, r0, D, a);
        }
        else  // close to both boundaries.  do the normal calculation.
        {
            const unsigned int maxi(guess_maxi(t));
            
            if (psurvTable.size() < maxi + 1)
            {
                IGNORE_RETURN getAlpha0(maxi);  // this updates the table
                this->createPsurvTable(psurvTable);
            }

            p = funcSum_all(boost::bind(&GreensFunction3DRadAbs::
                                          p_survival_i_exp_table, 
                                          this,
                                          _1, t, psurvTable),
                             maxi);
        }
    }
    */

    const uint maxi( guess_maxi(t) );
            
    if (psurvTable.size() < maxi + 1)
    {
        calculate_n_roots( maxi );  // this updates the table
        createPsurvTable( psurvTable );
    }
    
    p = funcSum_all(boost::bind(&GreensFunction1DAbsSinkAbs::p_survival_i, 
                                this, _1, t, psurvTable),
                    maxi);
    
    return p;
}


/* Calculates the i'th term of the p_survival sum. */
Real GreensFunction1DAbsSinkAbs::p_survival_i( uint i, 
                                               Real const& t, 
                                               RealVector const& table) const
{
    const Real root_i( get_root( i ) );
    return exp( - getD() * t * gsl_pow_2( root_i ) ) * table[ i ];
}


/* Calculates the part of the i'th term of p_surv not dependent on t */
Real GreensFunction1DAbsSinkAbs::p_survival_table_i( Real const& root_i ) const
{ 
    const Real D ( getD() );
    const Real k ( getk() ); 
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real L ( Lr + Ll );
    const Real LrmL0( Lr - L0 );

    const Real term1( sin( root_i * L ) - 
                      sin( root_i * LrmL0 ) - sin( root_i * (Ll + L0) ) );

    const Real term2( sin( root_i * Lr ) - sin( root_i * L0 ) 
                      - sin( root_i * LrmL0 ) );
        
    Real numerator( D * term1 + k * sin( root_i * Ll ) * term2 / root_i ); 
    numerator *= 2.0;

    return numerator / p_denominator_i( root_i );
}


/* Fills table with terms in the p_survival sum which don't depend on t. */
void GreensFunction1DAbsSinkAbs::createPsurvTable(RealVector& table) const
{
    uint const root_nbr( rootList_size() );
    uint i( table.size() );

    while( i <= root_nbr )
    {
        table.push_back( p_survival_table_i( get_root( i++ ) ) );
    }
}


/* Returns i'th term of prob_r (domain containing r0) function */
Real GreensFunction1DAbsSinkAbs::prob_r_r0_i(uint i, 
                                             Real const& rr, 
                                             Real const& t) const
{
    const Real root_i( get_root(i) );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    Real L0( getL0() );
    Real rr2( rr );
    
    /* Check if r is left or right of the starting position r0.
       If so, interchange rr with L0. */
    if( rr < L0 )
    {
        rr2 = L0;
        L0 = rr;
    }

    const Real LlpL0( Ll + L0 );
    const Real Lrmrr( Lr - rr2 );
    
    Real numerator( getD() * root_i * sin( root_i * LlpL0 ) + 
                    getk() * sin( root_i * Ll ) * sin( root_i * L0) );

    numerator *= sin( root_i * Lrmrr );
    
    return - 2.0 * p_exp_den_i(t, root_i, gsl_pow_2( root_i ) ) * numerator;
}


/* Returns i'th term of prob_r (domain not containing r0) function */
Real GreensFunction1DAbsSinkAbs::prob_r_nor0_i(uint i, 
                                               Real const& rr, 
                                               Real const&t) const
{
    const Real root_i( get_root(i) );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );

    const Real LrmL0( Lr - L0 );
    const Real Llprr( Ll + rr );

    const Real numerator( getD() * root_i * sin( root_i * Llprr ) * sin( root_i * LrmL0 ) );

    return - 2.0 * p_exp_den_i(t, root_i, gsl_pow_2( root_i ) ) * numerator;
}


/* Calculates the probability density of finding the particle at location r
   at time t. */
Real GreensFunction1DAbsSinkAbs::prob_r(Real r, Real t) const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    THROW_UNLESS( std::invalid_argument, (r-sigma) >= 0.0 && r <= a && (r0 - sigma) >= 0.0 && r0<=a );
    
    const Real D( getD() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L( Lr + Ll );    

    // if there was no time change or zero diffusivity => no movement
    if (t == 0 || D == 0)
    {
	    // the probability density function is a delta function
	    if (r == r0)
	    {
	        return INFINITY;
	    }
	    else
	    {
	        return 0.0;
	    }
    }

    // if r is at one of the the absorbing boundaries
    if ( fabs( a - r ) < EPSILON * L || fabs( r - sigma ) < EPSILON * L  )
    {
	    return 0.0;
    }

    const Real rr( getr0() - getrsink() >= 0 ? r - rsink : rsink - r  );

    Real p;  

    const uint maxi( guess_maxi( t ) );
    calculate_n_roots( maxi );

    /* Determine wether rr lies in the same sub-domain as r0.
       A different function is caculated when this is the case. */
    if( rr >= 0 )
    {   
        p = funcSum_all(boost::bind(&GreensFunction1DAbsSinkAbs::prob_r_r0_i, 
                                    this, _1, rr, t),
                        maxi);
	}
    else
    {
        p = funcSum_all(boost::bind(&GreensFunction1DAbsSinkAbs::prob_r_nor0_i, 
                                    this, _1, rr, t),
                        maxi);
    }
    
    return p;
}


/* Calculates the probability density of finding the particle at location r at
   time t, given that the particle is still in the domain. */
Real GreensFunction1DAbsSinkAbs::calcpcum(Real r, Real t) const
{
    return prob_r(r, t)/p_survival( t );
}


/* Function returns flux at absorbing bounday sigma. */
Real GreensFunction1DAbsSinkAbs::flux_leaves(Real t) const
{
    if( t == 0 || D == 0 )
    {
        return 0.0;
    }

    const Real maxi( guess_maxi( t ) );

    if( getr0() >= getrsink() )
        return flux_abs_Ll( t, maxi );
    else
        return - flux_abs_Lr( t, maxi );
}


/* Function returns flux at absorbing bounday a. */
Real GreensFunction1DAbsSinkAbs::flux_leavea(Real t) const
{
    if( t == 0 || D == 0 )
    {
        return 0.0;
    }

    const Real maxi( guess_maxi( t ) );

    if( getr0() < getrsink() )
        return - flux_abs_Ll( t, maxi );
    else
        return flux_abs_Lr( t, maxi );
}


/* Calculates the total probability flux leaving the domain at time t
   This is simply the negative of the time derivative of the survival prob.
   at time t [-dS(t')/dt' for t'=t]. */
Real GreensFunction1DAbsSinkAbs::flux_tot(Real t) const
{
    const Real D( getD() );
    Real p;

    const Real maxi( guess_maxi( t ) );

    p = funcSum_all(boost::bind(&GreensFunction1DAbsSinkAbs::flux_tot_i, 
                                this, _1, t),
                        maxi);

    return D * p;
}


/* Calculates i'th term of total flux leaving at time t. */
Real GreensFunction1DAbsSinkAbs::flux_tot_i(uint i, Real const& t) const
{
    const Real root_i( get_root( i ) );
    const Real root_i2( gsl_pow_2( root_i ) );

    return root_i2 * exp( - getD() * t * root_i2 ) * 
        p_survival_table_i( root_i );
}


/* Flux leaving throught absorbing boundary from sub-domain containing r0. */
Real GreensFunction1DAbsSinkAbs::flux_abs_Lr(Real t, uint const& maxi) const
{
    const Real D( getD() );
    Real p;
   
    p = funcSum_all(boost::bind(&GreensFunction1DAbsSinkAbs::flux_abs_Lr_i, 
                                this, _1, t),
                        maxi);

    return - D * 2 * p;
}


/* Calculates the i'th term of the flux at Lr. */
Real GreensFunction1DAbsSinkAbs::flux_abs_Lr_i(uint i, Real const& t) const
{
    const Real root_i( get_root( i ) );
    const Real D( getD() );
    const Real k( getk() );    
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real LlpL0( Ll + L0 );
        
    Real numerator( k * sin( root_i * Ll ) * sin ( root_i * L0 ) + 
                    D * root_i * sin( root_i * LlpL0 ) );

    numerator *= root_i;

    return p_exp_den_i(t, root_i, gsl_pow_2( root_i ) ) * numerator; 
}


/* Flux leaving throught absorbing boundary from sub-domain not containing r0. */
Real GreensFunction1DAbsSinkAbs::flux_abs_Ll(Real t, uint const& maxi) const
{
    const Real D2( gsl_pow_2( getD() ) );
    Real p;

    p = funcSum_all(boost::bind(&GreensFunction1DAbsSinkAbs::flux_abs_Ll_i, 
                                this, _1, t),
                        maxi);
    
    return 2 * D2 * p;
}


/* Calculates the i'th term of the flux at Ll. */
Real GreensFunction1DAbsSinkAbs::flux_abs_Ll_i(uint i, Real const& t) const
{
    const Real root_i( get_root( i ) );
    const Real root_i2( gsl_pow_2( root_i ) );
    const Real LrmL0( getLr() - getL0() );
        
    Real numerator( root_i2 * sin( root_i * LrmL0 ) );

    return p_exp_den_i(t, root_i, root_i2 ) * numerator;
}


/* Calculates the probability flux leaving the domain through the sink
   at time t */
Real GreensFunction1DAbsSinkAbs::flux_sink(Real t) const
{
    return getk() * prob_r(getrsink(), t);
}


/* Determine which event has occured, an escape or a reaction. Based on the
   fluxes through the boundaries and the sink at the given time. */
GreensFunction1DAbsSinkAbs::EventKind 
GreensFunction1DAbsSinkAbs::drawEventType( Real rnd, Real t ) const
{
    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, t > 0.0 );

    const Real a( geta() );
    const Real sigma( getsigma() );
    const Real r0( getr0() );
    const Real L( a - sigma );

    /* if the sink is impermeable (k==0) or
       the particle is at one the absorbing boundary (sigma or a) => IV_ESCAPE event */
    if ( k == 0 || fabs(a - r0) < EPSILON * L || fabs(sigma - r0) < EPSILON * L)
    {
	    return IV_ESCAPE;
    }

    /* The event is sampled from the flux ratios.
       Two possiblities:
       (1) Leave through left or right boundary - IV_ESCAPE
       (2) Leave through sink - IV_REACTION
    */

    rnd *= flux_tot( t );

    const Real p_sink( flux_sink( t ) );
    if (rnd < p_sink)
    {
      return IV_REACTION;
    }
    else
    {
      return IV_ESCAPE;
    }
}


/* This function is needed to cast the math. form of the function
   into the form needed by the GSL root solver. */
Real GreensFunction1DAbsSinkAbs::drawT_f(Real t, void *p)
{
    struct drawT_params *params = (struct drawT_params *)p;
    return params->rnd - params->gf->p_survival_table( t, params->table );
}


/* Draws the first passage time from the survival probability,
   using an assistance function drawT_f that casts the math. function
   into the form needed by the GSL root solver. */
Real GreensFunction1DAbsSinkAbs::drawTime(Real rnd) const
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
  
    const Real a( geta() );    
    const Real r0( getr0() );
    const Real D( getD() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real L( getLr() + getLl() );

    if ( D == 0.0 || L == INFINITY )
    {
	    return INFINITY;
    }

    if ( rnd <= EPSILON || L < 0.0 || fabs( a - r0 ) < EPSILON * L )
    {
	    return 0.0;
    }

    /* the structure to store the numbers to calculate the numbers for 1-S */
    RealVector psurvTable;
    drawT_params params = { this, psurvTable, rnd };

    /* Find a good interval to determine the first passage time.
       First we get the distance to one of the absorbing boundaries or the sink. */
    Real t_guess;
    Real dist( std::min(Lr - L0, Ll + L0) );
    dist = std::min( dist, L0 );

    t_guess = dist * dist / ( 2.0 * D );
    t_guess *= .1;

    const uint maxi( guess_maxi( t_guess ) );
    calculate_n_roots( maxi );

    /* Define the function for the rootfinder */
    gsl_function F;
    F.function = &GreensFunction1DAbsSinkAbs::drawT_f;
    F.params = &params;

    Real value( GSL_FN_EVAL( &F, t_guess ) );
    Real low( t_guess );
    Real high( t_guess );

    // scale the interval around the guess such that the function straddles
    if( value < 0.0 )
    {
        // if the guess was too low
        do
        {
            // keep increasing the upper boundary until the
            // function straddles
            high *= 10;
            value = GSL_FN_EVAL( &F, high );
            
            if( fabs( high ) >= t_guess * 1e6 )
            {
                std::cerr << "GF1DSink::drawTime Couldn't adjust high. F("
                          << high << ") = " << value << std::endl;
                throw std::exception();
            }
        }
        while ( value <= 0.0 );
    }
    else
    {
        // if the guess was too high
	    // initialize with 2 so the test below survives the first
	    // iteration
	    Real value_prev( 2 );
	    do
	    {
	        if( fabs( low ) <= t_guess * 1e-6 ||
	            fabs(value-value_prev) < EPSILON*1.0 )
	        {
		    std::cerr << "GF1DSink::drawTime Couldn't adjust low. F(" << low << ") = "
		              << value << " t_guess: " << t_guess << " diff: "
		              << (value - value_prev) << " value: " << value
		              << " value_prev: " << value_prev << " rnd: "
		              << rnd << std::endl;
		    return low;
	        }
	        value_prev = value;
	        // keep decreasing the lower boundary until the function straddles
	        low *= 0.1;
	        // get the accompanying value
	        value = GSL_FN_EVAL( &F, low );
	    }
	    while ( value >= 0.0 );
    }

    /* find the intersection on the y-axis between the random number and
       the function */
    
    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    // make a new solver instance
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    const Real t( findRoot( F, solver, low, high, t_scale*EPSILON, EPSILON,
                            "GreensFunction1DAbsSinkAbs::drawTime" ) );
    // return the drawn time
    return t;
}


/* Retrurns the c.d.f. with respect to the position at time t. 
   - Also a switchbox for which GF to choose. */
Real GreensFunction1DAbsSinkAbs::p_int_r_table(Real const& r, 
                                               Real const& t, 
                                               RealVector& table) const
{
    Real p;
    const Real rsink( getrsink() );

    /* when r0 lies left of the sink, mirror the domain around rsink
       : rr -> -rr. */
    const Real rr( getr0() - rsink >= 0 ? r - rsink : rsink - r  );

    const uint maxi( guess_maxi(t) );    
   
    if (table.size() < maxi + 1)
    {
        calculate_n_roots( maxi );  // this updates the table
        createP_int_r_Table( t, table );
    }

    /* Determine in which part of the domain rr lies, and
       thus which function to use. */
    Real (GreensFunction1DAbsSinkAbs::*p_int_r_i)
    (uint, Real const&, RealVector const& table) const = NULL;

    if( rr <= 0 )
        p_int_r_i = &GreensFunction1DAbsSinkAbs::p_int_r_leftdomain;
    else if( rr < getL0() )
            p_int_r_i = &GreensFunction1DAbsSinkAbs::p_int_r_rightdomainA;
        else
            p_int_r_i = &GreensFunction1DAbsSinkAbs::p_int_r_rightdomainB;

    p = funcSum_all(boost::bind(p_int_r_i, 
                                this, _1, rr, table),
                    maxi);

    return 2.0 * p;
}


Real GreensFunction1DAbsSinkAbs::p_int_r(Real const& r, 
                                         Real const& t) const
{
    THROW_UNLESS( std::invalid_argument, r >= getsigma() && r <= geta() );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    RealVector table;
    return p_int_r_table(r, t, table) / p_survival( t );
}


void GreensFunction1DAbsSinkAbs::createP_int_r_Table( Real const& t, 
                                                      RealVector& table ) const
{
    const uint root_nbr( rootList_size() );
    uint i( table.size() );

    while( i < root_nbr )
    {
        const Real root_i( get_root( i ) );
        table.push_back( p_exp_den_i(t, root_i, gsl_pow_2( root_i ) ) );
        i++;
    }
}


//Integrated Greens function for rr part of [-Ll, 0]
Real GreensFunction1DAbsSinkAbs::p_int_r_leftdomain(uint i, 
                                                    Real const& rr,
                                                    RealVector const& table) const
{
    const Real root_i( get_root( i ) );
    const Real LrmL0( getLr() - getL0() );
    const Real Llprr( getLl() + rr );
   
    const Real temp( getD() * sin( root_i * LrmL0 ) * 
                     ( cos( root_i * Llprr ) - 1.0 ) );

    return table[i] * temp;
}


//Integrated Greens function for rr part of (0, L0]
Real GreensFunction1DAbsSinkAbs::p_int_r_rightdomainA(uint i, 
                                                      Real const& rr, 
                                                      RealVector const& table) const
{
    const Real root_i( get_root( i ) );
    const Real LrmL0( getLr() - getL0() );
    const Real Llprr( getLl() + rr );
    const Real root_i_rr( root_i * rr );
    
    const Real temp( getD() * ( cos( root_i * Llprr ) - 1.0 ) + 
                     getk() / root_i * ( cos( root_i_rr ) - 1.0 ) 
                     * sin( root_i * getLl() ) );
    
    return table[i] * sin( root_i * LrmL0 ) * temp;
}


//Integrated Greens function for rr part of (L0, Lr]
Real GreensFunction1DAbsSinkAbs::p_int_r_rightdomainB(uint i, 
                                                      Real const& rr, 
                                                      RealVector const& table) const
{
    const Real root_i( get_root( i ) );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real L( Lr + Ll );
    const Real LrmL0( Lr - L0 );
    const Real Lrmrr( Lr - rr );
    const Real LlpL0( Ll + L0 );
                   
    const Real term1( sin( root_i * L ) - sin( root_i * LrmL0 ) - 
            sin( root_i * LlpL0 ) * cos( root_i * Lrmrr ) );
            
    const Real term2( sin( root_i * Lr ) - sin( root_i * LrmL0 ) -
            sin( root_i * L0 ) * cos( root_i * Lrmrr ));
        
    const Real temp( getD() * term1 + 
                     getk() * sin( root_i * Ll ) * term2 / root_i );

    return table[i] * temp;
}


/* Function for GFL rootfinder of drawR. */
Real GreensFunction1DAbsSinkAbs::drawR_f(Real r, void *p)
{
    struct drawR_params *params = (struct drawR_params *)p;
    return params->gf->p_int_r_table(r, params->t, params->table) - params->rnd;
}


Real GreensFunction1DAbsSinkAbs::drawR(Real rnd, Real t) const
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd <= 1.0 );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    
    const Real D( getD() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real r0( getr0() );
    const Real L( Lr + Ll );

    if (t == 0.0 || (D == 0.0 && v == 0.0) )
    {
	    // the trivial case
	    return r0;
    }

    if ( L < 0.0 )
    {
	    // if the domain had zero size
	    return 0.0;
    }

    if ( rnd <= EPSILON )
    {
        return getsigma();        
    }

    if( rnd >= ( 1 - EPSILON ) )
    {
        return geta();
    }

    // the structure to store the numbers to calculate r.
    const Real S( p_survival( t ) );

    RealVector drawR_table;
    drawR_params parameters = { this, t, drawR_table, rnd * S };

    // define gsl function for rootfinder
    gsl_function F;
    F.function = &drawR_f;
    F.params = &parameters;

    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );

    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    Real r( findRoot( F, solver, getsigma(), geta(), 
                      EPSILON*L, EPSILON, "GreensFunction1AbsSinkAbs::drawR" ) );

    // Convert the position rr to 'world' coordinates and return it.
    return r;
}


std::string GreensFunction1DAbsSinkAbs::dump() const
{
    std::ostringstream ss;
    ss << "D = " << getD() << ", sigma = " << getsigma() <<
        ", a = " << geta() <<
        ", r0 = " << getr0() <<
        ", rsink = " << getrsink() <<
        ", k = " << getk() << std::endl;
    return ss.str();
}

