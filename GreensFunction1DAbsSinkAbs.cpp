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
#include "GreensFunction1DAbsSinkAbs.hpp"

const unsigned int GreensFunction1DAbsSinkAbs::MAX_TERMS;
const unsigned int GreensFunction1DAbsSinkAbs::MIN_TERMS;

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
    uint i( rootList_size() );
    if( n <= i )
        return;

    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L( Lr + Ll );
    const Real Lm( Lr - Ll );
    const Real Lm_L( Lm / L );
    const Real h( getk() * L / ( 2 * getD() ) );
    
    Real root_i;
    real_pair lower_upper_pair;

    /* Define the root function. */
    gsl_function F;
    struct root_f_params params = { Lm_L, h };
     
    F.function = &GreensFunction1DAbsSinkAbs::root_f;
    F.params = &params;

    /* define and ad a new solver type brent */
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    /* If this is the first run, set parameters.*/
    if(i == 0)
    {
        lo_up_params.h = h;
        lo_up_params.Lm_L = Lm_L;
        lo_up_params.long_period = std::max( L/Lr * M_PI, L/Ll * M_PI );
        lo_up_params.short_period = std::min( L/Lr * M_PI, L/Ll * M_PI );    
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
        log_.warn("f(lower) has wrong sign at root# %6u, for h = %.5g, Lm/L = %.5g.",
                  rootList_size() + 1, lo_up_params.h, lo_up_params.Lm_L );

        log_.warn("f_low( %.16g ) =  %.16g , f_high( %.16g ) = %.16g",
                  lower, f_lower, upper, f_upper);
    }

    /* If the parity is incorrect -> correct it */
    if( f_upper * parity_op < 0)
    {
        int cntr( 0 );

        Real delta( .1 * safety * right_offset );
        const Real save_upper( upper );
        
        /* Assuming the upper point has overshoot the straddle region,
           subtract from the upper limit, until endpoits do straddle.
         */
        while(f_upper * parity_op < 0 && cntr++ < 10 )
        {
            //Correct for overshoot.
            upper -= delta;
            f_upper = root_f( upper, &p );
        }

        /* If upper point still doesn't straddle the root, increase the old estimate 
           of upper, until it does straddle*/
        if(cntr >= 10)
        {
            cntr = 0;
            upper = save_upper;
            f_upper = root_f( upper, &p );

            while(f_upper * parity_op < 0 && cntr++ < 10 )
            {
                //Correct for undershoot.
                upper += delta;
                f_upper = root_f( upper, &p );
            }

            // Still no straddle? => Crash is imminent.
            if(cntr >= 10)
            {
                log_.warn("Failed to straddle root# %6u. ",
                          rootList_size() + 1);
                log_.warn("next root est. = %.16g with stepsize: %.16g, ll: %.16g, ls: %.16g.", 
                          next_root_est, delta, lo_up_params.last_long_root, lo_up_params.last_short_root);
                log_.warn("f_low( %.16g ) =  %.16g, f_high( %.16g ) = %.16g.",
                          lower, f_lower, upper, f_upper);
            }
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
    const Real sigma(getsigma());
    const Real a( geta() );
    const Real L0( getL0() ); 

    if (t == 0.0 || D == 0.0 )
    {
	    //particle can't escape.
	    return 1.0;
    }

    /* First check if we need full solution. 
       Else we use approximation. */
    const Real maxDist(CUTOFF_H * sqrt(2.0 * D * t));
    // dist to nearest absorbing boundary.
    const Real distToAbs( std::min(a - r0, r0 -sigma) );

    if( L0 > maxDist ) //Sink not in sight
    {
        if( distToAbs > maxDist )
            return 1.0;
        else
            return XS10( t, distToAbs, D, 0.0 );  
    }
    else
    {
        if( distToAbs > maxDist ) //Only sink in sight.
        {
            return XS030( t, L0, getk(), D );
        }
    }

    const uint maxi( guess_maxi(t) );
    
    if( maxi == MAX_TERMS )
        log_.info("drawT: maxi was cut to MAX_TERMS for t = %.16g", t);
        
    if (psurvTable.size() < maxi )
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

    while( i < root_nbr )
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

    /* Determine wether rr lies in the same sub-domain as r0.
       A different function is caculated when this is the case. */
    if( rr >= 0 )
    {   
        p = funcSum(boost::bind(&GreensFunction1DAbsSinkAbs::prob_r_r0_i, 
                                this, _1, rr, t),
                    MAX_TERMS);
	}
    else
    {
        p = funcSum(boost::bind(&GreensFunction1DAbsSinkAbs::prob_r_nor0_i, 
                                this, _1, rr, t),
                    MAX_TERMS);
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

    if( t == 0 || ( D == 0 && getr0() != getrsink() ) )
    {
        return 0.0;
    }

    Real p;

    p = funcSum(boost::bind(&GreensFunction1DAbsSinkAbs::flux_tot_i, 
                            this, _1, t),
                MAX_TERMS);

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
   
    p = funcSum(boost::bind(&GreensFunction1DAbsSinkAbs::flux_abs_Lr_i, 
                            this, _1, t),
                MAX_TERMS);

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

    p = funcSum(boost::bind(&GreensFunction1DAbsSinkAbs::flux_abs_Ll_i, 
                            this, _1, t),
                MAX_TERMS);
    
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
    const Real D( getD() );

    if( t == 0 || ( D == 0 && getr0() != getrsink() ) )
    {
        return 0.0;
    }

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

    /* If the sink is impermeable (k==0) or
       the particle is at one the absorbing boundaries (sigma or a) => IV_ESCAPE event */
    if ( k == 0 || fabs(a - r0) < EPSILON * L || fabs(sigma - r0) < EPSILON * L)
    {
	    return IV_ESCAPE;
    }

    /* The event is sampled from the flux ratios.
       Two possiblities:
       (1) Leave through left or right boundary - IV_ESCAPE
       (2) Leave through sink - IV_REACTION
    */

    /* First check if we need to compare flux ratio's. 
       If particle is near only one boundary or sink, this is the escape event. */
    const Real maxDist(CUTOFF_H * sqrt(2.0 * getD() * t));
    // dist to nearest absorbing boundary.
    const Real distToAbs( std::min(a - r0, r0 -sigma) );

    if( getL0() > maxDist ) //Sink not in sight
    {
        if( distToAbs < maxDist )
            return IV_ESCAPE;
    }
    else
    {
        if( distToAbs > maxDist ) //Only sink in sight.
            return IV_REACTION;
    }

    /* Allready fill rootList with needed roots. */ 
    const uint maxi( guess_maxi( t ) );

    if( maxi == MAX_TERMS )
        log_.info("drawEventType: maxi was cut to MAX_TERMS for t = %.16g", t);

    calculate_n_roots( maxi );

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
    //const Real k( getk() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real L( getLr() + getLl() );

    if ( D == 0.0 || L == INFINITY )
    {
	    return INFINITY;
    }

    if ( rnd > (1 - EPSILON) || L < 0.0 || fabs( a - r0 ) < EPSILON * L )
    {
	    return 0.0;
    }

    /* the structure to store the numbers to calculate the numbers for 1-S */
    RealVector psurvTable;
    drawT_params params = { this, psurvTable, rnd };

    /* Find a good interval to determine the first passage time.
       First we get the distance to one of the absorbing boundaries or the sink. */
    Real dist( std::min(Lr - L0, Ll + L0) );
    Real t_guess;

    /* For a particle at contact, the time for the particle to be absorbed by the 
       radiation boundary is 4 D / (k*k). If the particle is at a distace x0 from the
       radiation boundary, we approximate the gues for the next event-time by:
       t_guess = D / (k*k) + x0 * x0 / D.
    */
    const Real t_Abs( gsl_pow_2( dist ) / D );
    const Real t_Rad( 4 * D / (k * k) + gsl_pow_2( L0 ) / D );

    t_guess = std::min(t_Abs, t_Rad );
    t_guess *= .1;

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
            if( fabs( high ) >= t_guess * 1e10 )
            {
                log_.error("drawTime: couldn't adjust high. F( %.16g ) = %.16g"
                          , high, value);
                throw std::exception();
            }
            // keep increasing the upper boundary until the
            // function straddles
            high *= 10;
            value = GSL_FN_EVAL( &F, high );
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
	        if( fabs( low ) <= t_guess * 1e-10 ||
	            fabs(value-value_prev) < EPSILON*1.0 )
	        {
                log_.warn("drawTime Couldn't adjust low. F( %.16g ) = %.16g"
                          , low, value);
                /*
                  std::cerr << "GF1DSink::drawTime Couldn't adjust low. F(" << low << ") = "
                  << value << " t_guess: " << t_guess << " diff: "
                  << (value - value_prev) << " value: " << value
                  << " value_prev: " << value_prev << " rnd: "
                  << rnd << std::endl;
                */
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

    /* First check if we need full solution. 
       Else we use approximation. */
    const Real maxDist(CUTOFF_H * sqrt(2.0 * D * t));
    const Real distToa( a - r0 );
    const Real distTos( r0 - sigma );

    if( getL0() > maxDist ) //Sink not in sight
    {
        if( distTos > maxDist )
        {
            if( distToa > maxDist )
                return XI00(r, t, r0, D, 0.0);  //Nothing in sight.
            else
                return XI10(a - r, t, distToa, D, 0.0); //Only a boundary in sight.
        }
        else
        {
            if( distToa > maxDist )
                return XI10( r - sigma, t, distTos, D, 0.0 ); //Only s boundary in sight.
        }
    }
    else
    {
        if( distToa > maxDist && distTos > maxDist ) //Only sink in sight.
        {
            return XI030( rr, t, getL0(), getk(), D );
        }
    }

    const uint maxi( guess_maxi(t) );

    if( maxi == MAX_TERMS )
        log_.warn("p_int_r_table: maxi was cut to MAX_TERMS for t = %.16g"
                  , t);
   
    if (table.size() < maxi )
    {
        calculate_n_roots( maxi );  // this updates the table
        create_p_int_r_Table( t, table );
    }

    /* Determine in which part of the domain rr lies, and
       thus which function to use. */
    Real (GreensFunction1DAbsSinkAbs::*p_int_r_i)
        (uint, Real const&, Real const&, RealVector& table) const = NULL;

    if( rr <= 0 )
        p_int_r_i = &GreensFunction1DAbsSinkAbs::p_int_r_leftdomain;
    else if( rr < getL0() )
            p_int_r_i = &GreensFunction1DAbsSinkAbs::p_int_r_rightdomainA;
        else
            p_int_r_i = &GreensFunction1DAbsSinkAbs::p_int_r_rightdomainB;

    p = funcSum(boost::bind(p_int_r_i, 
                            this, _1, rr, t, table),
                MAX_TERMS);

    return 2.0 * p;
}


Real GreensFunction1DAbsSinkAbs::p_int_r(Real const& r, 
                                         Real const& t) const
{
    THROW_UNLESS( std::invalid_argument, r >= getsigma() && r <= geta() );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    RealVector table;
    return p_int_r_table(r, t, table);
}


void GreensFunction1DAbsSinkAbs::create_p_int_r_Table( Real const& t, 
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
                                                    Real const& t,
                                                    RealVector& table) const
{
    const Real root_i( get_root( i ) );
    const Real LrmL0( getLr() - getL0() );
    const Real Llprr( getLl() + rr );
   
    const Real temp( getD() * sin( root_i * LrmL0 ) * 
                     ( cos( root_i * Llprr ) - 1.0 ) );

    return get_p_int_r_Table_i(i, t, table) * temp;
}


//Integrated Greens function for rr part of (0, L0]
Real GreensFunction1DAbsSinkAbs::p_int_r_rightdomainA(uint i, 
                                                      Real const& rr, 
                                                      Real const& t,
                                                      RealVector& table) const
{
    const Real root_i( get_root( i ) );
    const Real LrmL0( getLr() - getL0() );
    const Real Llprr( getLl() + rr );
    const Real root_i_rr( root_i * rr );
    
    const Real temp( getD() * ( cos( root_i * Llprr ) - 1.0 ) + 
                     getk() / root_i * ( cos( root_i_rr ) - 1.0 ) 
                     * sin( root_i * getLl() ) );
    
    return get_p_int_r_Table_i(i, t, table) * sin( root_i * LrmL0 ) * temp;
}


//Integrated Greens function for rr part of (L0, Lr]
Real GreensFunction1DAbsSinkAbs::p_int_r_rightdomainB(uint i, 
                                                      Real const& rr, 
                                                      Real const& t,
                                                      RealVector& table) const
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

    return get_p_int_r_Table_i(i, t, table) * temp;
}


/* Function for GFL rootfinder of drawR. */
Real GreensFunction1DAbsSinkAbs::drawR_f(Real r, void *p)
{
    struct drawR_params *params = (struct drawR_params *)p;
    return params->gf->p_int_r_table(r, params->t, params->table) 
        - params->rnd;
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

    if (t == 0.0 || D == 0.0 )
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
    RealVector drawR_table;
    drawR_params parameters = { this, t, drawR_table, rnd * p_survival( t ) };

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

Logger& GreensFunction1DAbsSinkAbs::log_(
        Logger::get_logger("GreensFunction1DAbsSinkAbs"));
