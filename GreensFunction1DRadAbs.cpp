#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>

#include <boost/bind.hpp>
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
#include "GreensFunction1DRadAbs.hpp"

// This is the appropriate definition of the function defining
// the roots of our Green's functions in GSL.
// Later needed by the rootfinder.
//
// It expects a reaction rate h=k/D already divided by D.
double
GreensFunction1DRadAbs::tan_f (double x, void *p)
{
    // casts the void to the struct pointer
    struct tan_f_params *params = (struct tan_f_params *)p;
    const Real a = (params->a);
    const Real h = (params->h);
    const Real h_a (h*a);

    if ( fabs( h_a ) < 1 )
    {
        // h = k/D
        return 1/tan(x) + (h_a)/x;
    }
    else
    {
        // h = k/D
        return tan(x) + x/(h_a);
    }
    
}

/* Fills the rootList with all the roots of tan(x*a)=-x/h up to n */
void GreensFunction1DRadAbs::calculate_n_roots(uint const& n) const
{
    uint i( rootList_size() );
    if( n <= i )
        return;

    const Real L( this->geta()-this->getsigma() );
    const Real h( (this->getk()+this->getv()/2.0) / this->getD() );
    // the drift v also comes into this constant, h=(k+v/2)/D
    Real upper, lower, root_i;

    //No drift, and k = 0, use reflective solution.
    if (getk() < EPSILON && fabs( getv() ) < EPSILON )
    {
        while(i < n)
        {
            ad_to_rootList( M_PI * ( i + 1.0/2 ) / L );
            i++;
        }
        return;
    }
    
    gsl_function F;
    struct tan_f_params params = { L, h };
     
    F.function = &GreensFunction1DRadAbs::tan_f;
    F.params = &params;

    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    /* Find all the roots up to the nth */
    if ( h*L < 1 )
    {
        lower = i*M_PI + 1E-10;
        upper = ( i + 1 ) * M_PI - 1E-10;
    }
    else
    {
        lower = i * M_PI + M_PI_2 + 1E-10;
        upper = ( i + 1 ) * M_PI + M_PI_2 - 1E-10;
    }

    while(i++ < n)
    {
        root_i = findRoot( F, solver, lower, upper, 
                           1.0*EPSILON, EPSILON, "GreensFunction1DRadAbs::root_tan" );

        ad_to_rootList( root_i / L );

        lower += M_PI;
        upper += M_PI;
    }

    gsl_root_fsolver_free( solver );
}


/* returns a guess for the number of terms needed for 
   the greensfunction to converge at time t */
uint GreensFunction1DRadAbs::guess_maxi(Real const& t) const
{
    const uint safety(2);

    if (t >= INFINITY)
    {
        return safety;
    }

    const Real D( getD() );
    const Real L( fabs( geta() - getsigma() ) );

    const Real root0( get_root( 0 ) );
    const Real Dt(D * t);

    const Real thr(exp(- Dt * root0 * root0) * EPSILON * 1e-1);

    if (thr <= 0.0)
    {
        return MAX_TERMS;
    }

    const Real max_root( sqrt(root0 * root0 - log(thr) / Dt) );

    const uint maxi(std::max( safety + 
                              static_cast<uint>
                              (max_root * L  / M_PI),
                              MIN_TERMS )
                    );

    return std::min(maxi, MAX_TERMS);
}


// This is the non-exponential factor in the Green's function sum, not
// including the factor containing the explicit r-dependency (The latter
// is given by the Bn's, see below).
//
// r0 is here still in the interval from 0 to a (and supposed to be the
// starting point of the particle at t0).
//
// The root a_n also must be the specific one for that interval, thus
// the one rescaled by a (see comments in function a_n(n) ).
//
// The factor calculated here is identical for the cases w. or w/o drift,
// only h changes.
Real
GreensFunction1DRadAbs::An (Real root_n) const
{
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    const Real sigma(this->getsigma());
    const Real L(this->geta()-this->getsigma());
    const Real r0(this->getr0());
    const Real rootn_r0_s = root_n*(r0-sigma);

    return (root_n*cos(rootn_r0_s) + h*sin(rootn_r0_s)) / (h + (root_n*root_n + h*h)*L);
}

// This factor appears in the survival prob.
Real
GreensFunction1DRadAbs::Bn (Real root_n) const
{
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    const Real k(this->getk());
    const Real D(this->getD());
    const Real v(this->getv());
    const Real sigma(this->getsigma());
    const Real a(this->geta());
    const Real L(this->geta()-this->getsigma());
    
    const Real rootnL(root_n*L);
    const Real rootn2(root_n*root_n);
    const Real h2(h*h);
    const Real v2D(v/2.0/D);

    if(v==0.0) return (h2 - (rootn2 + h2)*cos(rootnL)) / (h*root_n);
    else	return (exp(v2D*sigma)*h*k/D - exp(v2D*a)*(rootn2+h2)*cos(rootnL) ) / (h/root_n*(rootn2+v2D*v2D));
}

// This is the exponential factor in the Green's function sum, also
// appearing in the survival prob. and prop. function.
//
// Also here the root is the one refering to the interval of length L.
Real GreensFunction1DRadAbs::Cn (Real root_n, Real t)
const
{
    const Real D(this->getD());

    return exp(-D*root_n*root_n*t);
}


Real GreensFunction1DRadAbs::p_survival(Real t) const
{
    RealVector table;
    return p_survival_table(t, table);
}


/* Calculates survival probability using a table. 
   Switchbox for which greensfunction to use. */
Real GreensFunction1DRadAbs::p_survival_table(Real t, RealVector& psurvTable) const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    Real p;

    const Real a( geta() );
    const Real sigma( getsigma() );
    const Real L(a - sigma );
    const Real r0( getr0() );
    const Real D( getD() );
    const Real v( getv() );

    if ( fabs(a-r0) < L*EPSILON || L < 0.0 )
    {
        // The survival probability of a zero domain is zero
        return 0.0;
    }

    if (t == 0.0 || (D == 0.0 && v == 0.0) )
    {
	    //particle can't escape.
	    return 1.0;
    }

    /* First check if we need full solution. 
       Else we use approximation. */
    const Real distToa( a - r0 );
    const Real distTos( r0 - sigma );
    const Real maxDist( CUTOFF_H * ( sqrt(2.0 * D * t) + fabs(v * t) ) );

    if( distToa > maxDist ) //Absorbing boundary 'not in sight'.
    {
        if( distTos > maxDist ) //Radiation boundary 'not in sight'.
            return 1.0; //No prob. outflux.
        else
            return XS30(t, distTos, getk(), D, v); //Only radiation BCn.
    }
    else
    {
        if( distTos > maxDist )
            return XS10(t, distToa, D, -v); //Only absorbing BCn.
    }

    const uint maxi( guess_maxi(t) );
    
    if( maxi >= MAX_TERMS )
        log_.warn("drawT: maxi was cut to MAX_TERMS for t = %.16g", t);
    
    if ( psurvTable.size() < maxi )
    {
        calculate_n_roots( maxi );
        createPsurvTable( psurvTable );
    }

    p = funcSum_all(boost::bind(&GreensFunction1DRadAbs::p_survival_i, 
                                this, _1, t, psurvTable),
                    maxi);
    
    if( v == 0.0 )
    {   
        p *= 2.0;
    }
    else
    {       
        const Real vexpo(-v*v*t/4.0/D - v*r0/2.0/D);
        p *= 2.0 * exp( vexpo );
    }

    return p;
}


/* Calculates the i'th term of the p_survival sum */
Real GreensFunction1DRadAbs::p_survival_i( uint i, 
                                           Real const& t, 
                                           RealVector const& table) const
{
    return exp( - getD() * t * gsl_pow_2( get_root( i ) ) ) * table[ i ];
}


/* Calculates the part of the i'th term of p_surv not dependent on t, with drift */
Real GreensFunction1DRadAbs::p_survival_table_i_v( uint const& i ) const
{
    const Real sigma( getsigma() );
    const Real L( geta() - sigma );
    const Real r0( getr0() );
    const Real D( getD() );
    const Real v( getv() );
    const Real h( (getk() + v / 2.0) / D );

    const Real v2D(v/2.0/D);
    const Real exp_av2D(exp(a*v2D));
    const Real exp_sigmav2D(exp(sigma*v2D));
    
    const Real root_n( get_root( i ) );	
    const Real root_n2     = root_n * root_n;
    const Real root_n_r0_s = root_n * (r0-sigma);
    const Real root_n_L    = root_n * L;
    const Real h_root_n    = h / root_n;
    
    return ( h * sin( root_n_r0_s ) + root_n * cos( root_n_r0_s ) ) / 
        ( L * ( root_n2 + h * h) + h ) * ( exp_sigmav2D * h * k / D - 
                               exp_av2D * ( root_n2 + h * h )*
                               cos( root_n_L ) ) / 
            ( h_root_n * (root_n2 + v2D * v2D) );
}


/* Calculates the part of the i'th term of p_surv not dependent on t, without drift */
Real GreensFunction1DRadAbs::p_survival_table_i_nov( uint const& i ) const
{
    const Real sigma( getsigma() );
    const Real L( geta() - sigma );
    const Real r0( getr0() );
    const Real h( getk()/getD() );

    const Real root_n ( get_root( i ) );	
    const Real root_n2( root_n * root_n );
    const Real root_n_r0_s( root_n * (r0-sigma) );
    const Real root_n_L( root_n * L );
    const Real h_root_n( h / root_n );
    
    return (h*sin(root_n_r0_s) + root_n*cos(root_n_r0_s)) 
        / (L*(root_n2+h*h)+h) * ( h_root_n + sin(root_n_L) 
                                  - h_root_n*cos(root_n_L) );
}


/* Fills table with terms in the p_survival sum which don't depend on t */
void GreensFunction1DRadAbs::createPsurvTable( RealVector& table) const
{
    const uint root_nbr( rootList_size() );
    uint i( table.size() );

    if( getv() == 0.0 )
    {
        while( i < root_nbr )
        {
            table.push_back( p_survival_table_i_nov( i++ ) );
        }
    }
    else
    {
        while( i < root_nbr )
        {
            table.push_back( p_survival_table_i_v( i++ ) );
        }
    }
}


/* Calculates the probability density of finding the particle at location r
   at time t. */
Real GreensFunction1DRadAbs::prob_r (Real r, Real t) const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    THROW_UNLESS( std::invalid_argument, (r-sigma) >= 0.0 && r <= a 
                  && (r0 - sigma) >= 0.0 && r0<=a );
    
    const Real sigma(this->getsigma());
    const Real a(this->geta());
    const Real L(this->geta()-this->getsigma());
    const Real r0(this->getr0());
    const Real D(this->getD());
    const Real v(this->getv());
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    
    const Real vexpo(-v*v*t/D/4.0 + v*(r-r0)/D/2.0);

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

    // if r is at the absorbing boundary
    if ( fabs(a-r) < EPSILON*L )
    {
        return 0.0;
    }

    Real root_n, root_n_r_s;
    Real sum = 0, term = 0, prev_term = 0;

    const uint maxi( guess_maxi( t ) );
    calculate_n_roots( maxi );

    uint n = 0;
    do
    {
        if ( n >= MAX_TERMS )
        {
            log_.warn("Too many terms needed for prob_r. N: %5u", n);
            break;
        }

        root_n = this->get_root(n);
        root_n_r_s = root_n*(r-sigma);

        prev_term = term;
        term = Cn(root_n, t) * An(root_n) * (h*sin(root_n_r_s) + root_n*cos(root_n_r_s));
        sum += term;

        n++;
    }
    while (fabs(term/sum) > EPSILON*PDENS_TYPICAL || 
           fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
           n < MIN_TERMS );

    return 2.0*exp(vexpo)*sum;
}


/* Calculates the probability density of finding the particle at location z at
   timepoint t, given that the particle is still in the domain. */
Real GreensFunction1DRadAbs::calcpcum (Real r, Real t) const
{
    return prob_r(r, t)/p_survival(t);
}


/* Calculates the total probability flux leaving the domain at time t
   This is simply the negative of the time derivative of the survival prob.
   at time t [-dS(t')/dt' for t'=t]. */
Real GreensFunction1DRadAbs::flux_tot (Real t) const
{
    Real root_n;
    const Real D(this->getD());
    const Real v(this->getv());
    const Real vexpo(-v*v*t/4.0/D - v*r0/2.0/D);

    const Real D2 = D*D;
    const Real v2Dv2D = v*v/4.0/D2;
    double sum = 0, term = 0, prev_term = 0;

    const uint maxi( guess_maxi( t ) );
    calculate_n_roots( maxi );

    uint n = 0;
    do
    {
        if ( n >= MAX_TERMS )
        {
            log_.warn("Too many terms needed for flux_tot. N: %5u", n );
            break;
        }

        root_n = this->get_root(n);
        prev_term = term;
        term = (root_n * root_n + v2Dv2D) * Cn(root_n, t) * An(root_n) * Bn(root_n);
        n++;
        sum += term;
    }
    while (fabs(term/sum) > EPSILON*PDENS_TYPICAL ||
           fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
           n < MIN_TERMS );

    return 2.0*D*exp(vexpo)*sum;
}


/* Calculates the probability flux leaving the domain through the radiative
   boundary at time t */
Real GreensFunction1DRadAbs::flux_rad (Real t) const
{
    return this->getk() * prob_r(this->getsigma(), t);
}


/* Calculates the flux leaving the domain through the radiative boundary as a
   fraction of the total flux. This is the probability that the particle left
   the domain through the radiative boundary instead of the absorbing
   boundary. */
Real GreensFunction1DRadAbs::fluxRatioRadTot (Real t) const
{
    return flux_rad(t)/flux_tot(t);
}

/* Determine which event has occured, an escape or a reaction. Based on the
   fluxes through the boundaries at the given time. Beware: if t is not a
   first passage time you still get an answer! */
GreensFunction1DRadAbs::EventKind
GreensFunction1DRadAbs::drawEventType( Real rnd, Real t )
const
{
    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, t > 0.0 );
    // if t=0 nothing has happened => no event

    const Real a(this->geta());
    const Real sigma( this->getsigma() );
    const Real L(this->geta()-sigma);
    const Real r0(this->getr0());

    // if the radiative boundary is impermeable (k==0) or
    // the particle is at the absorbing boundary (at a) => IV_ESCAPE event
    if ( k == 0 || fabs(a-r0) < EPSILON*L )
    {
        return IV_ESCAPE;
    }

    /* First check if we need to compare flux ratio's.
       If only one boundary is 'visible' to the particle, use this boudnary as escape.*/
    const Real distToa( a - r0 );
    const Real distTos( r0 - sigma );
    const Real maxDist( CUTOFF_H * ( sqrt(2.0 * D * t) + fabs(v * t) ) );

    if( distToa > maxDist ) //Absorbing boundary 'not in sight'.
    {
        if( distTos < maxDist ) //Only radiation boundary 'in sight'.
            return IV_REACTION;
    }
    else
    {
        if( distTos > maxDist ) //Only absorbing boundary 'in sight'.
            return IV_ESCAPE;
    }

    // Else the event is sampled from the flux ratio
    const Real fluxratio (this->fluxRatioRadTot(t));

    if (rnd > fluxratio )
    {
        return IV_ESCAPE;
    }
    else
    {
        return IV_REACTION;
    }
}


/* This function is needed to cast the math. form of the function
   into the form needed by the GSL root solver. */
Real GreensFunction1DRadAbs::drawT_f (double t, void *p)
{
    struct drawT_params *params = (struct drawT_params *)p;
    return params->rnd - params->gf->p_survival_table( t, params->psurvTable );
}

/* Draws the first passage time from the survival probability,
   using an assistance function drawT_f that casts the math. function
   into the form needed by the GSL root solver. */
Real GreensFunction1DRadAbs::drawTime (Real rnd) const
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
  
    const Real sigma(this->getsigma());
    const Real a(this->geta());
    const Real L(this->geta()-this->getsigma());
    const Real r0(this->getr0());
    const Real k(this->getk());
    const Real D(this->getD());

    if ( D == 0.0 || L == INFINITY )
    {
        return INFINITY;
    }

    if ( rnd > 1 - EPSILON || L < 0.0 || fabs(a-r0) < EPSILON*L )
    {
        return 0.0;
    }

    /* Find a good interval to determine the first passage time. */
    Real t_guess;

    if( k != 0.0 )
    {
        const Real t_Abs( gsl_pow_2( a - r0 ) / D );
        const Real t_Rad( D / (k * k) + gsl_pow_2( r0 - sigma ) / D );

        t_guess = std::min(t_Abs, t_Rad );
    }
    else
    {
        t_guess = gsl_pow_2( a - r0 ) / D;
    }

    t_guess *= .1;

    // A different guess has to be made in case of nonzero drift to account for the displacement due to it
    // TODO: This does not work properly in this case yet, but we don't know why...
    // When drifting towards the closest boundary
    //if( (r0 >= a/2.0 && v > 0.0) || (r0 <= a/2.0 && v < 0.0) )	t_guess = sqrt(D*D/(v*v*v*v)+dist*dist/(v*v)) - D/(v*v);
    // When drifting away from the closest boundary
    //if( ( r0 < a/2.0 && v > 0.0) || ( r0 > a/2.0 && v < 0.0) )	t_guess = D/(v*v) - sqrt(D*D/(v*v*v*v)-dist*dist/(v*v));  

    /* Set params structure. */
    RealVector exponent_table;
    RealVector psurvTable;
    struct drawT_params parameters = {this, psurvTable, rnd};
    
    /* Define the function for the rootfinder. */
    gsl_function F;
    F.function = &GreensFunction1DRadAbs::drawT_f;
    F.params = &parameters;

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
                log_.warn("drawTime: couldn't adjust low. F( %.16g ) = %.16g"
                          , low, value);
                /*
                std::cerr << "GF1DRad::drawTime Couldn't adjust low. F(" << low << ") = "
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

    /* define a new solver type brent */
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    /* make a new solver instance */
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    const Real t( findRoot( F, solver, low, high, t_scale*EPSILON, EPSILON,
                            "GreensFunction1DRadAbs::drawTime" ) );

    // return the drawn time
    return t;
}

/* Returns c.d.f. for drawR */
Real GreensFunction1DRadAbs::p_int_r_table(Real const& r, 
                                           Real const& t,
                                           RealVector& table) const
{
    const Real a( geta() );
    const Real sigma( getsigma() );
    const Real L( a - sigma );
    const Real r0( getr0() );
    const Real D( getD() );
    const Real k( getk() );
    const Real v( getv() );

    /* If not all boundaries are 'visible' to the particle,
       use approximation. */
    const Real distToa( a - r0 );
    const Real distTos( r0 - sigma );
    const Real maxDist( CUTOFF_H * ( sqrt(2.0 * D * t) + fabs(v * t) ) );

    //TODO: include XI30 (c.d.f) with drift.
    if( distToa > maxDist ) //Absorbing boundary 'not in sight'.
    {
        if( distTos > maxDist ) //Radiation boundary 'not in sight'.
            return XI00(r, t, r0, D, v); //free particle.
        else
        {
            if( k != 0.0 && v == 0.0 )
                //Only radiation BCn.
                return XI30(r - sigma, t, distTos, getk(), D, 0.0);
            else if( k == 0.0 && v == 0.0 )
                //Only reflecting BCn.
                return XI20(r - sigma, t, distTos, D, 0.0);
        }
    }
    else
    {
        if( distTos > maxDist )
            //Only absorbing BCn.
            return XI10(a - r, t, distToa, D, -v);
    }

    Real p;
    const uint maxi( guess_maxi( t ) );

    const Real vexpo(-v*v*t/4.0/D - v*r0/2.0/D);
    const Real prefac( 2.0*exp(vexpo) );

    if( maxi >= MAX_TERMS )
    {
        log_.warn("drawR: maxi was cut to MAX_TERMS for t = %.16g", t);
        std::cerr << dump();
        std::cerr << "L: " <<  L << " r0: " << r0 - sigma << std::endl;
    } 

    if( table.size() < maxi )
    {
        calculate_n_roots( maxi );
        create_p_int_r_Table(t, table);
    }

    p = funcSum(boost::bind(&GreensFunction1DRadAbs::p_int_r_i, 
                            this, _1, r, t, table),
                MAX_TERMS);

    return prefac * p;
}

Real GreensFunction1DRadAbs::p_int_r_i(uint i, 
                                       Real const& r, 
                                       Real const& t, 
                                       RealVector& table) const
{
    const Real sigma( getsigma() );
    const Real D( getD() );
    const Real k( getk() );
    const Real v( getv() );
    const Real h(( k + v/2.0)/ D );
    const Real v2D( v/(2*D) );

    const Real costerm( k / D );	
    const Real sinterm( h * v2D );

    const Real expsigma(exp(sigma*v2D));
    const Real zs(r - sigma);

    Real root_n( get_root( i ) );
    
    Real term( ( expsigma*costerm - exp(v2D*r)*
                 ( costerm*cos(root_n*zs) - 
                   (root_n+sinterm/root_n)*sin(root_n*zs) )) );
    
    return get_p_int_r_Table_i(i, t, table) * term;
}

/* Fills table for p_int_r of factors independent of r. */
void GreensFunction1DRadAbs::create_p_int_r_Table( Real const& t,
                                                   RealVector& table ) const
{
    const uint root_nmbr( rootList_size() );
    uint n( table.size() );

    const Real sigma( getsigma() );
    const Real L( geta() - getsigma() );
    const Real r0( getr0() );
    const Real D( getD() );
    const Real v( getv() );
    const Real h(( k + v/2.0)/ D );

    const Real v2D(v/2.0/D);
    const Real v2Dv2D(v2D*v2D);

    Real term;
    Real root_n2, root_n_r0_s, root_n;

    while( n < root_nmbr )
    {
        root_n = get_root(n);
       	root_n2 = root_n * root_n;
       	root_n_r0_s = root_n * (r0-sigma);

       	term = exp(-D*root_n2*t)
            * (root_n*cos(root_n_r0_s) + h*sin(root_n_r0_s)) / (L*(root_n2 + h*h) + h)
            * root_n / (root_n2 + v2Dv2D);

        table.push_back( term );
        n++;
    }
}

/* Function for GSL rootfinder of drawR. */
Real GreensFunction1DRadAbs::drawR_f(Real r, void *p)
{
    struct drawR_params *params = (struct drawR_params *)p;
    return params->gf->p_int_r_table(r, params->t, params->table) 
        - params->rnd;
}

/* Return new position */
Real GreensFunction1DRadAbs::drawR (Real rnd, Real t) const
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    
    const Real sigma( getsigma() );
    const Real a( geta() );
    const Real L( a - sigma );
    const Real r0( getr0() );
    const Real D( getD() );
    const Real v( getv() );
    
    if (t == 0.0 || (D == 0.0 && v == 0.0) )
    {
        return r0;
    }
    if ( a < 0.0 )
    {
        return 0.0;
    }

    // the structure to store the numbers to calculate the numbers for 1-S
    RealVector pintTable;
    struct drawR_params parameters = {this, t, pintTable, rnd * p_survival( t )};

    // define gsl function for rootfinder
    gsl_function F;
    F.function = &GreensFunction1DRadAbs::drawR_f;
    F.params = &parameters;

    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    Real r( findRoot( F, solver, sigma, a, EPSILON*L, EPSILON,
                            "GreensFunction1DRadAbs::drawR" ) );

    // return the drawn position
    return r;
}

std::string GreensFunction1DRadAbs::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << ", sigma = " << this->getsigma() <<
        ", a = " << this->geta() <<
        ", r0 = " << this->getr0() <<
        ", k = " << this->getk() << std::endl;
    return ss.str();
}

Logger& GreensFunction1DRadAbs::log_(
        Logger::get_logger("GreensFunction1DRadAbs"));
