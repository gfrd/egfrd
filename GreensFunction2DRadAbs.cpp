// Greens function class for 2d Green's Function for 2d annulus with radial and
// axial dependence. Inner boundary is radiative (rad) (reaction event), outer 
// boundary is absorbing (abs) (escape event). Different "draw" functions 
// provide a way to draw certain values from the Green's Function, e.g. an
// escape angle theta ("drawTheta" function).
// 
// Based upon code from Riken Institute. 
// Written by Laurens Bossen, Adapted by Martijn Wehrens. FOM Institute AMOLF.

//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include "compat.h"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <sstream>

#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_integration.h>

#include "factorial.hpp"
#include "funcSum.hpp"
#include "findRoot.hpp"
#include "freeFunctions.hpp"
#include "CylindricalBesselGenerator.hpp"
#include "GreensFunction2DRadAbs.hpp"

const Real GreensFunction2DRadAbs::MIN_T_FACTOR;
const unsigned int GreensFunction2DRadAbs::MAX_ORDER;
const unsigned int GreensFunction2DRadAbs::MAX_ALPHA_SEQ;


// This is the constructor
GreensFunction2DRadAbs::
GreensFunction2DRadAbs( const Real D, 
				    const Real kf,
				    const Real r0, 
				    const Real Sigma,
                    const Real a )
    :
    PairGreensFunction(D, kf, r0, Sigma),   
    h( kf / (D * 2.0 * M_PI * Sigma) ),       
    a( a ),                                    
    estimated_alpha_root_distance_(M_PI/(a-Sigma)) // observed convergence of
                                                   // distance roots f_alpha().                                                                                     
    // ^: Here parent "constructors" are specified that are executed, 
    // also a constructor initialization list can be (and is) specified. 
    // These variables will be set before the contents of the constructor are 
    // called.
{
    const Real sigma(this->getSigma());

    // Check wether input makes sense, outer boundary a should be > inner 
    // boundary sigma.    
    if (a < sigma)
    {
        throw std::invalid_argument((boost::format("a >= sigma : a=%.16g, sigma=%.16g") % a % sigma).str());
    }
    
    // Clear AlphaTables
    GreensFunction2DRadAbs::clearAlphaTable();
    

}

GreensFunction2DRadAbs::~GreensFunction2DRadAbs()
{
    ; // do nothing
}

//
// Alpha-related methods
//

// Resets the alpha-tables
void GreensFunction2DRadAbs::clearAlphaTable() const
{
    const Real estimated_alpha_root_distance(getestimated_alpha_root_distance_());

    // Clears all vectors in the alphaTable
    std::for_each( this->alphaTable.begin(), this->alphaTable.end(),
		   boost::mem_fn( &RealVector::clear ) );
		   
	// Sets all values of the alpha_x_scan_table_ to zero.
    std::fill( this->alpha_x_scan_table_.begin(), 
               this->alpha_x_scan_table_.end(),
	           SCAN_START*estimated_alpha_root_distance); // DEBUG TODO; actual number here might
                     	            // be important for functioning of Bessel
                     	            // functions.
	           
	// Sets all values of the alpha_correctly_estimated_ table to zero.
    std::fill( this->alpha_correctly_estimated_.begin(), 
               this->alpha_correctly_estimated_.end(),
	           0 );	           
	           
}

// The method evaluates the equation for finding the alphas for given alpha. This
// is needed to find the alpha's at which the expression is zero -> alpha is the root.
const Real 
GreensFunction2DRadAbs::f_alpha0( const Real alpha ) const
{
	const Real a( this->geta() );
	const Real sigma( getSigma() );
	const Real h( this->geth() );
	const Real s_An( sigma * alpha );
	const Real a_An( a * alpha );
//	const Real h_s( h * sigma);

	const double J0_s_An (gsl_sf_bessel_J0(s_An));
	const double J1_s_An (gsl_sf_bessel_J1(s_An));
	const double J0_a_An (gsl_sf_bessel_J0(a_An));

	const double Y0_s_An (gsl_sf_bessel_Y0(s_An));
	const double Y1_s_An (gsl_sf_bessel_Y1(s_An));
	const double Y0_a_An (gsl_sf_bessel_Y0(a_An));	

//  const double rho1 ( ( (h_s * J0_s_An) + (s_An * J1_s_An) ) * Y0_a_An );
//	const double rho2 ( ( (h_s * Y0_s_An) + (s_An * Y1_s_An) ) * J0_a_An );
//	return (rho1 - rho2); 
//
//  Sigma can be divided out, roots will remain same:

	const double rho1 ( ( (h * J0_s_An) + (alpha * J1_s_An) ) * Y0_a_An );
	const double rho2 ( ( (h * Y0_s_An) + (alpha * Y1_s_An) ) * J0_a_An );
	return rho1 - rho2;
}


// Wraps function of which we need to find roots alpha.
//
// Params contains pointer to gf object (params.gf), where f_alpha is the member
// function of which we have to find the roots. This function is thus a mere
// wrapper of f_alpha.  
const Real 
GreensFunction2DRadAbs::f_alpha0_aux_F( const Real alpha,
						const f_alpha0_aux_params* const params )
{
    // create pointer to Green's Function object
    const GreensFunction2DRadAbs* const gf( params->gf ); 	
    
    return gf->f_alpha0( alpha );
}


// f_alpha() Calculates the value of the mathematical function f_alpha(). The
// roots (y=0) of this function are constants in the Green's Functions.
const Real GreensFunction2DRadAbs::f_alpha( const Real alpha,
						    const Integer n ) const
{
	const Real a( this->geta() );
	const Real sigma( getSigma() );
	const Real h( this->geth() );
	const Real s_An( sigma * alpha );
	const Real a_An( a * alpha );	
    const Real realn( static_cast<Real>( n ) );

	const double Jn_s_An (gsl_sf_bessel_Jn(n,s_An));
	const double Jn1_s_An (gsl_sf_bessel_Jn(n+1,s_An));
	const double Jn_a_An (gsl_sf_bessel_Jn(n,a_An));

	const double Yn_s_An (gsl_sf_bessel_Yn(n,s_An));
	const double Yn1_s_An (gsl_sf_bessel_Yn(n+1,s_An));
	const double Yn_a_An (gsl_sf_bessel_Yn(n,a_An));

    const double rho1 ( ( (h*sigma * Jn_s_An) + (s_An * Jn1_s_An) - realn*Jn_s_An ) * Yn_a_An );
	const double rho2 ( ( (h*sigma * Yn_s_An) + (s_An * Yn1_s_An) - realn*Yn_s_An ) * Jn_a_An );
	return (rho1 - rho2); 

//  Or..? -:	
//  const double rho1 ( ( ((h*sigma-realn) * Jn_s_An) + (s_An * Jn1_s_An) ) * Yn_a_An );
//	const double rho2 ( ( ((h*sigma-realn) * Yn_s_An) + (s_An * Yn1_s_An) ) * Jn_a_An );
//	return (rho1 - rho2); 		
       
} 

// Simply a wrapper for f_alpha().
const Real 
GreensFunction2DRadAbs::f_alpha_aux_F( const Real alpha,
	       					const f_alpha_aux_params* const params )     					
{
    // Params contains pointer to gf object (params.gf), which has f_alpha() as 
    // a member.
    const GreensFunction2DRadAbs* const gf( params->gf ); 
    const Integer n( params->n );

    return gf->f_alpha( alpha, n );	
}


// calculates the constant part of the i-th term for the survival probability
const Real 
GreensFunction2DRadAbs::p_survival_i( const Real alpha) const
{
	const Real a( geta() );		// get the needed parameters
	const Real sigma( getSigma() );

	const Real s_An (sigma*alpha);
	const Real a_An (a*alpha);
    const Real alpha_sq (alpha*alpha);

	// calculate all the required Bessel functions
    const Real J1_sAn  (gsl_sf_bessel_J1(s_An));
    const Real J0_aAn  (gsl_sf_bessel_J0(a_An));
    const Real Y0_aAn  (gsl_sf_bessel_Y0(a_An));
	const Real Y1_sAn  (gsl_sf_bessel_Y1(s_An));

	// calculate C0,n
        const Real C_i_0 (calc_A_i_0(alpha)); 

	// calculate the integral over Bn,0
	const Real dB_n_0dr (J1_sAn*Y0_aAn - Y1_sAn*J0_aAn);	// this is only the 
	                                   // part without alpha of dB0,n(sigma)/dr
	const Real B_n_0_int (Real(2.0)/(M_PI*alpha_sq) - (sigma/alpha)*dB_n_0dr);

	// return the total result
	const Real result (C_i_0 * B_n_0_int);
	return result;
}

// Calculates the factor An,0 for (for example) determination of the flux 
// through the outer interface
const Real 
GreensFunction2DRadAbs::calc_A_i_0( const Real alpha) const
{
    // Get the required parameters
    const Real a( geta() );    		        
    const Real sigma( getSigma() );
    const Real h (this->geth());

    const Real s_An (sigma*alpha);
    const Real a_An (a*alpha);
    const Real r0An (r0*alpha);

    // Calculate all the required Bessel functions
	const Real J0_sAn  (gsl_sf_bessel_J0(s_An));    
    const Real J1_sAn  (gsl_sf_bessel_J1(s_An));
    const Real J0_aAn  (gsl_sf_bessel_J0(a_An));

    const Real J0_r0An (gsl_sf_bessel_J0(r0An));
    const Real Y0_aAn  (gsl_sf_bessel_Y0(a_An));
    const Real Y0_r0An (gsl_sf_bessel_Y0(r0An));
    
    // Calculate return value
    const Real alpha_sq (alpha*alpha);

    const Real rho (h*J0_sAn + alpha*J1_sAn);
    const Real rho_sq (rho*rho);

    const Real B_n_0 (J0_r0An*Y0_aAn - Y0_r0An*J0_aAn);

    const Real A_i_0 ((alpha_sq * rho_sq * B_n_0)/( rho_sq - (J0_aAn*J0_aAn)*(h*h + alpha_sq)));

	return A_i_0;
}

// Calculates the n-th term of the summation for calculating the flux through 
// the inner interface (reaction)
const Real 
GreensFunction2DRadAbs::leaves_i( const Real alpha) const
{
    const Real a( geta() );         // get the needed parameters
    const Real sigma( getSigma() );

    const Real s_An (sigma*alpha);
    const Real a_An (a*alpha);

    // calculate all the required Bessel functions
    const Real J1_sAn  (gsl_sf_bessel_J1(s_An));
    const Real Y1_sAn  (gsl_sf_bessel_Y1(s_An));
    const Real J0_aAn  (gsl_sf_bessel_J0(a_An));
    const Real Y0_aAn  (gsl_sf_bessel_Y0(a_An));
        

    // calculate An,0
    const Real A_i_0 (calc_A_i_0(alpha));		// calculate the coefficient A0,n

	// calculate dBn,0(sigma)/dr
	const Real dB_n_0dr (-alpha*(J1_sAn*Y0_aAn - Y1_sAn*J0_aAn));

	// calculate the total result
	const Real result (A_i_0 * dB_n_0dr);
        return result;
}


// calculates a table with all the constant factors for the survival probability
void 
GreensFunction2DRadAbs::createPsurvTable( RealVector& table) const
{
    const RealVector& alphaTable_0( this->getAlphaTable( 0 ) );	// get the roots for the survival probability

    table.clear();				// empty the table
    table.reserve( alphaTable_0.size() );	// and get the nescessary memory

    std::transform( alphaTable_0.begin(), alphaTable_0.end(),
		    std::back_inserter( table ),
		    boost::bind( &GreensFunction2DRadAbs::p_survival_i,
				 this, _1) );	// This gets all the roots from 'begin' to 'end'
						// passes them as an argument to p_survival_i and
						// the result is passed to back_inserter
}

// Creates the tables with various Bessel functions used in drawR, the table is used to speed things up
void
GreensFunction2DRadAbs::createY0J0Tables( RealVector& Y0_Table,
							RealVector& J0_Table,
							RealVector& Y0J1J0Y1_Table,
							const Real t ) const
{
	const RealVector& alphaTable_0( this->getAlphaTable( 0 ) );
							// get the roots for the survival probability
	Y0_Table.clear();				// empty the table
	J0_Table.clear();
	Y0J1J0Y1_Table.clear();

	Y0_Table.reserve( alphaTable_0.size() );	// and get the nescessary memory
	J0_Table.reserve( alphaTable_0.size() );
	Y0J1J0Y1_Table.reserve( alphaTable_0.size() );

	boost::tuple<Real,Real,Real> result;
	
	for (unsigned int count = 0; count < alphaTable_0.size(); count++)
	{	result = Y0J0J1_constants(alphaTable_0[count], t);

		Y0_Table.push_back (result.get<0>());
		J0_Table.push_back (result.get<1>());
		Y0J1J0Y1_Table.push_back (result.get<2>());
	}

}


// Creates the values for in the tables Y0, J0 and Y0J1J0Y1
const boost::tuple<Real,Real,Real>
GreensFunction2DRadAbs::Y0J0J1_constants ( const Real alpha,
							const Real t) const
{	const Real D(this->getD());
	const Real sigma(this->getSigma());
	const Real a(this->geta());

	const Real s_An (sigma*alpha);
	const Real a_An (a*alpha);
	const Real alpha_sq (alpha*alpha);

	// calculate all the required Bessel functions	
	const Real J0_aAn  (gsl_sf_bessel_J0(a_An));
	const Real Y0_aAn  (gsl_sf_bessel_Y0(a_An));
	const Real J1_sAn  (gsl_sf_bessel_J1(s_An));
	const Real Y1_sAn  (gsl_sf_bessel_Y1(s_An));
	

    // calculate An,0
	const Real A_i_0 (calc_A_i_0(alpha)); //_sq * rho_sq * B_n_0)/( rho_sq - J0_bAn*J0_bAn*(h*h + alpha_sq)));

	// calculate the exponent with the time
	const Real expT( std::exp(-D*alpha_sq*t));
	// and the product
	const Real Ai0_expT (A_i_0 * expT / alpha);

	// calculate the large constant term in the intergral of Bn,0
	const Real Y0J1_J0Y1 (Y0_aAn*sigma*J1_sAn - J0_aAn*sigma*Y1_sAn);

	return boost::make_tuple (Ai0_expT*Y0_aAn, Ai0_expT*J0_aAn, Ai0_expT*Y0J1_J0Y1);
}


// =============================================================================
// Root finding algorithm for alpha function. 
// Roots (y=0) are necessary to calculate 
// Modified by wehrens@amolf.nl. nov, 2011
// =============================================================================

// Scans for next interval where a sign change is observed, and thus where a 
// root is expected. 
const void GreensFunction2DRadAbs::GiveRootInterval(
            Real& low,             // Variable to return left boundary interval
            Real& high,            // Variable to return right boundary interval
            const Integer n) const // Order of Bessel functions
{
    // Variables for function values @ resp. left and right boundary interval.
    Real f_low , f_high;
    
    // # Get/calculate boundaries and determine width of interval to check for
    // sign change.
    const Real estimated_alpha_root_distance_( this->getestimated_alpha_root_distance_() );
    const Real interval( FRACTION_SCAN_INTERVAL * estimated_alpha_root_distance_);
                                                        
    // If order (n) is zero, the offset is zero, otherwhise take n-1 offset
    // value as first estimate.
    //      (This can be done because the offsets only get bigger, the roots 
    // shift to the right with the order of the Bessel functions n.)   
    if (alpha_x_scan_table_[n] == 0) { // which implies i == 0   
	    if (n > 0)
	    {	
	        alpha_x_scan_table_[n] = ( this->alphaTable[n-1][0] ); 
	    }
	}	    						

    // Define new interval as between x1=offset ("low") and x2=(offset+interval)
    // ("high"). Make sure "low" is > 0.
    low = alpha_x_scan_table_[n];		
    high = alpha_x_scan_table_[n] + interval;	
	if (low <= 0)			        // TODO this check should be redundant
	{	std::cerr << "Left alpha search interval boundary < 0.";//low = EPSILON/L_TYPICAL;
    	throw std::exception();
	}

    // # Look for the sign change:
    
    // Get the values of the function at x "low" and x "high".
    // Different for n=0 because this results in a simpler function.
    //      (Note: This code could be optimized by duplicating all involved  
    // functions and removing this if-statement. 

    if (n == 0) {
        f_low  = f_alpha0(low);	
        f_high = f_alpha0(high);
    } else {
        f_low  = f_alpha(low,n);	
        f_high = f_alpha(high,n);
    }	

    // Continue shifting the search interval until a sign change is detected.
    while( f_low * f_high > 0 )		
    {	      				  
        low =  high;			
        f_low = f_high;			

        high += interval;
        f_high = f_alpha( high, n );
    }
    
    // When above loop has finished, low and high have values inbetween which
    // a root of the alpha function should lie have been found. Make sure that 
    // scanning for the next root starts at the end of the domain [low, high] 
    // found here.
    alpha_x_scan_table_[n] = high;
    
    return;
}

// Simply returns an interval based upon previous root, estimated interval 
// inbetween roots and INTERVAL_MARGIN (see .hpp).
const void GreensFunction2DRadAbs::GiveRootIntervalSimple(
            Real& low,          // Variable to return left boundary interval
            Real& high,         // Variable to return right boundary interval
            const Integer n,    // Order of Bessel functions
            const Real i) const // ith root 
{
    // Offset is simply based on previous root, the interval in which the first 
    // root (i=0) lies is never calculated with this function.
    const Real previous_root (getAlpha(n, i-1));

    // get estimated interval
    const Real estimated_alpha_root_distance_( this->getestimated_alpha_root_distance_() );
 
    // Calculates interval [low, high] where root is expected based on the 
    // assumption where in a converging regime, where the deviation from this
    // estimate is not more than INTERVAL_MARGIN.     
    low  = previous_root + 
                estimated_alpha_root_distance_ * (1 - INTERVAL_MARGIN);
    high = previous_root + 
                estimated_alpha_root_distance_ * (1 + INTERVAL_MARGIN);
    
    return;
}

// This function calls the GSL root finder, for roots for which n = 0. (This is
// a special case for which the function simplifies.)
const Real 
GreensFunction2DRadAbs::getAlphaRoot0( const Real low,  // root lies between low
                                      const Real high   // .. and high 
                                    ) const
{                                    
    // Reinterpret_cast converts any pointer type to any other pointer type, 
    // even of unrelated classes. Hence below code converts pointer to 
    // function ::f_alpha_aux_F() to function pointer, which then points to 
    // memory location of params (&params).

    // f_alpha0_aux_params is a struct: {gf, value}
    f_alpha0_aux_params params = { this, 0 }; // TODO: purpose of this zero is unclear!!!	

    gsl_function F = 
    {       	    
        reinterpret_cast<typeof(F.function)>
        ( &GreensFunction2DRadAbs::f_alpha0_aux_F ),
        &params 
    };				

    // Define solvertype.
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );  
    // Initialize the solver.
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );   

    // Call rootfinder, find root called "alpha" and return it.
    const Real alpha ( findRoot( F, solver, low, high,
		 EPSILON/L_TYPICAL, EPSILON, "GreensFunction2DRadAbs::getAlphaRoot0" ) );
    gsl_root_fsolver_free( solver );
    
    return alpha;
}    

// This function calls the GSL root finder, for roots for which n > 0. (n = 0 is
// a special case for which the function simplifies.)
const Real 
GreensFunction2DRadAbs::getAlphaRootN( const Real low,  // root lies between low
                                      const Real high,  // .. and high 
                                      const Integer n   // nth order Bessel
                                    ) const
{                                    

    // f_alpha_aux_params is a struct: {gf, n, value}
    // n is the summation index (the order of the Bessel functions used
    f_alpha_aux_params params = { this, n, 0 };	// TODO: purpose of this zero is unclear!!!
    
    gsl_function F = 
    {	    
        reinterpret_cast<typeof(F.function)>
        ( &GreensFunction2DRadAbs::f_alpha_aux_F ),
        &params
    };

    // Define solvertype.
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );  
    // Initialize the solver.
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );   

    // Call rootfinder, find root called "alpha" and return it.
    const Real alpha ( findRoot( F, solver, low, high,
		 EPSILON/L_TYPICAL, EPSILON, "GreensFunction2DRadAbs::getAlphaRootN" ) );
    gsl_root_fsolver_free( solver );
    
    return alpha;
}

// Simply calls the correct function to call the rootfinder (either 
// getAlphaRoot0 or getAlphaRootN).
const Real 
GreensFunction2DRadAbs::getAlphaRoot( const Real low,   // root lies between low
                                      const Real high,  // .. and high 
                                      const Integer n   // nth order Bessel
                                    ) const
{
    Real alpha;

    if (n == 0) {
        alpha = getAlphaRoot0(low, high);
    } else {
        alpha = getAlphaRootN(low, high, n);
    }

    return alpha;
}

// Subfunction of getAlpha
// 
// The function just counts the number of times the root lies in the interval 
// that can be used to estimate the next root. If this happens enough, then it
// edits alpha_correctly_estimated_ to a negative value to signify that it can 
// be assumed that for all next roots the distance inbetween roots will
// equal the interval +/- margin.
// TODO: This could be more sophisticated.
const void GreensFunction2DRadAbs::decideOnMethod2(size_t n, 
                                                   RealVector::size_type i
                                                  ) const
{    

    // Since the function can only decide with two alpha's already calculated, 
    // it can't do anything if i = 0.
    if (i > 0) { 

        const Real dx(getAlpha(n, i)-getAlpha(n, i-1)); // note the recursiveness!
        const Real estimated_alpha_root_distance_( this->getestimated_alpha_root_distance_() );

        // If the relative deviation from the expected difference is smaller 
        // than the expected margin, increase number of would-be correct 
        // guesses.
        if (fabs(1-dx/estimated_alpha_root_distance_) < INTERVAL_MARGIN) {
            ++alpha_correctly_estimated_[n];            
        } else {
            alpha_correctly_estimated_[n] = 0;
        }   
        
        // If guessing would have worked for last CONVERGENCE_ASSUMED roots, 
        // assume it will for all following roots.
        if (alpha_correctly_estimated_[n] > CONVERGENCE_ASSUMED) {
            alpha_x_scan_table_[n] = -2; // permanently switch
        }
    }
    
    return;
}

// This function searches for roots (y=0) of the so-called alpha-function 
// (::f_alpha()). It either moves a small search interval along the x-axis to
// check for sign-change (which would indicate a root), and calls the GSL
// root finder, or directly calls the root finder if the spacing between 
// roots is found to be converged. 
const Real GreensFunction2DRadAbs::getAlpha( size_t n,               // order
                                             RealVector::size_type i // ith root
                                           ) const 
{
    const Real sigma(getSigma());
    Real current_root_, low, high;   
    
    // # "Administration"

    // Gets the line of table that corresponds to order (n),     
    RealVector& alphaTable( this->alphaTable[n] );		        
    
    // Gets it's size 
    const RealVector::size_type oldSize( alphaTable.size() );	
    
    // Resizes and fills with zeroes if the vector is too small
    if( i >= oldSize )
    {
        alphaTable.resize( i+1, 0 );	
    }

    // # Calculating the root

    // Only calculates the root is this has not been done already
    if (alphaTable[i] == 0)		
    {        
    
        // Method 1. SCANNING. If the roots are not expected to lie close enough 
        // to the estimate, use the "scanning" procedure to find an interval 
        // that contains a root. (More robust method.)
        //      If it is established that we can use method 2, 
        // alpha_x_scan_table_[n] will contain a value < 0.
	    if (alpha_x_scan_table_[n] >= 0)
	    {		        
	        // ### Gets estimate of interval by sign-change-searching
	        //      high and low are the return-values of GiveRootInterval.
	        GiveRootInterval(low, high, n);	        
	        
	        // ### Finds the root using the GSL rootfinder
            current_root_ = getAlphaRoot(low, high, n);
            
            // ### Puts the found root in the table.
		    alphaTable[i]= current_root_;
		    		    
		    // Check if we can use method 2 for next roots
		    decideOnMethod2(n, i);		    
	    }
	    // Method 2. ASSUMING ROOTS AT ~FIXED INTERVAL. If next root is expected 
	    // to lie at distance close enough to estimated distance the root finder 
	    // can be called using a simple estimated interval.
	    else 
	    {
            // ### Get interval by simple extrapolation
	        GiveRootIntervalSimple(low, high, n, i);	        
	    
	        // ### Finds the root using the GSL rootfinder
            current_root_ = getAlphaRoot(low, high, n);

            // ### Puts the found root in the table.
		    alphaTable[i] = current_root_;		    
	    }
    }
    
    return alphaTable[i];

}

// =============================================================================
// End modification
// =============================================================================

// calculates the ith term with exponent and time for the survival probability
const Real 
GreensFunction2DRadAbs::p_survival_i_exp_table( const unsigned int i,
							const Real t,
							const RealVector& table ) const
{
    const Real alpha( this->getAlpha( 0, i ) );
    return std::exp( - getD() * t * alpha * alpha ) * table[i];
}

// adds the exponential with the time to the sum. Needed for the calculation of the flux throught the outer
// interface
const Real 
GreensFunction2DRadAbs::leavea_i_exp( const unsigned int i,
					      const Real t) const
{
    const Real alpha( this->getAlpha( 0,i ) );
    return std::exp( - getD() * t * alpha * alpha ) * calc_A_i_0( alpha );
}

// adds the exponential with the time to the sum. Needed for the inner interface (reaction)
const Real 
GreensFunction2DRadAbs::leaves_i_exp( const unsigned int i,
					      const Real t) const
{
    const Real alpha( this->getAlpha( 0, i ) );

    return std::exp( - getD() * t * alpha * alpha ) * leaves_i( alpha );
}


// calculates the Bossen function for a given r
const Real
GreensFunction2DRadAbs::p_int_r_i_exp_table( const unsigned int i,
				                     const Real r,
                				     const RealVector& Y0_aAnTable,
               					     const RealVector& J0_aAnTable,
							const RealVector& Y0J1J0Y1Table ) const
{
	const Real alpha( this->getAlpha( 0, i ) );	// get the root An
	const Real r_An( r*alpha);
   
	const Real J1_rAn (gsl_sf_bessel_J1(r_An));
	const Real Y1_rAn (gsl_sf_bessel_Y1(r_An));

	const Real result (Y0_aAnTable[i]*r*J1_rAn - J0_aAnTable[i]*r*Y1_rAn - Y0J1J0Y1Table[i]);
	return result;
}


// This tries to guess the maximum number of n iterations it needs for calculating the survival probability
// Not really sure yet how this works

const unsigned int
GreensFunction2DRadAbs::guess_maxi( const Real t ) const
{
    const unsigned int safety( 2 );

    if( t >= INFINITY )
    {
        return safety;
    }

    const Real D( getD() );
    const Real sigma( getSigma() );
    const Real a( geta() );

    const Real alpha0( getAlpha( 0, 0 ) );
    const Real Dt( D * t );
    const Real thr( ( exp( - Dt * alpha0 * alpha0 ) / alpha0 ) * this->EPSILON * 1e-1 );
    const Real thrsq( thr * thr );

    if( thrsq <= 0.0 )
    {
        return this->MAX_ALPHA_SEQ;
    }

    const Real max_alpha( 1.0 / ( sqrt( exp( gsl_sf_lambert_W0( 2 * Dt / thrsq ) ) * thrsq ) ) );
    const unsigned int maxi( safety + static_cast<unsigned int>( max_alpha * ( a - sigma ) / M_PI ) );

    return std::min( maxi, this->MAX_ALPHA_SEQ );
}


// Calculates the survival probability at a given time.
// This is a little wrapper for the p_survival_table so that you can easily calculate the survival probability
// at a given time
const Real 
GreensFunction2DRadAbs::p_survival( const Real t) const
{
    RealVector psurvTable;

    const Real p( p_survival_table( t, psurvTable ) );

    return p;
}

// This actually calculates the Survival probability at time t given the particle was at r0 at time 0
// It uses the pSurvTable for efficiency (so you don't have to calculate all the constant factors all
// the time)
const Real 
GreensFunction2DRadAbs::p_survival_table( const Real t,
		  				  RealVector& psurvTable ) const
{
	Real p;
	const unsigned int maxi( guess_maxi( t ) );	// guess the maximum number of iterations required
            
        if( psurvTable.size() < maxi + 1 )		// if the dimensions are good then this means
        {						// that the table is filled
        	IGNORE_RETURN getAlpha( 0, maxi );	// this updates the table of roots
                this->createPsurvTable( psurvTable);	// then the table is filled with data
        }
        p = funcSum_all( boost::bind( &GreensFunction2DRadAbs::p_survival_i_exp_table,
                                          this,
                                          _1, t, psurvTable ),
                         maxi );	// calculate the sum at time t

        return p*M_PI*M_PI_2;
}

// calculates the flux leaving through the inner interface at a given moment
// FIXME: This is inaccurate for small t's!!
const Real 
GreensFunction2DRadAbs::leaves( const Real t) const
{
    const Real sigma(this->getSigma());
    const Real D(this->getD() );

    const Real p( funcSum( boost::bind( &GreensFunction2DRadAbs::leaves_i_exp,
                                        this,
                                        _1, t),
                           this->MAX_ALPHA_SEQ ) );

    return M_PI_2*M_PI*D*sigma*p;	// The minus is not there because the flux is in the negative r
					// direction, and the direction is already accounted for in the derivative of B0,n(r)
					// See also leaves_i
}

// calculates the flux leaving through the outer interface at a given moment
const Real 
GreensFunction2DRadAbs::leavea( const Real t) const
{
    const Real D(this->getD() );

    const Real p( funcSum( boost::bind( &GreensFunction2DRadAbs::leavea_i_exp,
                                        this,
                                        _1, t),
                           this->MAX_ALPHA_SEQ ) );
    return M_PI*D*p;
}


// calculates the sum of the sequence for drawR based upon the values in the tables and r
const Real
GreensFunction2DRadAbs::p_int_r_table( const Real r,
               					const RealVector& Y0_aAnTable,
               					const RealVector& J0_aAnTable,
               					const RealVector& Y0J1J0Y1Table ) const
{
    const Real p( funcSum( boost::bind( &GreensFunction2DRadAbs::p_int_r_i_exp_table,
                                        this,
                                        _1, r, Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table ),
                           Y0_aAnTable.size() ) );
    return p*M_PI*M_PI_2;
}

// Used by drawTime
// Wrapper for p_survival_table for the interator to find the root for drawTime
const Real
GreensFunction2DRadAbs::p_survival_table_F( const Real t,
                    					const p_survival_table_params* params )
{
    const GreensFunction2DRadAbs* const gf( params->gf ); // the current gf (not sure why this is here)
    RealVector& table( params->table );		// table is empty but will be filled in p_survival_table
    const Real rnd( params->rnd );

    return rnd - gf->p_survival_table( t, table );
}


// a wrapper to make p_int_r_table available to the iterator calculating the root
const Real
GreensFunction2DRadAbs::p_int_r_F( const Real r,
					   const p_int_r_params* params )
{
    const GreensFunction2DRadAbs* const gf( params->gf ); 
    const RealVector& Y0_aAnTable( params->Y0_aAnTable );
    const RealVector& J0_aAnTable( params->J0_aAnTable );
    const RealVector& Y0J1J0Y1Table( params->Y0J1J0Y1Table );
    const Real rnd( params->rnd );

    return gf->p_int_r_table( r, Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table) - rnd;
}


// Draws a first passage time, this could be an escape (through the outer boundary) or a reaction (through
// the inner boundary)
Real GreensFunction2DRadAbs::drawTime( const Real rnd) const
{
    const Real D( this->getD() );
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );
    const Real kf( this->getkf() );
    const Real r0( this->getr0() );

    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
    THROW_UNLESS( std::invalid_argument, sigma <= r0 && r0 <= a );

    Real dist;


    if( r0 == a || a == sigma )		// when the particle is at the border or if the PD has no real size
    {
	return 0.0;
    }

    if ( kf == 0.0 )			// if there was only one absorbing boundary
    {	
        dist = a - r0;
    }
    else
    {	
        dist = std::min( a - r0, r0 - sigma );	// take the shortest distance to a boundary
    }
    Real t_guess = dist * dist / ( 4.0 * D );	// get some initial guess for the time, dr=sqrt(2dDt) with d
						// the dimensionality (2 in this case)


	const Real minT( std::min( sigma * sigma / D * this->MIN_T_FACTOR,
                               t_guess * 1e-7 ) );	// something with determining the lowest possible t


	RealVector psurvTable;		// this is still empty as of now->not used
	p_survival_table_params params = { this, psurvTable, rnd };
	gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_survival_table_F ),
	    &params 
	};

	// put in a upper and lower limit (the picked time cannot be infinite!)
	Real low( t_guess );	
	Real high( t_guess );

	// adjust high and low to make sure that f( low ) and f( high ) straddle.
	Real value( GSL_FN_EVAL( &F, t_guess ) );

	if( value < 0.0 )			// if the function is below zero at the guess the upper
	{					// boundary should be moved (passed the zero point)
		do
		{
			high *= 10;
			value = GSL_FN_EVAL( &F, high );
            
			if( fabs( high ) >= 1e10 )	// if high time is way too high forget about it
			{
				std::cerr << "Couldn't adjust high. F(" << high <<
					") = " << GSL_FN_EVAL( &F, high ) << "; r0 = " << r0 << 
					", " << dump() << std::endl;
				throw std::exception();
			}
		}
		while ( value < 0.0 );
	}
	else				// if the function is over zero (or at zero!) then the lower
	{					// boundary should be moved
                Real value_prev( value );
                do
                {       low *= .1;      // keep decreasing the lower boundary until the function straddles
                        value = GSL_FN_EVAL( &F, low );     // get the accompanying value

                        if( fabs( low ) <= minT || fabs( value - value_prev ) < EPSILON*this->T_TYPICAL )
                        {
                                std::cerr << "Couldn't adjust low. F(" << low <<
                                        ") = " << value << std::endl;
                                return low;
                        }
                        value_prev = value;
                }
                while ( value >= 0.0 );
    }

    // find the intersection of the cummulative survival probability and the randomly generated number
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );	// initialize the solver
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    const Real t( findRoot( F, solver, low, high,		// find the intersection between the random
		EPSILON*T_TYPICAL, EPSILON, "GreensFunction2DRadAbs::drawTime" ) );
								// number and the cumm probability
    gsl_root_fsolver_free( solver );

    return t;
}

// This determines based on the flux at a certain time, if the 'escape' was a reaction or a proper escape
GreensFunction2DRadAbs::EventKind
GreensFunction2DRadAbs::drawEventType( const Real rnd, 
					       const Real t ) const
{
    const Real D( this->getD() );
    const Real sigma( this->getSigma() );
    const Real kf( this->getkf() );
    const Real a( this->geta() );
    const Real r0( this->getr0() );

    THROW_UNLESS( std::invalid_argument, 0 <= rnd && rnd < 1.0 );
    THROW_UNLESS( std::invalid_argument, sigma <= r0 && r0 < a );
    THROW_UNLESS( std::invalid_argument, t > 0.0 );

    if( kf == 0.0 )	// if there cannot be any flow through the radiating boundary it is always an escape
    {
        return IV_ESCAPE;
    }
    
    // First, check if r0 is close only either to a or sigma relative
    // to Dt.  In such cases, the event type is always ESCAPE or REACTION,
    // respectively.   This avoids numerical instability in calculating
    // leavea() and/or leaves().

    // Here, use a rather large threshold for safety.
    const unsigned int H( 6 ); 				// 6 times the msd travelled as threshold
    const Real max_dist( H * sqrt( 4.0 * D * t ) );
    const Real a_dist( a - r0 );
    const Real s_dist( r0 - sigma );


    if( a_dist > max_dist )
    {
        if( s_dist < max_dist )
        {
            return IV_REACTION;
        }
    }
    else // a_dist < max_dist
    {
        if( s_dist > max_dist )
        {
            return IV_ESCAPE;
        }
    }

    const Real reaction( leaves( t ) );	// flux through rad boundary
    const Real escape( leavea( t ) );	// flux through abs boundary
    const Real value( reaction / ( reaction + escape ) );

    if( rnd <= value )  
    {
	return IV_REACTION;   // leaves -> return 0
    }
    else 
    {
	return IV_ESCAPE;     // leavea -> return 1
    }
}

// This draws a radius R at a given time, provided that the particle was at r0 at t=0
Real GreensFunction2DRadAbs::drawR( const Real rnd, 
						  const Real t ) const
{
    const Real D( this->getD() );
    const Real sigma( getSigma() );
    const Real a( this->geta() );
    const Real r0( this->getr0() );

    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 < a );

    if( t == 0.0 )		// if no time has passed
    {
	return r0;
    }

    const Real psurv( p_survival( t ) );	// calculate the survival probability at this time
						// this is used as the normalization factor
						// BEWARE!!! This also produces the roots An and therefore
						// SETS THE HIGHEST INDEX -> side effect
						// VERY BAD PROGRAMMING PRACTICE!!

    RealVector Y0_aAnTable;
    RealVector J0_aAnTable;
    RealVector Y0J1J0Y1Table;
    createY0J0Tables( Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table, t);

    p_int_r_params params = { this, t, Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table, rnd * psurv };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_int_r_F ),
	    &params 
	};

    // adjust low and high starting from r0.
    // this is necessary to avoid root finding in the long tails where
    // numerics can be unstable.
    Real low( r0 );				// start with the initial position as the first guess
    Real high( r0 );
    Real value (0);
    unsigned int H( 3 );

    const Real msd( sqrt( 4.0 * D * t ) );
    if( GSL_FN_EVAL( &F, r0 ) < 0.0 )
    {
	do
        {
            high = r0 + H * msd;
            if( high > a )
            {
                if( GSL_FN_EVAL( &F, a ) < 0.0 )	// something is very wrong, this should never happen
                {
                    printf( "drawR: p_int_r_table( a ) < 0.0. returning a.\n" );
                    return a;
                }
                high = a;
                break;
            }
            value = GSL_FN_EVAL( &F, high );
            ++H;
        }
	while (value < 0.0);

    }
    else
    {
	do
        {
            low = r0 - H * msd;
            if( low < sigma )
            {
                if( GSL_FN_EVAL( &F, sigma ) > 0.0 )
                {
                    printf( "drawR: p_int_r_table( sigma ) > 0.0. "
                            "returning sigma.\n" );
                    return sigma;
                }

                low = sigma;
                break;
            }

            value = GSL_FN_EVAL( &F, low );
            ++H;
        }
	while ( value > 0.0 );
    }


    // root finding by iteration.

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    const Real r( findRoot( F, solver, low, high,		// find the intersection between the random
		 L_TYPICAL*EPSILON, EPSILON, "GreensFunction2DRadAbs::drawR" ) );
								// number and the cumm probability
    gsl_root_fsolver_free( solver );

    return r;
}



// The calculates constant factor m,n for the drawing of theta. These factors are summed later.
const Real GreensFunction2DRadAbs::p_m_alpha( const unsigned int n,
						      const unsigned int m,
						      const Real r,
						      const Real t ) const
{
	const Real sigma( this->getSigma() );
	const Real h( this->geth() );
	const Real a( this->geta() );
	const Real D( this->getD() );
	const Real alpha( this->getAlpha( m, n ) ); // Gets the n-th root using the 
	                                            // besselfunctions of order m.
	const Real r0( this->getr0() );

	const Real alpha_sq( alpha * alpha );
	const Real realm( static_cast<Real>( m ) );
	const Real msq( realm * realm);
	const Real ssq( sigma * sigma);

	const Real s_Anm (sigma * alpha);
	const Real a_Anm (a * alpha);
	const Real r0Anm (r0 * alpha);
	const Real r_Anm (r * alpha);

//	const CylindricalBesselGenerator& s(CylindricalBesselGenerator::instance());

	// calculate the needed bessel functions
	const Real Jm_sAnm   (gsl_sf_bessel_Jn(m, s_Anm));
	const Real Jmp1_sAnm (gsl_sf_bessel_Jn(m+1, s_Anm));	// prime
//	const Real Jm_sAnm   (s.J(m, s_Anm));
//	const Real Jmp1_sAnm (s.J(m+1, s_Anm));	// prime

	const Real Jm_aAnm   (gsl_sf_bessel_Jn(m, a_Anm));
	const Real Ym_aAnm   (gsl_sf_bessel_Yn(m, a_Anm));
//	const Real Jm_aAnm   (s.J(m, a_Anm));
//	const Real Ym_aAnm   (s.Y(m, a_Anm));

	const Real Jm_r0Anm  (gsl_sf_bessel_Jn(m, r0Anm));
	const Real Ym_r0Anm  (gsl_sf_bessel_Yn(m, r0Anm));
//	const Real Jm_r0Anm  (s.J(m, r0Anm));
//	const Real Ym_r0Anm  (s.Y(m, r0Anm));

	const Real Jm_rAnm   (gsl_sf_bessel_Jn(m, r_Anm));
	const Real Ym_rAnm   (gsl_sf_bessel_Yn(m, r_Anm));
//	const Real Jm_rAnm   (s.J(m, r_Anm));
//	const Real Ym_rAnm   (s.Y(m, r_Anm));

	// calculating An,m
	const Real h_ma (h - realm/sigma);
	const Real rho (h_ma*Jm_sAnm + alpha*Jmp1_sAnm);
	const Real rho_sq (rho*rho);
	// calculating Bn,m(r')
	const Real B_n_m_r0 (Jm_r0Anm * Ym_aAnm  -  Ym_r0Anm * Jm_aAnm);

	const Real A_n_m ((alpha_sq * rho_sq * B_n_m_r0)/( rho_sq - (Jm_aAnm*Jm_aAnm)*(h*h + alpha_sq - msq/ssq)));


	// calculating Bn,m(r*)
	const Real B_n_m_r (Jm_rAnm * Ym_aAnm  -  Ym_rAnm * Jm_aAnm);


	// calculating the result
	const Real result( A_n_m * B_n_m_r * exp(-D*alpha_sq*t) );

	return result;
}


// This calculates the m-th constant factor for the drawTheta method. 
const Real 
GreensFunction2DRadAbs::p_m( const Integer m,
				     const Real r,
				     const Real t ) const
{
    const Real p( funcSum( boost::bind( &GreensFunction2DRadAbs::p_m_alpha,
					this,
					_1, m, r, t ),	// The m-th factor is a summation over n
			   MAX_ALPHA_SEQ, EPSILON ) );
    return p;
}

// this should make the table of constants used in the iteration for finding the root for drawTheta
// The index of the array is consistent with the index of the summation
void
GreensFunction2DRadAbs::makep_mTable( RealVector& p_mTable,
					      const Real r, 
					      const Real t ) const
{
	p_mTable.clear();

	const Real p_0 ( this->p_m( 0, r, t ) );	// This is the p_m where m is 0, for the denominator
	p_mTable.push_back( p_0 );			// put it in the table

	const Real p_1 ( this->p_m( 1, r, t ) / p_0 );
	p_mTable.push_back( p_1 );			// put the first result in the table

	if( p_1 == 0 )
	{
		return;					// apparantly all the terms are zero? We are finished
		                        // TODO: is this assumption correct??
	}

	const Real threshold( fabs( EPSILON * p_1  ) );	// get a measure for the allowed error, is this correct?

	Real p_m_abs (fabs (p_1));
	Real p_m_prev_abs;	
	unsigned int m( 1 );
	do
	{
		m++;
		if( m >= this->MAX_ORDER )		// If the number of terms is too large
		{
			std::cerr << "p_m didn't converge (m=" << m << "), continuing..." << std::endl;
			break;
		}

		p_m_prev_abs = p_m_abs;					// store the previous term
		const Real p_m( this->p_m( m, r, t ) / p_0 );	// get the next term

		if( ! std::isfinite( p_m ) )		// if the calculated value is not valid->exit
		{
			std::cerr << "makep_mTable: invalid value; " <<
				p_m << "( m= " << m << ")." << std::endl;
			break;
		}

		p_mTable.push_back( p_m );				// put the result in the table
		p_m_abs = fabs( p_m );					// take the absolute value
        }
	while (p_m_abs >= threshold || p_m_prev_abs >= threshold || p_m_abs >= p_m_prev_abs );

	// truncate when converged enough.
	// if the current term is smaller than threshold
	// AND the previous term is also smaller than threshold
	// AND the current term is smaller than the previous
		
}

// This method calculates the constants for the drawTheta method when the particle is at the boundary
const Real 
GreensFunction2DRadAbs::dp_m_alpha_at_a( const unsigned int n,
						 const unsigned int m,
						 const Real t ) const
{
        const Real sigma( this->getSigma() );
        const Real h( this->geth() );
        const Real a( this->geta() );
        const Real D( this->getD() );
        const Real r0( this->getr0() );

        const Real alpha( this->getAlpha( m, n ) ); // get the n-th root using the besselfunctions of order m
       
        const Real alpha_sq( alpha * alpha );
        const Real realm( static_cast<Real>( m ) );
        const Real msq( realm * realm);
        const Real ssq( sigma * sigma);

        const Real s_Anm (sigma*alpha);
        const Real a_Anm (a*alpha);
        const Real r0Anm (r0*alpha);

        const Real Jm_sAnm   (gsl_sf_bessel_Jn(m, s_Anm));
        const Real Jmp1_sAnm (gsl_sf_bessel_Jn(m+1, s_Anm));    
        const Real Jm_aAnm   (gsl_sf_bessel_Jn(m, a_Anm));
        const Real Ym_aAnm   (gsl_sf_bessel_Yn(m, a_Anm)); 

        const Real Jm_r0Anm  (gsl_sf_bessel_Jn(m, r0Anm));
        const Real Ym_r0Anm  (gsl_sf_bessel_Yn(m, r0Anm)); 

        // calculating An,m
        const Real h_ma (h - realm/sigma);
        const Real rho (h_ma*Jm_sAnm + alpha*Jmp1_sAnm);
        const Real rho_sq (rho*rho);

        // calculating Bn,m(r')
        const Real B_n_m_r0 (Jm_r0Anm * Ym_aAnm  -  Ym_r0Anm * Jm_aAnm);

        const Real A_n_m ((alpha_sq * rho_sq * B_n_m_r0)/( rho_sq - (Jm_aAnm*Jm_aAnm)*(h*h + alpha_sq - msq/ssq)));

        // calculating the result
        const Real result( A_n_m * exp(-D*alpha_sq*t) );

	return result;
}

// Makes the sum over n for order m for the constants for the drawtheta Method
const Real 
GreensFunction2DRadAbs::dp_m_at_a( const Integer m,
					   const Real t ) const
{

    const Real p( funcSum( 
                        boost::bind( 
                            &GreensFunction2DRadAbs::dp_m_alpha_at_a,
			                    this,
			                    _1, 
                                m,
                                t 
                        ),
		                MAX_ALPHA_SEQ, 
                        EPSILON 
                  ) 
                );

               // boost::bind 
               // explanation by example:
               // "bind(f, _1, 5)(x)"  is equivalent to "f(x, 5)"
               // this means funcsum receives f(x, m=.., t=..) as
               // input, with m and t already determined.
    
                // Arguments of dp_m_alpha_at_a:
                // n, m, t

    return p;
}

// creates a tables of constants for drawTheta when the particle is at the edge of the domain
void
GreensFunction2DRadAbs::makedp_m_at_aTable( RealVector& p_mTable,
						    const Real t ) const
{
    p_mTable.clear();

    const Real p_0 ( this->dp_m_at_a( 0, t ) );	// This is the p_m where m is 0, for the denominator
    p_mTable.push_back( p_0 );                      // put it in the table

    const Real p_1 ( this->dp_m_at_a( 1, t ) / p_0 );
    p_mTable.push_back( p_1 );                      // put the first result in the table


    if( p_1 == 0 )
    {
            return;                 // apparantly all the terms are zero? We are finished
    }

    const Real threshold( fabs( EPSILON * p_1  ) );	// get a measure for the allowed error

    Real p_m_abs (fabs (p_1));
    Real p_m_prev_abs;      
    unsigned int m( 1 );
    do
    {
            m++;
            if( m >= this->MAX_ORDER )				// If the number of terms is too large
            {
                    std::cerr << "dp_m didn't converge (m=" << m << "), continuing..." << std::endl;
                    break;
            }


            p_m_prev_abs = p_m_abs;					// store the previous term
            const Real p_m( this->dp_m_at_a( m, t ) / p_0 );	// get the next term

            // DEBUG (something to check in the future?)
            // if (p_m_abs == 0) {
            //    std::cerr << "Zero valued term found, but convergence is:" <<
            //        p_mTable[p_mTable.size()-1-1]/p_mTable[p_mTable.size()-2-1];
            // }
            // END DEBUG

            if( ! std::isfinite( p_m ) )			// if the calculated value is not valid->exit
            {
                    std::cerr << "makedp_m_at_aTable: invalid value; " <<
                            p_m << "( m= " << m << ")." << "( t= " << t << ")."<< "( p_0= " << p_0 << ")." << std::endl;
                    break;
            }

            p_mTable.push_back( p_m );                              // put the result in the table
            p_m_abs = fabs( p_m );                                  // take the absolute value
    }            
    while (p_m_abs >= threshold || p_m_prev_abs >= threshold || p_m_abs >= p_m_prev_abs );
    // truncate when converged enough.
    // if the current term is smaller than threshold
    // AND the previous term is also smaller than threshold
    // AND the current term is smaller than the previous
}

// This calculates the m-th term of the summation for the drawTheta calculation
// Note that m here starts at 0 and in the equations the sum starts at 1!
const Real 
GreensFunction2DRadAbs::ip_theta_n( const unsigned int m,
						const Real theta,
						const RealVector& p_nTable) const
{
	const unsigned int m_p1 (m+1);		// artificial increase of m to make sure m starts at 1
        return p_nTable[m_p1] * sin (m_p1*theta)/m_p1;
}


// calculates the cummulative probability of finding the particle at a certain theta
// It is used by the drawTheta method
// It uses the p_nTable for it to speed things up
const Real 
GreensFunction2DRadAbs::ip_theta_table( const Real theta,
						const RealVector& p_nTable ) const
{
    const unsigned int maxm( p_nTable.size()-1 );	// get the length of the sum
							// it is shifted one because the first entry should
							// be used (m=0)

    const Real p( funcSum_all( boost::bind( &GreensFunction2DRadAbs::ip_theta_n,
                                             this,
                                             _1, theta, p_nTable ),
                                maxm ) );
    return p;
}

// function to iterate when drawing the theta
const Real
GreensFunction2DRadAbs::ip_theta_F( const Real theta,
					    const ip_theta_params* params )
{
    const GreensFunction2DRadAbs* const gf( params->gf ); 
    const RealVector& p_nTable( params->p_nTable );	// table with useful constants
    const Real value( params->value );

    return theta/(M_PI*2) + (gf->ip_theta_table( theta, p_nTable )/M_PI) - value;
}


// This method draws a theta given a certain r and time (and intial condition of course)
Real 
GreensFunction2DRadAbs::drawTheta( const Real rnd,
					   const Real r, 
					   const Real t ) const
{
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );
    const Real D( this->getD() );
    const Real r0( this->getr0() );

	// input parameter range checks.
	THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
	THROW_UNLESS( std::invalid_argument, sigma <= r0 && r0 <= a );
	THROW_UNLESS( std::invalid_argument, sigma <= r && r <= a);
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	// t == 0 means no move.
	if( t <= T_TYPICAL*EPSILON || D == 0 || fabs(r0 - a) <= EPSILON*L_TYPICAL || rnd <= EPSILON)
	{
		return 0.0;
	}
	else if (r == sigma)	// a reaction has occured, the angle is irrelevant
	{
		return 0.0;
	}

	// making the tables with constants

	RealVector p_mTable;			// a table with constants to make calculations much faster
	if( fabs(r - a) <= EPSILON*L_TYPICAL )	// If the r is at the outer boundary
	{
		makedp_m_at_aTable( p_mTable, t );	// making the table if particle on the outer boundary
	}
	else
	{
		makep_mTable( p_mTable, r, t );	// making the table of constants for the regular case
	}

	// preparing the function
	
	// ip_theta_params is a struct.
	ip_theta_params params = { this, r, t, p_mTable, rnd*0.5 };	// r, r0, t are not required	    
	                                                        // 0.5 is not even used
    // F is a struct of type "gsl_function" that contains a function pointer
	// and the required parameters.
	gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &ip_theta_F ), // ip_theta_F is 
	                                                    // theta pdf function.
	    &params 
        // reinterpret_cast converts any pointer type to any other pointer type.
        // reinterpret_cast<new_type> variable
	}; 
	    
	// finding the root
	const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
	gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
	const Real theta( findRoot( F, solver, 0, M_PI, EPSILON, EPSILON,
                "GreensFunction2DRadAbs::drawTheta" ) );
	gsl_root_fsolver_free( solver );

	return theta;
}


//
// debug
//

// Debug functionality
std::string GreensFunction2DRadAbs::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << ", sigma = " << this->getSigma() <<
	", a = " << this->geta() <<
	", kf = " << this->getkf() <<
	", r0 = " << this->getr0() <<
	", h = " << this->geth() << std::endl;
    return ss.str();
}    

// Debug functionality
// Directly outputs probability distribution function value of leaving angle
// for given theta, r and t.
Real 
GreensFunction2DRadAbs::givePDF( const Real theta,
					   const Real r, 
					   const Real t ) const
{
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );
    const Real D( this->getD() );
    const Real r0( this->getr0() );

	// input parameter range checks.
	THROW_UNLESS( std::invalid_argument, sigma <= r0 && r0 <= a );
	THROW_UNLESS( std::invalid_argument, sigma <= r && r <= a);
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	// making the tables with constants
	RealVector p_mTable;			// a table with constants to make calculations much faster
	if( fabs(r - a) <= EPSILON*L_TYPICAL )	// If the r is at the outer boundary
	{
		makedp_m_at_aTable( p_mTable, t );	// making the table if particle on the outer boundary
	}
	else
	{
		makep_mTable( p_mTable, r, t );	// making the table of constants for the regular case
	}

    // Return the pdf
  	ip_theta_params params = { this, r, t, p_mTable, 0.0 };	// r, r0, t are not required	    
    Real PDF (ip_theta_F( theta, &params ));
    return PDF;
	
}

// It is used by the drawTheta method
// It uses the p_nTable for it to speed things up
/*
const Real 
GreensFunction2DRadAbs::debug_ip_theta_table( const Real theta) const
{
    const RealVector& p_nTable( params->p_nTable );	// table with useful constants
    
    const unsigned int maxm( p_nTable.size()-1 );	// get the length of the sum
							// it is shifted one because the first entry should
							// be used (m=0)

    const Real p( funcSum_all( boost::bind( &GreensFunction2DRadAbs::ip_theta_n,
                                             this,
                                             _1, theta, p_nTable ),
                                maxm ) );
    return p;
}
*/
