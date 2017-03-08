// Greens function class for 2d Green's Function for 2d annulus with radial and
// axial dependence. Inner boundary is radiative (rad) (reaction event), outer 
// boundary is absorbing (abs) (escape event). Different "draw" functions 
// provide a way to draw certain values from the Green's Function, e.g. an
// escape angle theta ("drawTheta" function).
// 
// Written by Laurens Bossen. Adapted by Martijn Wehrens.
// FOM Institute AMOLF.

#if !defined( __FIRSTPASSAGEPAIRGREENSFUNCTION2D_HPP )
#define __FIRSTPASSAGEPAIRGREENSFUNCTION2D_HPP 

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>
#include <boost/array.hpp>

#include <gsl/gsl_roots.h>
#include "Logger.hpp"

#include "PairGreensFunction.hpp"


class GreensFunction2DRadAbs
    :
    public PairGreensFunction
{

public:
    // Defines vector from template defined in standard template library,
    // gives it's membervalues type Real. std:: clarifies the std namespace.
    typedef std::vector<Real> RealVector;

private:
    // Error tolerance used by default.
    static const Real TOLERANCE = 1e-8;

    // SphericalBesselGenerator's accuracy, used by some
    // theta-related calculations.

    static const Real MIN_T_FACTOR = 1e-8;

    static const Real L_TYPICAL = 1E-7; // typical length scale
    static const Real T_TYPICAL = 1E-5; // typical time scale
    static const Real EPSILON   = 1E-12; // relative numeric error  // TESTING temporarily increased; was 1e-12

    // DEFAULT = 30
    static const unsigned int MAX_ORDER = 30;         // The maximum number of m 
                                                      // terms
    static const unsigned int MAX_ALPHA_SEQ = 500;    // The maximum number of n 
                                                      // terms
    
    // Parameters for alpha-root finding
    // ======
    // See getAlpha() in cpp file for more information. 
    //
    // Parameters for scanning method
    // Left boundary of 1st search interval 1st root
    static const Real SCAN_START = 0.001;     
    // Length of the scanning interval relative to estimated interval
    static const Real FRACTION_SCAN_INTERVAL = .5; // TODO CHANGED THIS FROM .5 to .2
    
    // Other paramters
    // After CONVERGENCE_ASSUMED subsequent roots that lay within +/- 
    // INTERVAL_MARGIN from the distance to which the distance is known to 
    // converge, it is assumed all following roots have a distances inbetween
    // that don't deviate for more than INTERVAL_MARGIN from the distance to 
    // which the roots are known to converge (Pi/(a-sigma)).
    static const Real CONVERGENCE_ASSUMED = 25;
    static const Real INTERVAL_MARGIN = .33; 



public:

    const char* getName() const
    {
        return "GreensFunction2DRadAbs";
    }
    
    GreensFunction2DRadAbs( const Real D, 
				    const Real kf,
               			    const Real r0, 
				    const Real Sigma,
		                    const Real a );
  
    virtual ~GreensFunction2DRadAbs();

    const Real geth() const
    {
	return this->h;
    }

    const Real geta() const
    {
	return this->a;
    }

    const Real getestimated_alpha_root_distance_() const
    {
    return this->estimated_alpha_root_distance_;
    }

    virtual Real drawTime( const Real rnd) const;

    virtual EventKind drawEventType( const Real rnd, 
				   const Real t ) const;
    
    virtual Real drawR( const Real rnd, 
		      const Real t ) const;
    
    virtual Real drawTheta( const Real rnd,
			  const Real r, 
			  const Real t ) const;
      
    const Real f_alpha0( const Real alpha ) const;
  
    const Real f_alpha( const Real alpha, const Integer n ) const;

    
    const Real p_survival( const Real t) const;

    const Real p_survival_table( const Real t,
				 RealVector& table ) const;


    const Real leaves( const Real t) const;

    const Real leavea( const Real t) const;

    const Real p_m( const Integer n, const Real r, const Real t ) const;

    const Real dp_m_at_a( const Integer m, const Real t ) const;


    const Real p_m_alpha( const unsigned int n,
			  const unsigned int m,
			  const Real r, 
			  const Real t ) const;

    const Real dp_m_alpha_at_a( const unsigned int n,
				const unsigned int m,
				const Real t ) const;

    // methods below are kept public for debugging purpose.

    std::string dump() const;

                                                          
    const void GiveRootInterval(  Real& low,          
                                  Real& high,         
                                  const Integer n
                               ) const;
                                                                    

    const void GiveRootIntervalSimple(  Real& low, Real& high, 
                                        const Integer n, const Integer i
                                     ) const;

    const Real getAlphaRoot0( const Real low,
                              const Real high
                            ) const;
                                    
    const Real getAlphaRootN( const Real low,  
                              const Real high,  
                              const Integer n  
                            ) const;                                 
                                                            
    const Real getAlphaRoot(  const Real high, const Real low, const Integer n
                           ) const;
                                                    
    const void decideOnMethod2(size_t n, RealVector::size_type i) const;

    const void 
        needToSwitchBackMethod1( size_t n, RealVector::size_type i ) const;

    const Real getAlpha( size_t n, RealVector::size_type i ) const; 
    

    const Real p_survival_i( const Real alpha) const;

    const Real calc_A_i_0( const Real alpha) const;

    const Real leaves_i( const Real alpha) const;

    const boost::tuple<Real,Real,Real> Y0J0J1_constants ( const Real alpha,
                                                          const Real t) const;

//    const Real getAlpha( const size_t n, const RealVector::size_type i ) const;

//    const Real getAlpha0( const RealVector::size_type i ) const;

    Real givePDFTheta( const Real theta,
					   const Real r, 
					   const Real t ) const;    
					   
    Real givePDFR( const Real r, const Real t ) const;
    
    void dumpRoots( int n );

protected:

    void clearAlphaTable() const;


    RealVector& getAlphaTable( const size_t n ) const
    {
	return this->alphaTable[n];
    }

    const Real p_int_r_table( const Real r,
				const RealVector& Y0_aAnTable,
				const RealVector& J0_aAnTable,
				const RealVector& Y0J1J0Y1Table ) const;

    const Real ip_theta_table( const Real theta,
			       const RealVector& p_nTable ) const;

    const Real p_survival_i_exp_table( const unsigned int i,
				       const Real t,
				       const RealVector& table ) const;

    const Real leavea_i_exp( const unsigned int i,
			     const Real alpha) const;

    const Real leaves_i_exp( const unsigned int i,
			     const Real alpha) const;

    const Real ip_theta_n( const unsigned int m,
			   const Real theta,
			   const RealVector& p_nTable ) const;


    const Real p_int_r_i_exp_table( const unsigned int i,
					const Real r,
					const RealVector& Y0_aAnTable,
					const RealVector& J0_aAnTable,
					const RealVector& Y0J1J0Y1Table ) const;

    void createPsurvTable( RealVector& table) const; 

    void createY0J0Tables( RealVector& Y0_Table, RealVector& J0_Table, RealVector& Y0J1J0Y1_Table,
				const Real t ) const;

    void makep_mTable( RealVector& p_mTable,
		       const Real r, 
		       const Real t ) const;
    
    void makedp_m_at_aTable( RealVector& p_mTable,
			     const Real t ) const;

    const unsigned int guess_maxi( const Real t ) const;

    struct f_alpha0_aux_params
    {
        const GreensFunction2DRadAbs* const gf;
        const Real value;
    };

    static const Real
        f_alpha0_aux_F( const Real alpha,
                const f_alpha0_aux_params* const params );


    struct f_alpha_aux_params
    {
        const GreensFunction2DRadAbs* const gf;
        const Integer n;
        Real value;
    };

    static const Real 
    f_alpha_aux_F( const Real alpha,
		   const f_alpha_aux_params* const params );

    struct p_survival_table_params
    { 
	const GreensFunction2DRadAbs* const gf;
//	const Real r0;
	RealVector& table;
	const Real rnd;
    };

    static const Real 
    p_survival_table_F( const Real t,
                        const p_survival_table_params* const params );

    struct p_int_r_params
    { 
	const GreensFunction2DRadAbs* const gf;
	const Real t;
//	const Real r0;
	const RealVector& Y0_aAnTable;
	const RealVector& J0_aAnTable;
	const RealVector& Y0J1J0Y1Table;
	const Real rnd;
    };

    static const Real 
    p_int_r_F( const Real r,
	       const p_int_r_params* const params );

    struct ip_theta_params
    { 
	const GreensFunction2DRadAbs* const gf;
	const Real r;
//	const Real r0;
	const Real t;
	const RealVector& p_nTable;
	const Real value;
    };

    static const Real 
    ip_theta_F( const Real theta,
		const ip_theta_params* const params );

private:
    
    const Real h;
    const Real a;

    // Tables that hold calculated roots (y=0) of "alpha" function for each    
    // order n.
    mutable boost::array<RealVector,MAX_ORDER+1> alphaTable;

    // Constants used in the roots of f_alpha() finding algorithm.
    // ====
    //
    // This constant will simply be M_PI/(a-Sigma), the value to which the 
    // distance between roots of f_alpha() should converge.
    const Real estimated_alpha_root_distance_;        
    //
    // Table which tells us at which x we're left with scanning the alpha 
    // function for a sign change, for a given order n. (A sign change would 
    // indicate a root (y=0) lies between the boundaries of the "scanned" 
    // interval.) 
    //      If x_scan[n] < 0, this indicates scanning is no longer required
    // because the distance between the roots is converging and within 
    // boundaries that allow the direct use of the estimate interval width 
    // pi/(sigma-a).
    //      Initial values are set by constructor.
    mutable boost::array<Real,MAX_ORDER+1> alpha_x_scan_table_;     
    //
    // Table that keeps track of the number of previous subsequent roots that 
    // we're within margin of the distance to which they're expected to 
    // converge.
    mutable boost::array<int,MAX_ORDER+1> alpha_correctly_estimated_;
    
    static Logger& log_;

};



#endif // __FIRSTPASSAGEPAIRGREENSFUNCTION2D_HPP
