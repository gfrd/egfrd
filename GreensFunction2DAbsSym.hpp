#if !defined( __FIRSTPASSAGEGREENSFUNCTION2D_HPP )
#define __FIRSTPASSAGEGREENSFUNCTION2D_HPP

#include <vector>
#include <boost/multi_array.hpp>

#include "Defs.hpp"


class GreensFunction2DAbsSym
{

public:

    GreensFunction2DAbsSym( const Real D, const Real a )
	:
	D( D ),
        a( a )
    {
	;
    }

    virtual ~GreensFunction2DAbsSym()
    {
	;
    }

    const char* getName() const
    {
        return "GreensFunction2DAbsSym";
    }

    const Real getD() const
    {
	return this->D;
    }

    const Real geta() const
    {
        return this->a;
    }

    const Real p_survival( const Real t ) const; 

    const Real drawTime( const Real rnd ) const;

    const Real drawR( const Real rnd, const Real t ) const;

    const Real p_int_r( const Real r, const Real t ) const;
    const Real p_int_r_free( const Real r, const Real t ) const;

    const std::string dump() const;

// private methods
private:

    struct p_survival_params
    {
	const GreensFunction2DAbsSym* const gf;
	const Real rnd;
    };

    static const Real p_survival_F( const Real t, 
				    const p_survival_params* params );

    struct p_r_params
    {
	const GreensFunction2DAbsSym* const gf;
	const Real t;
	const Real target;
    };

    static const Real p_r_free_F( const Real r, 
                                  const p_r_params* params );

    static const Real p_r_F( const Real r, 
			     const p_r_params* params );

// private variables
private:

    static const Real CUTOFF = 1e-10;

    // H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    // 5.6: ~1e-8, 6.0: ~1e-9
    static const Real CUTOFF_H = 6.0;

    const Real D;

    const Real a;
};



#endif // __PAIRGREENSFUNCTION2D_HPP
