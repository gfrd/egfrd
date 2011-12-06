#ifndef FREE_FUNTIONS_HPP
#define FREE_FUNTIONS_HPP

#include "Defs.hpp"

Real expxsq_erfc(Real x);

Real W(Real a, Real b);

/* Functions needed in 1D Green's functions. For
   explanation, see the cpp file. */

Real XP00( Real r, Real t, Real r0, Real D, Real v );

Real XI00( Real r, Real t, Real r0, Real D, Real v );

Real XS00( Real t, Real r0, Real D, Real v );

Real XP10( Real r, Real t, Real r0, Real D, Real v );

Real XI10( Real r, Real t, Real r0, Real D, Real v );

Real XS10( Real t, Real r0, Real D, Real v );

Real XP20( Real r, Real t, Real r0, Real D, Real v );

Real XI20( Real r, Real t, Real r0, Real D, Real v );

Real XS20( Real t, Real r0, Real D, Real v );

Real XP30term_nov( Real r, Real t, Real r0, Real ka, Real D );

Real XP30term_v( Real r, Real t, Real r0, Real ka, Real D, Real v );

Real XP30( Real r, Real t, Real r0, Real ka, Real D, Real v );

Real XI30term_nov( Real r, Real t, Real r0, Real ka, Real D );

Real XI30( Real r, Real t, Real r0, Real ka, Real D, Real v );

Real XS30( Real t, Real r0, Real ka, Real D, Real v );

Real XP030( Real r, Real t, Real r0, Real ka, Real D );

Real XI030( Real r, Real t, Real r0, Real ka, Real D );

Real XS030( Real t, Real r0, Real ka, Real D );

/* 3D functions below. */

Real p_irr(Real r, Real t, Real r0, Real kf, Real D, Real sigma);

Real __p_irr(Real r, Real t, Real r0, Real kf, Real D, Real sigma, Real alpha);

Real p_free(Real r, Real r0, Real theta, Real t); 

Real p_survival_irr(Real t, Real r0, Real kf, Real D, Real sigma);

Real __p_reaction_irr(Real t, Real r0, Real kf, Real D, Real sigma, Real alpha, Real kD);

Real __p_reaction_irr_t_inf(Real r0, Real kf, Real sigma, Real kD);

Real p_survival_nocollision(Real t, Real r0, Real D, Real a);

Real dp_survival_nocollision(Real t, Real r0, Real D, Real a);

Real p_theta_free(Real theta, Real r, Real r0, Real t, Real D);

Real ip_theta_free(Real theta, Real r, Real r0, Real t, Real D);

/* Functions used in old Brownian Dynamic scheme. */

Real g_bd_3D(Real r0, Real sigma, Real t, Real D);

Real I_bd_3D(Real sigma, Real t, Real D);

Real I_bd_r_3D(Real r, Real sigma, Real t, Real D);

Real drawR_gbd_3D(Real rnd, Real sigma, Real t, Real D);

Real g_bd_1D(Real r0, Real sigma, Real t, Real D, Real v);

Real I_bd_1D(Real sigma, Real t, Real D, Real v);

Real I_bd_r_1D(Real r, Real sigma, Real t, Real D, Real v);

Real drawR_gbd_1D(Real rnd, Real sigma, Real t, Real D, Real v);

#endif /* FREE_FUNTIONS_HPP */
