#ifndef FREE_FUNTIONS_HPP
#define FREE_FUNTIONS_HPP

#include "Defs.hpp"

Real expxsq_erfc(Real x);

Real W(Real a, Real b);

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

Real g_bd_3D(Real r0, Real sigma, Real t, Real D);

Real I_bd_3D(Real sigma, Real t, Real D);

Real I_bd_r_3D(Real r, Real sigma, Real t, Real D);

Real drawR_gbd_3D(Real rnd, Real sigma, Real t, Real D);

Real g_bd_1D(Real r0, Real sigma, Real t, Real D, Real v);

Real I_bd_1D(Real sigma, Real t, Real D, Real v);

Real I_bd_r_1D(Real r, Real sigma, Real t, Real D, Real v);

Real drawR_gbd_1D(Real rnd, Real sigma, Real t, Real D, Real v);

#endif /* FREE_FUNTIONS_HPP */
