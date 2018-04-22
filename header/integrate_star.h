#ifndef INTEGRATE_STAR
#define INTEGRATE_STAR

std::pair<double, double> tov_integrate(double rho_C);
std::pair<double, double> RK4_step(double mass_old, double press_old, double r, double dr);

#endif // INTEGRATE_STAR

