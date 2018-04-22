#ifndef EOS_ANALYTICAL
#define EOS_ANALYTICAL


double Kf(double x);
double Power(double base, double exp);
double Log(double x);

double p_of_rho_analytical(double rho);
double rho_of_p_analytical(double p);

double dp_drho(double rho,double pressy);
double drho_dp_analytical(double pressy);
double ddp_drhorho(double rho,double pressy);
double ddrho_dPP_analytical(double pressy);

double psi_of_eps(double eps);

double dpsi_deps(double eps);
double deps_dpsi(double psi,double rho);
double ddeps_dpsipsi(double psi,double rho);
double ddpsi_depseps(double eps);

#endif // EOS_ANALYTICAL

