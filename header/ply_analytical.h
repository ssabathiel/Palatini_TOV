#ifndef PLY_ANALYTICAL
#define PLY_ANALYTICAL



double Kf(double x);
double Power(double base, double exp);
double Log(double x);

double p_of_rho_analytical_ply(double rho);
double rho_of_p_analytical_ply(double p);

double dp_drho_ply(double rho,double pressy);
double drho_dp_analytical_ply(double pressy);
double ddp_drhorho_ply(double rho,double pressy);
double ddrho_dPP_analytical_ply(double pressy);

double psi_of_eps_ply(double eps);

double dpsi_deps_ply(double eps);
double deps_dpsi_ply(double psi,double rho);
double ddeps_dpsipsi_ply(double psi,double rho);
double ddpsi_depseps_ply(double eps);





#endif // PLY_ANALYTICAL

