#ifndef AP4_ANALYTICAL
#define AP4_ANALYTICAL

double psi_of_eps_ap4(double eps);
double p_of_rho_analytical_ap4(double rho);
double rho_of_p_analytical_ap4(double p);
double dp_drho_ap4(double rho,double pressy);
double drho_dp_analytical_ap4(double pressy);
double ddp_drhorho_ap4(double rho,double pressy);
double ddrho_dPP_analytical_ap4(double pressy);
double dpsi_deps_ap4(double eps);
double deps_dpsi_ap4(double psi,double rho);
double ddeps_dpsipsi_ap4(double psi,double rho);
double ddpsi_depseps_ap4(double eps);

#endif // AP4_ANALYTICAL

