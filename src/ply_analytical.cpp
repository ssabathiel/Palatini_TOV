
#include "header/ply_analytical.h"
#include <vector>
#include <math.h>
#include <cmath>
#include "linalg_library/spline.h"
#include <iomanip>
#include "linalg_library/interpolation.h"
#include "header/glob_variables.h"



tk::spline pf;       //pf = p(rho)-function
tk::spline rhof;     //rhof = rho(p)-function

/*
double E = exp(1);


double Kf(double x)
{
    double k_local;
    k_local = 1.0/(exp(x) + 1);
    return k_local;
}

*/

double psi_of_eps_ply(double eps)
{
    double psi = 2.0*eps + 5.29355;
    return psi;
}

double eps_of_psi_ply(double psi)
{
    double eps = (psi - 5.29355)/2.0;
    return eps;
}


double p_of_rho_analytical_ply(double rho)
{
    double eps=log10(rho);
    double psi = psi_of_eps_ply(eps);
    double pressy =  pow(10,psi); //exp(psi);          //--> should be 10^psi
    return pressy;
}





double rho_of_p_analytical_ply(double p)
{

    /*
    pf.set_points(rhos,press);
    rhof.set_points(press,rhos);
    rho = rhof(p_copy);
    */
    //return alglib::spline1dcalc(rho_of_press_alg, p);
    double psi=log10(p);
    double eps = eps_of_psi_ply(psi);
    double rho =  pow(10.0,eps); //exp(psi);          //--> should be 10^psi
    return rho;

}


/*
double Power(double base, double exp)
{
    return pow(base,exp);
}

double Log(double x)
{
    return log10(x);
}
*/
double dp_drho_ply(double rho,double pressy)
{
    double dpdrho;
    dpdrho=pressy*dpsi_deps_ply(log10(rho))*(1.0/rho);
    //dpdrho=pressy*deps_dpsi(log10(rho))*(1.0/rho);
    return dpdrho;
}

double drho_dp_analytical_ply(double pressy)
{
    double rho = rho_of_p_analytical_ply(pressy);
    double drhodp;
    drhodp=rho*deps_dpsi_ply(log10(pressy),rho)*(1.0/pressy);
    //drhodp=rho*1.0/dpsi_deps_ply(log10(rho))*(1.0/pressy);
    return drhodp;
}


double ddp_drhorho_ply(double rho,double pressy)
{
    double ddpddrho;
    ddpddrho=dp_drho_ply(rho,pressy)*dpsi_deps_ply(log10(rho))*(1.0/rho) + pressy*ddpsi_depseps_ply(log10(rho))*(1.0/(pow(rho,2)*log(10.0))) + pressy*deps_dpsi_ply(log10(pressy),rho)*(-1.0/pow(rho,2)) ;
    return ddpddrho;
}

double ddrho_dPP_analytical_ply(double pressy)
{
    double rho = rho_of_p_analytical_ply(pressy);
    double ddpddrho = ddp_drhorho_ply(rho, pressy);
    double dpdrho = dp_drho_ply(rho,pressy);
    double ddrhodPP;
    ddrhodPP = -ddpddrho/pow(dpdrho,3);
    return ddrhodPP;
}


double dpsi_deps_ply(double eps)
{
    double dpsideps;
    dpsideps = 2.0;
     return dpsideps;
}


double deps_dpsi_ply(double psi,double rho)
{
    //double eps = eps_of_psi(psi);
    double pressy = pow(10,psi); //exp(psi);
    //double rho = rho_of_press(pressy);
    double eps = log10(rho);

    double dpsideps = dpsi_deps_ply(eps);
    double depsdpsi=1.0/dpsideps;

    return depsdpsi;
}


double ddeps_dpsipsi_ply(double psi,double rho)
{
    //double eps = eps_of_psi(psi);
    double pressy = pow(10,psi); //exp(psi);
    //double rho = rho_of_press(pressy);
    double eps = log10(rho);
    double dpsideps = dpsi_deps_ply(eps);
    double ddpsidepseps = ddpsi_depseps_ply(eps);

    double ddepsdpsipsi = -ddpsidepseps/pow(dpsideps,3);

    return ddepsdpsipsi;
}


double ddpsi_depseps_ply(double eps)
{
    double ddpsidepseps = 0;

     return ddpsidepseps;
}

