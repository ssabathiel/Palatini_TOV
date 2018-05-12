#include "header/eos_analytical.h"
#include <vector>
#include <math.h>
#include <cmath>
#include "linalg_library/spline.h"
#include <iomanip>
#include "linalg_library/interpolation.h"
#include "header/glob_variables.h"



double a1 =6.22;
double a2 =6.121;
double a3 =0.005925;
double a4 =0.16326;
double a5 =6.48;
double a6 =11.4971;
double a7 =19.105;
double a8 =0.8938;
double a9 =6.54;
double a10 = 11.4950;
double a11 = -22.775;
double a12 = 1.5707;
double a13 = 4.3;
double a14 = 14.08;
double a15 = 27.80;
double a16 = -1.653;
double a17 = 1.50;
double a18 = 14.67;

double E = exp(1);

tk::spline pf;       //pf = p(rho)-function
tk::spline rhof;     //rhof = rho(p)-function



double Kf(double x)
{
    double k_local;
    k_local = 1.0/(exp(x) + 1);
    return k_local;
}


double psi_of_eps(double eps)
{
    double psi =  ((a1 + a2*eps + a3*pow(eps,3))/(1+a4*eps))*Kf(a5*(eps-a6)) + (a7+a8*eps)* Kf(a9*(a10-eps)) + (a11 + a12*eps)*Kf(a13*(a14-eps)) + (a15 + a16*eps)*Kf(a17*(a18-eps))    ;
    return psi;
}


double p_of_rho_analytical(double rho)
{
    double eps=log10(rho);
    double psi = psi_of_eps(eps);
    double pressy =  pow(10,psi); //exp(psi);          //--> should be 10^psi
    return pressy;
}





double rho_of_p_analytical(double p)
{
    double rho=-1;
    double p_copy=p;
    /*
    pf.set_points(rhos,press);
    rhof.set_points(press,rhos);
    rho = rhof(p_copy);
    */
    return alglib::spline1dcalc(rho_of_press_alg, p);

}



double Power(double base, double exp)
{
    return pow(base,exp);
}

double Log(double x)
{
    return log10(x);
}

double dp_drho(double rho,double pressy)
{
    double dpdrho;    
    dpdrho=pressy*dpsi_deps(log10(rho))*(1.0/rho);
    //dpdrho=pressy*deps_dpsi(log10(rho))*(1.0/rho);
    return dpdrho;
}

double drho_dp_analytical(double pressy)
{
    double rho = rhof(pressy);
    double drhodp;
    drhodp=rho*deps_dpsi(log10(pressy),rho)*(1.0/pressy);
    return drhodp;
}


double ddp_drhorho(double rho,double pressy)
{
    double ddpddrho;
    ddpddrho=dp_drho(rho,pressy)*dpsi_deps(log10(rho))*(1.0/rho) + pressy*ddpsi_depseps(log10(rho))*(1.0/pow(rho,2)) + pressy*deps_dpsi(log10(pressy),rho)*(-1.0/pow(rho,2)) ;
    return ddpddrho;
}

double ddrho_dPP_analytical(double pressy)
{
    double rho = rho_of_p_analytical(pressy);
    double ddpddrho = ddp_drhorho(rho, pressy);
    double dpdrho = dp_drho(rho,pressy);
    double ddrhodPP;
    ddrhodPP = -ddpddrho/pow(dpdrho,3);
    return ddrhodPP;
}


double dpsi_deps(double eps)
{
    double dpsideps;
    dpsideps = a8/(1 + Power(E,a9*(a10 - eps))) + a12/(1 + Power(E,a13*(a14 - eps))) + a16/(1 + Power(E,a17*(a18 - eps))) + (a13*Power(E,a13*(a14 - eps))*(a11 + a12*eps))/Power(1 + Power(E,a13*(a14 - eps)),2) +
            (a17*Power(E,a17*(a18 - eps))*(a15 + a16*eps))/Power(1 + Power(E,a17*(a18 - eps)),2) + (a9*Power(E,a9*(a10 - eps))*(a7 + a8*eps))/Power(1 + Power(E,a9*(a10 - eps)),2) +
            (a2 + 3*a3*Power(eps,2))/((1 + Power(E,a5*(-a6 + eps)))*(1 + a4*eps)) - (a4*(a1 + a2*eps + a3*Power(eps,3)))/((1 + Power(E,a5*(-a6 + eps)))*Power(1 + a4*eps,2)) -
            (a5*Power(E,a5*(-a6 + eps))*(a1 + a2*eps + a3*Power(eps,3)))/(Power(1 + Power(E,a5*(-a6 + eps)),2)*(1 + a4*eps));
     return dpsideps;
}


double deps_dpsi(double psi,double rho)
{
    //double eps = eps_of_psi(psi);
    double pressy = pow(10,psi); //exp(psi);
    //double rho = rho_of_press(pressy);
    double eps = log10(rho);

    double dpsideps = dpsi_deps(eps);
    double depsdpsi=1.0/dpsideps;

    return depsdpsi;
}


double ddeps_dpsipsi(double psi,double rho)
{
    //double eps = eps_of_psi(psi);
    double pressy = pow(10,psi); //exp(psi);
    //double rho = rho_of_press(pressy);
    double eps = log10(rho);
    double dpsideps = dpsi_deps(eps);
    double ddpsidepseps = ddpsi_depseps(eps);

    double ddepsdpsipsi = -ddpsidepseps/pow(dpsideps,3);

    return ddepsdpsipsi;
}


double ddpsi_depseps(double eps)
{
    double ddpsidepseps;
    ddpsidepseps = (2*a8*a9*Power(E,a9*(a10 - eps)))/Power(1 + Power(E,a9*(a10 - eps)),2) + (2*a12*a13*Power(E,a13*(a14 - eps)))/Power(1 + Power(E,a13*(a14 - eps)),2) +
            (2*Power(a13,2)*Power(E,2*a13*(a14 - eps))*(a11 + a12*eps))/Power(1 + Power(E,a13*(a14 - eps)),3) - (Power(a13,2)*Power(E,a13*(a14 - eps))*(a11 + a12*eps))/Power(1 + Power(E,a13*(a14 - eps)),2) +
            (6*a3*eps)/((1 + Power(E,a5*(-a6 + eps)))*(1 + a4*eps)) + (2*Power(a9,2)*Power(E,2*a9*(a10 - eps))*(a7 + a8*eps))/Power(1 + Power(E,a9*(a10 - eps)),3) -
            (Power(a9,2)*Power(E,a9*(a10 - eps))*(a7 + a8*eps))/Power(1 + Power(E,a9*(a10 - eps)),2) - (2*a4*(a2 + 3*a3*Power(eps,2)))/((1 + Power(E,a5*(-a6 + eps)))*Power(1 + a4*eps,2)) -
            (2*a5*Power(E,a5*(-a6 + eps))*(a2 + 3*a3*Power(eps,2)))/(Power(1 + Power(E,a5*(-a6 + eps)),2)*(1 + a4*eps)) +
            (2*Power(a4,2)*(a1 + a2*eps + a3*Power(eps,3)))/((1 + Power(E,a5*(-a6 + eps)))*Power(1 + a4*eps,3)) +
            (2*a4*a5*Power(E,a5*(-a6 + eps))*(a1 + a2*eps + a3*Power(eps,3)))/(Power(1 + Power(E,a5*(-a6 + eps)),2)*Power(1 + a4*eps,2)) +
            (2*Power(a5,2)*Power(E,2*a5*(-a6 + eps))*(a1 + a2*eps + a3*Power(eps,3)))/(Power(1 + Power(E,a5*(-a6 + eps)),3)*(1 + a4*eps)) -
            (Power(a5,2)*Power(E,a5*(-a6 + eps))*(a1 + a2*eps + a3*Power(eps,3)))/(Power(1 + Power(E,a5*(-a6 + eps)),2)*(1 + a4*eps));
     return ddpsidepseps;
}

