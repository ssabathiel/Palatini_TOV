#include "header/ap4_analytical.h"
#include "header/eos_analytical.h"
#include <vector>
#include <math.h>
#include <cmath>
#include "linalg_library/spline.h"
#include <iomanip>
#include "linalg_library/interpolation.h"
#include "header/glob_variables.h"

using namespace constants;






double psi_of_eps_ap4(double eps)
{
    double psi =  ((a1 + a2*eps + a3*pow(eps,3))/(1+a4*eps))*Kf(a5*(eps-a6)) + (a7+a8*eps)* Kf(a9*(a10-eps)) + (a11 + a12*eps)*Kf(a13*(a14-eps)) + (a15 + a16*eps)*Kf(a17*(a18-eps))    ;
    return psi;
}


double p_of_rho_analytical_ap4(double rho)
{
    double eps=log10(rho);
    double psi = psi_of_eps_ap4(eps);
    double pressy =  pow(10,psi); //exp(psi);          //--> should be 10^psi
    return pressy;
}





double rho_of_p_analytical_ap4(double p)
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





double dp_drho_ap4(double rho,double pressy)
{
    double dpdrho;
    dpdrho=pressy*dpsi_deps_ap4(log10(rho))*(1.0/rho);
    //dpdrho=pressy*deps_dpsi(log10(rho))*(1.0/rho);
    return dpdrho;
}

double drho_dp_analytical_ap4(double pressy)
{
    double rho = alglib::spline1dcalc(rho_of_press_alg, pressy);
    //double psi = log10(pressy);
    //rho = pow(10,(psi - 5.29355)/2.0);
    double drhodp;
    //std::cout << "rho= " << rho << std::endl;
    //std::cout << "p  = " << pressy << std::endl;
    drhodp=rho*deps_dpsi_ap4(log10(pressy),rho)*(1.0/pressy);
    //drhodp=rho*0.5*(1.0/pressy);
    //std::cout << "hey heyy" << std::endl;
    return drhodp;
}


double ddp_drhorho_ap4(double rho,double pressy)
{
    double ddpddrho;
    ddpddrho=dp_drho_ap4(rho,pressy)*dpsi_deps_ap4(log10(rho))*(1.0/rho) + pressy*ddpsi_depseps_ap4(log10(rho))*(1.0/(pow(rho,2)*log(10.0))) + pressy*deps_dpsi_ap4(log10(pressy),rho)*(-1.0/pow(rho,2)) ;
    return ddpddrho;
}

double ddrho_dPP_analytical_ap4(double pressy)
{
    double rho = rho_of_p_analytical_ap4(pressy);
    double ddpddrho = ddp_drhorho_ap4(rho, pressy);
    double dpdrho = dp_drho_ap4(rho,pressy);
    double ddrhodPP;
    ddrhodPP = -ddpddrho/pow(dpdrho,3);
    return ddrhodPP;
}


double dpsi_deps_ap4(double eps)
{
    double dpsideps;
    dpsideps = a8/(1 + Power(E,a9*(a10 - eps))) + a12/(1 + Power(E,a13*(a14 - eps))) + a16/(1 + Power(E,a17*(a18 - eps))) + (a13*Power(E,a13*(a14 - eps))*(a11 + a12*eps))/Power(1 + Power(E,a13*(a14 - eps)),2) +
            (a17*Power(E,a17*(a18 - eps))*(a15 + a16*eps))/Power(1 + Power(E,a17*(a18 - eps)),2) + (a9*Power(E,a9*(a10 - eps))*(a7 + a8*eps))/Power(1 + Power(E,a9*(a10 - eps)),2) +
            (a2 + 3*a3*Power(eps,2))/((1 + Power(E,a5*(-a6 + eps)))*(1 + a4*eps)) - (a4*(a1 + a2*eps + a3*Power(eps,3)))/((1 + Power(E,a5*(-a6 + eps)))*Power(1 + a4*eps,2)) -
            (a5*Power(E,a5*(-a6 + eps))*(a1 + a2*eps + a3*Power(eps,3)))/(Power(1 + Power(E,a5*(-a6 + eps)),2)*(1 + a4*eps));
/*
    dpsideps = 0.8938/(1 + Power(E,6.54*(11.495 - eps))) + 1.5707/(1 + Power(E,4.3*(14.08 - eps))) -
            1.653/(1 + Power(E,1.5*(14.67 - eps))) +
            (1.5*Power(E,1.5*(14.67 - eps))*(27.8 - 1.653*eps))/
             Power(1 + Power(E,1.5*(14.67 - eps)),2) +
            (6.54*Power(E,6.54*(11.495 - eps))*(19.105 + 0.8938*eps))/
             Power(1 + Power(E,6.54*(11.495 - eps)),2) +
            (4.3*Power(E,4.3*(14.08 - eps))*(-22.775 + 1.5707*eps))/
             Power(1 + Power(E,4.3*(14.08 - eps)),2) +
            (6.121 + 0.017775*Power(eps,2))/((1 + Power(E,6.48*(-11.4971 + eps)))*(1 + 0.16326*eps)) -
            (0.16326*(6.22 + 6.121*eps + 0.005925*Power(eps,3)))/
             ((1 + Power(E,6.48*(-11.4971 + eps)))*Power(1 + 0.16326*eps,2)) -
            (6.48*Power(E,6.48*(-11.4971 + eps))*(6.22 + 6.121*eps + 0.005925*Power(eps,3)))/
             (Power(1 + Power(E,6.48*(-11.4971 + eps)),2)*(1 + 0.16326*eps));
             */
    return dpsideps;
}


double deps_dpsi_ap4(double psi,double rho)
{
    //double eps = eps_of_psi(psi);
    double pressy = pow(10,psi); //exp(psi);
    //double rho = rho_of_press(pressy);
    double eps = log10(rho);

    double dpsideps = dpsi_deps_ap4(eps);
    double depsdpsi=1.0/dpsideps;

    return depsdpsi;
}


double ddeps_dpsipsi_ap4(double psi,double rho)
{
    //double eps = eps_of_psi(psi);
    double pressy = pow(10,psi); //exp(psi);
    //double rho = rho_of_press(pressy);
    double eps = log10(rho);
    double dpsideps = dpsi_deps_ap4(eps);
    double ddpsidepseps = ddpsi_depseps_ap4(eps);

    double ddepsdpsipsi = -ddpsidepseps/pow(dpsideps,3);

    return ddepsdpsipsi;
}


double ddpsi_depseps_ap4(double eps)
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
    ddpsidepseps = (11.690904*Power(E,6.54*(11.495 - eps)))/Power(1.0 + Power(E,6.54*(11.495 - eps)),2) +
       (13.50802*Power(E,4.3*(14.08 - eps)))/Power(1.0 + Power(E,4.3*(14.08 - eps)),2) -
       (4.959*Power(E,1.5*(14.67 - eps)))/Power(1 + Power(E,1.5*(14.67 - eps)),2) +
       (4.5*Power(E,3.*(14.67 - eps))*(27.8 - 1.653*eps))/Power(1 + Power(E,1.5*(14.67 - eps)),3) -
       (2.25*Power(E,1.5*(14.67 - eps))*(27.8 - 1.653*eps))/
        Power(1.0 + Power(E,1.5*(14.67 - eps)),2) +
       (85.5432*Power(E,13.08*(11.495 - eps))*(19.105 + 0.8938*eps))/
        Power(1.0 + Power(E,6.54*(11.495 - eps)),3) -
       (42.7716*Power(E,6.54*(11.495 - eps))*(19.105 + 0.8938*eps))/
        Power(1.0 + Power(E,6.54*(11.495 - eps)),2) +
       (0.03555*eps)/((1 + Power(E,6.48*(-11.4971 + eps)))*(1 + 0.16326*eps)) +
       (36.98*Power(E,8.6*(14.08 - eps))*(-22.775 + 1.5707*eps))/
        Power(1.0 + Power(E,4.3*(14.08 - eps)),3) -
       (18.49*Power(E,4.3*(14.08 - eps))*(-22.775 + 1.5707*eps))/
        Power(1.0 + Power(E,4.3*(14.08 - eps)),2) -
       (0.32652*(6.121 + 0.017775*Power(eps,2)))/
        ((1.0 + Power(E,6.48*(-11.4971 + eps)))*Power(1.0 + 0.16326*eps,2)) -
       (12.96*Power(E,6.48*(-11.4971 + eps))*(6.121 + 0.017775*Power(eps,2)))/
        (Power(1.0 + Power(E,6.48*(-11.4971 + eps)),2)*(1.0 + 0.16326*eps)) +
       (0.05330765519999999*(6.22 + 6.121*eps + 0.005925*Power(eps,3)))/
        ((1.0 + Power(E,6.48*(-11.4971 + eps)))*Power(1.0 + 0.16326*eps,3)) +
       (2.1158495999999998*Power(E,6.48*(-11.4971 + eps))*
          (6.22 + 6.121*eps + 0.005925*Power(eps,3)))/
        (Power(1.0 + Power(E,6.48*(-11.4971 + eps)),2)*Power(1.0 + 0.16326*eps,2)) +
       (83.98080000000002*Power(E,12.96*(-11.4971 + eps))*
          (6.22 + 6.121*eps + 0.005925*Power(eps,3)))/
        (Power(1 + Power(E,6.48*(-11.4971 + eps)),3)*(1.0 + 0.16326*eps)) -
       (41.99040000000001*Power(E,6.48*(-11.4971 + eps))*
          (6.22 + 6.121*eps + 0.005925*Power(eps,3)))/
        (Power(1.0 + Power(E,6.48*(-11.4971 + eps)),2)*(1.0 + 0.16326*eps));

    ddpsidepseps = (2*a8*a9*Power(E,a9*(a10 - eps)))/Power(1 + Power(E,a9*(a10 - eps)),2) +
       (2*a12*a13*Power(E,a13*(a14 - eps)))/Power(1 + Power(E,a13*(a14 - eps)),2) +
       (2*a16*a17*Power(E,a17*(a18 - eps)))/Power(1 + Power(E,a17*(a18 - eps)),2) +
       (2*Power(a13,2)*Power(E,2*a13*(a14 - eps))*(a11 + a12*eps))/
        Power(1 + Power(E,a13*(a14 - eps)),3) -
       (Power(a13,2)*Power(E,a13*(a14 - eps))*(a11 + a12*eps))/
        Power(1 + Power(E,a13*(a14 - eps)),2) +
       (2*Power(a17,2)*Power(E,2*a17*(a18 - eps))*(a15 + a16*eps))/
        Power(1 + Power(E,a17*(a18 - eps)),3) -
       (Power(a17,2)*Power(E,a17*(a18 - eps))*(a15 + a16*eps))/
        Power(1 + Power(E,a17*(a18 - eps)),2) +
       (6*a3*eps)/((1 + Power(E,a5*(-a6 + eps)))*(1 + a4*eps)) +
       (2*Power(a9,2)*Power(E,2*a9*(a10 - eps))*(a7 + a8*eps))/
        Power(1 + Power(E,a9*(a10 - eps)),3) -
       (Power(a9,2)*Power(E,a9*(a10 - eps))*(a7 + a8*eps))/Power(1 + Power(E,a9*(a10 - eps)),2) -
       (2*a4*(a2 + 3*a3*Power(eps,2)))/((1 + Power(E,a5*(-a6 + eps)))*Power(1 + a4*eps,2)) -
       (2*a5*Power(E,a5*(-a6 + eps))*(a2 + 3*a3*Power(eps,2)))/
        (Power(1 + Power(E,a5*(-a6 + eps)),2)*(1 + a4*eps)) +
       (2*Power(a4,2)*(a1 + a2*eps + a3*Power(eps,3)))/
        ((1 + Power(E,a5*(-a6 + eps)))*Power(1 + a4*eps,3)) +
       (2*a4*a5*Power(E,a5*(-a6 + eps))*(a1 + a2*eps + a3*Power(eps,3)))/
        (Power(1 + Power(E,a5*(-a6 + eps)),2)*Power(1 + a4*eps,2)) +
       (2*Power(a5,2)*Power(E,2*a5*(-a6 + eps))*(a1 + a2*eps + a3*Power(eps,3)))/
        (Power(1 + Power(E,a5*(-a6 + eps)),3)*(1 + a4*eps)) -
       (Power(a5,2)*Power(E,a5*(-a6 + eps))*(a1 + a2*eps + a3*Power(eps,3)))/
        (Power(1 + Power(E,a5*(-a6 + eps)),2)*(1 + a4*eps));



     return ddpsidepseps;
}


