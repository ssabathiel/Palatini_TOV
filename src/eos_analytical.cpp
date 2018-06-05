#include "header/eos_analytical.h"
#include <vector>
#include <math.h>
#include <cmath>
#include "linalg_library/spline.h"
#include <iomanip>
#include "linalg_library/interpolation.h"
#include "header/glob_variables.h"

double c1 = 10.6557;
double c2 = 3.7863;
double c3 = 0.8124;
double c4 = 0.6823;
double c5 = 3.5279;
double c6 = 11.8100;
double c7 = 12.0584;
double c8 = 1.4663;
double c9 = 3.4952;
double c10 = 11.8007;
double c11 = 14.4114;
double c12 = 14.4081;


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
    double psi;
    if(analytical_EOS==1 || fps_EOS==1)
    {
        psi =  ((a1 + a2*eps + a3*pow(eps,3))/(1+a4*eps))*Kf(a5*(eps-a6)) + (a7+a8*eps)* Kf(a9*(a10-eps)) + (a11 + a12*eps)*Kf(a13*(a14-eps)) + (a15 + a16*eps)*Kf(a17*(a18-eps))    ;
    }
    if(ap4_EOS==1)
    {
        psi =((10.6557 + 3.7863*Power(-0.8124 + eps,0.6823))/(1 + Power(E,3.5279*(-11.81 + eps))) +
              (12.0584 + 1.4663*eps)/(1 + Power(E,3.4952*(11.8007 - eps))))/
            (1 + Power(E,4.329*(-14.4114 + eps))) +
           ((9.1131 - 0.475*eps)/(1 + Power(E,3.4614*(14.88 - eps))) +
              (21.3141 + 0.1023*eps + 0.0495*Power(eps,2))/(1 + Power(E,4.9401*(10.2957 - eps))))/
            (1 + Power(E,4.3622*(14.4081 - eps)));
    }
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
    double rho = alglib::spline1dcalc(rho_of_press_alg, pressy);
    //double psi = log10(pressy);
    //rho = pow(10,(psi - 5.29355)/2.0);
    double drhodp;
    //std::cout << "rho= " << rho << std::endl;
    //std::cout << "p  = " << pressy << std::endl;
    drhodp=rho*deps_dpsi(log10(pressy),rho)*(1.0/pressy);
    //drhodp=rho*0.5*(1.0/pressy);
    //std::cout << "hey heyy" << std::endl;
    return drhodp;
}


double ddp_drhorho(double rho,double pressy)
{
    double ddpddrho;
    ddpddrho=dp_drho(rho,pressy)*dpsi_deps(log10(rho))*(1.0/rho) + pressy*ddpsi_depseps(log10(rho))*(1.0/(pow(rho,2)*log(10.0))) + pressy*deps_dpsi(log10(pressy),rho)*(-1.0/pow(rho,2)) ;
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
    if(analytical_EOS==1 || fps_EOS==1)
    {
        dpsideps = a8/(1 + Power(E,a9*(a10 - eps))) + a12/(1 + Power(E,a13*(a14 - eps))) + a16/(1 + Power(E,a17*(a18 - eps))) + (a13*Power(E,a13*(a14 - eps))*(a11 + a12*eps))/Power(1 + Power(E,a13*(a14 - eps)),2) +
                (a17*Power(E,a17*(a18 - eps))*(a15 + a16*eps))/Power(1 + Power(E,a17*(a18 - eps)),2) + (a9*Power(E,a9*(a10 - eps))*(a7 + a8*eps))/Power(1 + Power(E,a9*(a10 - eps)),2) +
                (a2 + 3*a3*Power(eps,2))/((1 + Power(E,a5*(-a6 + eps)))*(1 + a4*eps)) - (a4*(a1 + a2*eps + a3*Power(eps,3)))/((1 + Power(E,a5*(-a6 + eps)))*Power(1 + a4*eps,2)) -
                (a5*Power(E,a5*(-a6 + eps))*(a1 + a2*eps + a3*Power(eps,3)))/(Power(1 + Power(E,a5*(-a6 + eps)),2)*(1 + a4*eps));

    }
    if(ap4_EOS==1)
    {

        dpsideps =(1.4663/(1 + Power(E,3.4952*(11.8007 - eps))) -
                   (3.5279*Power(E,3.5279*(-11.81 + eps))*(10.6557 + 3.7863*Power(-0.8124 + eps,0.6823)))/
                    Power(1 + Power(E,3.5279*(-11.81 + eps)),2) +
                   2.58339249/((1 + Power(E,3.5279*(-11.81 + eps)))*Power(-0.8124 + eps,0.3177)) +
                   (3.4952*Power(E,3.4952*(11.8007 - eps))*(12.0584 + 1.4663*eps))/
                    Power(1 + Power(E,3.4952*(11.8007 - eps)),2))/(1 + Power(E,4.329*(-14.4114 + eps))) -
                (4.329*Power(E,4.329*(-14.4114 + eps))*
                   ((10.6557 + 3.7863*Power(-0.8124 + eps,0.6823))/(1 + Power(E,3.5279*(-11.81 + eps))) +
                     (12.0584 + 1.4663*eps)/(1 + Power(E,3.4952*(11.8007 - eps)))))/
                 Power(1 + Power(E,4.329*(-14.4114 + eps)),2) +
                (-0.475/(1 + Power(E,3.4614*(14.88 - eps))) +
                   (3.4614*Power(E,3.4614*(14.88 - eps))*(9.1131 - 0.475*eps))/
                    Power(1 + Power(E,3.4614*(14.88 - eps)),2) +
                   (0.1023 + 0.099*eps)/(1 + Power(E,4.9401*(10.2957 - eps))) +
                   (4.9401*Power(E,4.9401*(10.2957 - eps))*(21.3141 + 0.1023*eps + 0.0495*Power(eps,2)))/
                    Power(1 + Power(E,4.9401*(10.2957 - eps)),2))/(1 + Power(E,4.3622*(14.4081 - eps))) +
                (4.3622*Power(E,4.3622*(14.4081 - eps))*
                   ((9.1131 - 0.475*eps)/(1 + Power(E,3.4614*(14.88 - eps))) +
                     (21.3141 + 0.1023*eps + 0.0495*Power(eps,2))/(1 + Power(E,4.9401*(10.2957 - eps)))))/
                 Power(1 + Power(E,4.3622*(14.4081 - eps)),2);
    }

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


    if(analytical_EOS==1 || fps_EOS==1)
    {
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

    }
    if(ap4_EOS==1)
    {
        ddpsidepseps = ((10.25002352*Power(E,3.4952*(11.8007 - eps)))/Power(1 + Power(E,3.4952*(11.8007 - eps)),2) +
                        (24.892156819999997*Power(E,7.0558*(-11.81 + eps))*
                           (10.6557 + 3.7863*Power(-0.8124 + eps,0.6823)))/
                         Power(1 + Power(E,3.5279*(-11.81 + eps)),3) -
                        (12.446078409999998*Power(E,3.5279*(-11.81 + eps))*
                           (10.6557 + 3.7863*Power(-0.8124 + eps,0.6823)))/
                         Power(1 + Power(E,3.5279*(-11.81 + eps)),2) -
                        0.820743794073/
                         ((1 + Power(E,3.5279*(-11.81 + eps)))*Power(-0.8124 + eps,1.3176999999999999)) -
                        (18.227900730942*Power(E,3.5279*(-11.81 + eps)))/
                         (Power(1 + Power(E,3.5279*(-11.81 + eps)),2)*Power(-0.8124 + eps,0.3177)) +
                        (24.43284608*Power(E,6.9904*(11.8007 - eps))*(12.0584 + 1.4663*eps))/
                         Power(1 + Power(E,3.4952*(11.8007 - eps)),3) -
                        (12.21642304*Power(E,3.4952*(11.8007 - eps))*(12.0584 + 1.4663*eps))/
                         Power(1 + Power(E,3.4952*(11.8007 - eps)),2))/(1 + Power(E,4.329*(-14.4114 + eps))) -
                     (8.658*Power(E,4.329*(-14.4114 + eps))*
                        (1.4663/(1 + Power(E,3.4952*(11.8007 - eps))) -
                          (3.5279*Power(E,3.5279*(-11.81 + eps))*(10.6557 + 3.7863*Power(-0.8124 + eps,0.6823)))/
                           Power(1 + Power(E,3.5279*(-11.81 + eps)),2) +
                          2.58339249/((1 + Power(E,3.5279*(-11.81 + eps)))*Power(-0.8124 + eps,0.3177)) +
                          (3.4952*Power(E,3.4952*(11.8007 - eps))*(12.0584 + 1.4663*eps))/
                           Power(1 + Power(E,3.4952*(11.8007 - eps)),2)))/
                      Power(1 + Power(E,4.329*(-14.4114 + eps)),2) +
                     (37.480481999999995*Power(E,8.658*(-14.4114 + eps))*
                        ((10.6557 + 3.7863*Power(-0.8124 + eps,0.6823))/(1 + Power(E,3.5279*(-11.81 + eps))) +
                          (12.0584 + 1.4663*eps)/(1 + Power(E,3.4952*(11.8007 - eps)))))/
                      Power(1 + Power(E,4.329*(-14.4114 + eps)),3) -
                     (18.740240999999997*Power(E,4.329*(-14.4114 + eps))*
                        ((10.6557 + 3.7863*Power(-0.8124 + eps,0.6823))/(1 + Power(E,3.5279*(-11.81 + eps))) +
                          (12.0584 + 1.4663*eps)/(1 + Power(E,3.4952*(11.8007 - eps)))))/
                      Power(1 + Power(E,4.329*(-14.4114 + eps)),2) +
                     (0.099/(1 + Power(E,4.9401*(10.2957 - eps))) -
                        (3.2883299999999998*Power(E,3.4614*(14.88 - eps)))/
                         Power(1 + Power(E,3.4614*(14.88 - eps)),2) +
                        (23.962579919999996*Power(E,6.9228*(14.88 - eps))*(9.1131 - 0.475*eps))/
                         Power(1 + Power(E,3.4614*(14.88 - eps)),3) -
                        (11.981289959999998*Power(E,3.4614*(14.88 - eps))*(9.1131 - 0.475*eps))/
                         Power(1 + Power(E,3.4614*(14.88 - eps)),2) +
                        (9.8802*Power(E,4.9401*(10.2957 - eps))*(0.1023 + 0.099*eps))/
                         Power(1 + Power(E,4.9401*(10.2957 - eps)),2) +
                        (48.80917602*Power(E,9.8802*(10.2957 - eps))*
                           (21.3141 + 0.1023*eps + 0.0495*Power(eps,2)))/
                         Power(1 + Power(E,4.9401*(10.2957 - eps)),3) -
                        (24.40458801*Power(E,4.9401*(10.2957 - eps))*
                           (21.3141 + 0.1023*eps + 0.0495*Power(eps,2)))/
                         Power(1 + Power(E,4.9401*(10.2957 - eps)),2))/(1 + Power(E,4.3622*(14.4081 - eps))) +
                     (8.7244*Power(E,4.3622*(14.4081 - eps))*
                        (-0.475/(1 + Power(E,3.4614*(14.88 - eps))) +
                          (3.4614*Power(E,3.4614*(14.88 - eps))*(9.1131 - 0.475*eps))/
                           Power(1 + Power(E,3.4614*(14.88 - eps)),2) +
                          (0.1023 + 0.099*eps)/(1 + Power(E,4.9401*(10.2957 - eps))) +
                          (4.9401*Power(E,4.9401*(10.2957 - eps))*(21.3141 + 0.1023*eps + 0.0495*Power(eps,2)))/
                           Power(1 + Power(E,4.9401*(10.2957 - eps)),2)))/
                      Power(1 + Power(E,4.3622*(14.4081 - eps)),2) +
                     (38.057577679999994*Power(E,8.7244*(14.4081 - eps))*
                        ((9.1131 - 0.475*eps)/(1 + Power(E,3.4614*(14.88 - eps))) +
                          (21.3141 + 0.1023*eps + 0.0495*Power(eps,2))/(1 + Power(E,4.9401*(10.2957 - eps)))))/
                      Power(1 + Power(E,4.3622*(14.4081 - eps)),3) -
                     (19.028788839999997*Power(E,4.3622*(14.4081 - eps))*
                        ((9.1131 - 0.475*eps)/(1 + Power(E,3.4614*(14.88 - eps))) +
                          (21.3141 + 0.1023*eps + 0.0495*Power(eps,2))/(1 + Power(E,4.9401*(10.2957 - eps)))))/
                      Power(1 + Power(E,4.3622*(14.4081 - eps)),2);
    }

     return ddpsidepseps;
}

