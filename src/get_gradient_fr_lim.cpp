
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <cstdlib>  // for atof string to double conversion
#include "linalg_library/spline.h"
#include <iomanip>
#include "linalg_library/interpolation.h"
#include "header/plotradmassrho.h"
#include "header/eos_filetoarray.h"
#include "header/eos_functions.h"
#include "header/glob_variables.h"
#include "header/eos_analytical.h"
#include "header/meta_functions.h"
#include "header/get_gradient_fr.h"

using namespace constants;
using namespace std;

//double p_m = +1;

pair<double, double> get_gradients_fR_lim(double m, double press, double r)
{


    if(press<min_press){press=min_press;}
    else{press = press;}

    double rho = rho_of_p[num_an](press);


    ////////////////////////////////////////
    //PRECALCULATE TERMS NECESSARY FOR THE GRADIENTS
    /////////////////////////////////////////



        ////////////////////
        //dPdr PREPARATION
        //////////////////////
        double drhodP = drho_dp[num_an](press);
        double ddrhodPP = ddrho_dPP[num_an](press);

        double T=(-rho+3*press/pow(clight,2));
        double R=-kappa_2*T;
        double f=R+a*pow(R,2)/Rp;

       double f_tilde = R+a*pow(R,2)/Rp;
       double f_tildeR = 1 + 2*a*R/Rp;
       double fR=1+2*a*R/Rp;
       double b=1.0;

       double dRdP = -kappa_2*(   -drhodP + (3/pow(clight,2) )    );
       double dfdP = fR*dRdP;
       double dR_2dP = 2*R*dRdP;
       double ddRdPP = kappa_2*ddrhodPP;
       double ddR_2dPP = 2*R*dRdP; //18* kappa_2;


       double dfRdP = (2*a/Rp)*dRdP;					// multiplied Rp artificially!!!!!!!!!!!!!!!!
       double ddfRdPP = (2*a/Rp)*(-kappa_2*(-ddrhodPP));
       double dfR_2dP = 2*fR*dfRdP; //4*a*dRdP/Rp + 4*pow(a,2)*dR_2dP/pow(Rp,2);
       double ddfR_2dPP = 2*(dfRdP*dfRdP + fR*ddfRdPP); //4*a*a/pow(Rp,2) * kappa_2*kappa_2*18;

       double sigma_1=fR;
       double sigma_2=fR;

       double S=fR;
       double Omega=fR;

       double Ar=1-2*m*ggrav/(r*pow(clight,2));

       double det_sigma=pow(fR,4);           //1/math.sqrt(sigma1 * sigma2 * sigma2 * sigma2)
       double sqrt_det_sigma=pow(det_sigma,0.5);

       double tau_rr=(1/sqrt_det_sigma)*(f*0.5+kappa_2*press/pow(clight,2));    //T_rr=pressure for a perfect fluid
       double tau_tt=(1/sqrt_det_sigma)*(f*0.5+kappa_2*(-rho));     //press cancel each other

       double dSigma1dP = dfRdP;
       double dSigma2dP = dfRdP;

       double dOmegadP = dfRdP;
       double dSdP = dfRdP;

        //dp_dr_limt
        double P_r_0 = ggrav*((rho+press/pow(clight,2))/(r*(r-2*ggrav*m/pow(clight,2))))*(m-(tau_rr+tau_tt)*pow(r,3)*0.25*pow(clight,2)/ggrav);
        double alpha_r = (rho+press/pow(clight,2))*dfRdP/fR;
        double beta_r = 2*r*(dfRdP/fR)*(1-(rho+press/pow(clight,2))*(3*dfRdP)/(4*fR) );
        double dPdr = - (P_r_0*2)/( (1-alpha_r)*( 1+ pow(1-beta_r*P_r_0,0.5) )     );



        double dfRdr = dfRdP*dPdr;
        double ddSdPP = ddfRdPP;
        double ddOmegadPP = ddfRdPP;
        double dOmegaPDOmegadP = dOmegadP*((-1)/pow(Omega,2))*dOmegadP + ddOmegadPP/Omega;
        double dSPDSdP = dSdP*((-1)/pow(S,2))*dSdP + ddSdPP/S;


        double dalpha_rdr = 0.5*(drho_dp[num_an](press)+ 1/pow(clight,2) )*(dOmegadP/Omega + dSdP/S) + (rho + press/pow(clight,2))*0.5*(dOmegaPDOmegadP + dSPDSdP);
        dalpha_rdr = dalpha_rdr*dPdr;
        double dbeta_rdr = 2*dOmegadP/Omega + 2*r*dOmegaPDOmegadP;
        dbeta_rdr = dbeta_rdr + 2*dOmegadP/Omega*(rho+press/pow(clight,2))*0.5*(0.5*dOmegadP/Omega-dSPDSdP/S)
                              + 2*r*dOmegaPDOmegadP*(rho+press/pow(clight,2))*0.5*(0.5*dOmegadP/Omega-dSPDSdP/S)
                              + 2*r*dOmegadP/Omega*0.5*(0.5*dOmegadP/Omega-dSPDSdP/S)
                              + 2*r*dOmegadP/Omega*(rho+press/pow(clight,2))*0.5*(0.5*dOmegaPDOmegadP - dSPDSdP);
        dbeta_rdr = dbeta_rdr*dPdr;



        double ss=1;
        double aa= 1+                        0.5*ss* ( (beta_r*P_r_0)/(pow(1-beta_r*P_r_0,0.5)*(1+ss*pow(1-beta_r*P_r_0,0.5)))  );
        double bb=(dalpha_rdr/(1-alpha_r)) + 0.5*ss* ( (dbeta_rdr*P_r_0)/(pow(1-beta_r*P_r_0,0.5)*(1+ss*pow(1-beta_r*P_r_0,0.5)))  );   //same 2nd term as in a

        double big_phi = (tau_rr + (Omega/S)*tau_tt);

        double drhodp = drho_dp[num_an](press);//(1.0/gamma_now)*pow(press/K_now, (1.0/gamma_now)-1 )*(1.0/K_now);
        double d1Dsqrt_det_sigmadP = -0.5*pow(det_sigma,-3./2.)*(pow(sigma_2,3)*dSigma1dP + 3*sigma_1*pow(sigma_2,2)*dSigma2dP);
        double dtau_rrdP = d1Dsqrt_det_sigmadP*(f*0.5+kappa_2*press)+ 1/sqrt_det_sigma*(dfdP*0.5+kappa_2);
        double dtau_ttdP = d1Dsqrt_det_sigmadP*(f*0.5-kappa_2*rho)+ 1/sqrt_det_sigma*(dfdP*0.5);

        double dPhidP = dtau_rrdP + dtau_ttdP*(Omega/S) + tau_tt*(dOmegadP/S - Omega/pow(S,2)*dSdP);
        double cc = ( (1+drhodp)/(rho+press/pow(clight,2)) - (dPhidP*pow(r,3)*0.25)/(m-big_phi * pow(r,3)*0.25)  )*dPdr -
        (  (2*(r-m*ggrav/pow(clight,2)))/(r*(r-2*m*ggrav/pow(clight,2))) + (3*big_phi*r*r*0.25 )/(m-big_phi*pow(r,3)*0.25)   );
        double dd = (2/(r-2*m*ggrav/pow(clight,2)) + 1/(m-big_phi*r*r*r*0.25));



        double xx = ((dfRdr/fR) + (2.0/r) )*(1.0/r);
        double yy = (3*tau_rr - tau_tt)/(2.0);
        double zz = (dfRdr/fR)*( (2*r-3*ggrav*m/pow(clight,2))/(r*(r-2*ggrav*m/pow(clight,2))   )  - 0.75*(dfRdr)/(fR)   );

        double ww = Ar*(1/fR)*kappa_2*2*a*(1/Rp)*(3/pow(clight,2)-drhodP)*dd*aa*dPdr;




        double P_rr_0 = dPdr*(r*(r-2*m*ggrav/pow(clight,2)))*(m-(tau_rr+Omega*(1/S)*tau_tt)*pow(r,3)*0.25)
                               +(rho+press/pow(clight,2))*(2*r-2*m*ggrav/pow(clight,2))*(m-(tau_rr+Omega*(1/S)*tau_tt)*pow(r,3)*0.25)
                               +(rho+press/pow(clight,2))*(r*(r-2*m*ggrav/pow(clight,2)))*(  (3*r*r/4.)*m - dtau_rrdP*dPdr*pow(r,3)/4 - tau_rr*3+pow(r,2)/4
                                                            - dOmegadP*dPdr*(1/S)*tau_tt*pow(r,3)*0.25
                                                            - Omega*(-1)*pow(S,-2)*dSdP*dPdr*tau_tt*pow(r,3)*0.25
                                                            - Omega*(1/S)*dtau_ttdP*dPdr*pow(r,3)*0.25
                                                            - Omega*(1/S)*tau_tt*3*pow(r,2)*0.25);

        double P_rr= -P_rr_0*(1/(1-alpha_r))*2*1/(1+pow(1-beta_r*P_r_0,0.5))
                             -P_r_0*(-1)*pow(1-alpha_r,-2)*dalpha_rdr*2*1/(1+pow(1-beta_r*P_r_0,0.5))
                             -P_rr_0*(1/(1-alpha_r))*2*(-1)*pow(1+pow(1-beta_r*P_r_0,0.5),-2)
                                                           *0.5*pow(1-beta_r*P_r_0,-0.5)
                                                           *(-1)*(dbeta_rdr*P_r_0 + beta_r*P_rr_0);







        double ddrhodrr = ddrhodPP*pow(dPdr,2) + drhodP*P_rr;                                                                                   //MODIFIED 2nd derivative
        double vv = -Ar*( (1/fR)*kappa_2*(2*a/Rp)*(-ddrhodrr  + 3*(cc*aa+bb)*dPdr/pow(clight,2) )   );

        double dm_dr_lim = (yy+vv+Ar*zz)/(xx+ww)*( (pow(clight,2))/(ggrav) );

        if(r>1*pow(10,-6))     // does not matter how I cange this value: always at this value dm_dr gets negative
        {
            dPdr = dPdr;
            dm_dr_lim = dm_dr_lim;
        //dp_dr = (m+4.0*pi*r*r*r*press/pow(clight,2))/(r*(r-2.0*ggrav*m/pow(clight,2)));
        }
        else
        {
            dPdr = (4.0*pi*r*press/pow(clight,2));
            dPdr = -ggrav*( rho+press/pow(clight,2) )*dPdr;
            dm_dr_lim = 4*pi*rho*pow(r,2);
        }

    double dp_dr = dPdr;
    double dm_dr = dm_dr_lim;

    return make_pair(dm_dr, dp_dr);
}






