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

#include "header/get_gradient_frq.h"
#include "header/get_gradient_fr.h"
#include "header/get_gradient_gr.h"
#include "header/get_gradient_fr_lim.h"

#include "header/eos_analytical.h"
#include "header/meta_functions.h"

using namespace constants;
using namespace std;

/////////////////////////////
//TOV EQUATION: dm/dr, dp/dr
/////////////////////////////



get_gradients_functions get_gradients2[]=
    {
        get_gradients_GR,
        get_gradients_fR,
        get_gradients_fRQ,
        get_gradients_fR_lim
    };



pair<double, double> get_gradients_fRQ(double m, double press, double r)
{
    if(press<min_press){press=min_press;}
    else{press = press;}

    //Get dm/dr=4*pi*rho², rho need to be get from p first
    double rho = rho_of_p[num_an](press);
    double drhodP = drho_dp[num_an](press);

    double ddrhodPP = ddrho_dPP[num_an](press);

    ////////////////////////////////////////
    //PRECALCULATE TERMS NECESSARY FOR THE GRADIENTS
    /////////////////////////////////////////


    ////////////////////
    //dPdr PREPARATION
    //////////////////////
    double T=(-rho+3*press/pow(clight,2));
    double R=-kappa_2*T;
    double f_tilde = R+a*pow(R,2)/Rp;
    double f_tildeR = 1 + 2*a*R/Rp;
    double b=1.0;
    if(Rp/(32.0) * pow(     -  pow( (R)/(Rp) + f_tildeR,2 ) -(4*kappa_2*(rho+press/pow(clight,2) ) )/(Rp), 1  )<0){b=-1.0;}
    double Q=(3*pow(Rq,2)/8)*(1-(2*kappa_2*(rho+press/pow(clight,2))/Rq)+(2*pow(kappa_2,2)*pow((rho-3*press/pow(clight,2)),2)/(3*pow(Rq,2)))-pow(1-(4*kappa_2*(rho+press/pow(clight,2)))/Rq,0.5));     //make Energy-dependent??

    Q = 2*Rp* (-(kappa_2*press/pow(clight,2) + f_tilde*0.5 + (Rp/(8.0*b))*pow(f_tildeR,2) )
    + Rp/(32.0*b) * pow(    3*((R*b)/(Rp) + f_tildeR) - pow( pow( (R*b)/(Rp) + f_tildeR,2 ) -(4*b*kappa_2*(rho+press/pow(clight,2) ) )/(Rp), 0.5  ),2) ); //From Helios--> Olmo --> p/m-paper

    double f=R+a*pow(R,2)/Rp + Q/Rq;
    double fR=1+2*a*R/Rp;
    double fRR = 2*a/Rp;
    double fQ=1/Rq;

    double dRdT = -kappa_2;
    double dTdP = -drhodP + 3.0/pow(clight,2);
    double dRdP = dRdT*dTdP;
    double dR_2dP = 2*R*dRdP;
    double ddRdPP = kappa_2*ddrhodPP;
    double ddR_2dPP = 2*R*dRdP;

    double dQdP = (3*pow(Rq,2)/8)*(-2*kappa_2*(drhodP + 1/pow(clight,2))/Rq + (2*pow(kappa_2,2)*2*(rho-3*press/pow(clight,2))*(drhodP-3/pow(clight,2)))/(3*pow(Rq,2))-0.5*pow(1-4*kappa_2*(rho+press/pow(clight,2))/Rq,-0.5)*(-4*kappa_2*(drhodP + 1/pow(clight,2))/Rq));
    //ddQdPP = (3*pow(Rq,2)/8) * ( 36*pow(kappa_2,2)/(3*Rq*Rq) +0.25*pow(1-4*kappa_2*(rho+press/pow(clight,2))/Rq,-3./2.)*pow(-4*kappa_2/Rq,2) );
    double ddQdPP = 3*Rq*Rq*(1./8.)* (  - (2*kappa_2*(ddrhodPP))/(Rq)
                                + 2*kappa_2*2*( (drhodP - 3/pow(clight,2))*(drhodP - 3/pow(clight,2))* (rho-3*press/pow(clight,2))*ddrhodPP   )/(3*Rq*Rq)
                                - 0.5*(  -0.5*pow(1-(4*kappa_2*(rho+press/pow(clight,2) ))/(Rq), -1.5 ) * (-4*kappa_2*(drhodP + 1/pow(clight,2)/Rq))* (-4*kappa_2*(drhodP + 1/pow(clight,2)/Rq))
                                          + (pow(1-(4*kappa_2*(rho+press/pow(clight,2) ))/(Rq), -0.5 ))*(4*kappa_2*(ddrhodPP)/Rq)   )
                            );

    double dfdP = fR*dRdP;
    double ddfdPP = ddRdPP*fR + a*ddR_2dPP/Rp + ddQdPP/Rq;
    double dfRdP = fRR*dRdP;
    double ddfRdPP = (2*a/Rp)*(-kappa_2*(-ddrhodPP));
    double dfR_2dP = 2*fR*dfRdP;
    double ddfR_2dPP = 2*(dfRdP*dfRdP + fR*ddfRdPP);

    double lambda=pow(kappa_2*press/pow(clight,2) + f/2.0 + pow(fR,2)/(8*fQ),0.5);
    double lambda_2=pow(lambda,2);

    double sigma_1=fR*0.5+ p_m*pow(2*fQ,0.5)*pow(lambda_2-kappa_2*(rho+press/pow(clight,2)),0.5);
    double sigma_2=fR*0.5+pow(2*fQ,0.5)*lambda;

    double S=pow(sigma_2,2)/pow(sigma_1*sigma_2,0.5);
    double Omega=pow(sigma_1*sigma_2,0.5);

    double Ar=1-2*m*ggrav/(r*pow(clight,2));

    double det_sigma=sigma_1*pow(sigma_2,3);
    double sqrt_det_sigma=pow(det_sigma,0.5);

    double tau_rr=(1/sqrt_det_sigma)*(f*0.5+kappa_2*press/pow(clight,2))*fR;
    double tau_tt=(1/sqrt_det_sigma)*(f*0.5+kappa_2*(-rho)) * fR;

    double dLambda_2dP = kappa_2/pow(clight,2) + dfdP*0.5 + (1/(8*fQ))*dfR_2dP;
    double dLambdadP = 0.5*pow(lambda_2,-0.5)*dLambda_2dP;

    double dSigma1dP = dfRdP*0.5 + pow(2*fQ,0.5)*pow(lambda_2-kappa_2*(rho+press/pow(clight,2)),-0.5)*0.5*(dLambda_2dP-kappa_2*(drho_dp[num_an](press)    + 1/pow(clight,2)));
    double dSigma2dP = 0.5*dfRdP+pow(2*fQ,0.5)*dLambdadP;

    double  dOmegadP = 0.5*pow(sigma_1*sigma_2,-0.5)*sigma_2*dSigma1dP + 0.5*pow(sigma_1*sigma_2,-0.5)*sigma_1*dSigma2dP;
    dOmegadP = pow(sigma_1*sigma_2,-0.5) * 0.5 * (dSigma1dP*sigma_2 + sigma_1* dSigma2dP);

    double dSdSigma1 = -0.5*pow(sigma_2,2)*pow(sigma_1*sigma_2,-3./2.)*sigma_2;
    //double dSdSigma2 = 3*sigma_2*sigma_2/(2*sigma_1)*pow(sigma_2*sigma_2*sigma_2/sigma_1,-0.5);
    double dSdSigma2 = 2*sigma_2*pow(sigma_1*sigma_2,-0.5) + pow(sigma_2,2)*(-0.5)*pow(sigma_1*sigma_2,-1.5)*sigma_1;
    double dSdP = dSdSigma1*dSigma1dP + dSdSigma2*dSigma2dP;
    dSdP = (2*sigma_2*dSigma2dP/(pow(sigma_1*sigma_2,0.5) ) )   - 0.5 * sigma_2*sigma_2*pow(sigma_1*sigma_2,-1.5)* (dSigma1dP*sigma_2 + sigma_1* dSigma2dP);     //skip dSdSigmaversion

    double alpha_r=(rho+press/pow(clight,2))*0.5*(dOmegadP/Omega+dSdP/S);
    double beta_r= (2*r)*dOmegadP/Omega*  (1/pow(clight,2)-(rho+press/pow(clight,2))*0.5*((3/2)*(dOmegadP/Omega)-(dOmegadP/Omega-dSdP/S) ));

    double P_r_0 = ggrav*((rho+press/pow(clight,2))/(r*(r-2*ggrav*m/pow(clight,2))))*(m-(tau_rr+(Omega/S)*tau_tt)*pow(r,3)*0.25*pow(clight,2)/ggrav);
    //double P_r_0 = ((rho+press/pow(clight,2))/(r*(r-2*ggrav*m/pow(clight,2))))*(m*ggrav/pow(clight,2)-(tau_rr+(Omega/S)*tau_tt)*pow(r,3)*0.25)*pow(clight,2);

    P_r_0 =  ((rho+press/pow(clight,2)) )/(r*(r-2*ggrav*m/pow(clight,2)))*(m*ggrav/pow(clight,2) -  (f+kappa_2*(press/pow(clight,2) - rho))/(fR)*pow(r,3)*0.25   )*pow(clight,2);
    alpha_r = (rho+press/pow(clight,2) ) * dfRdP/fR;
    beta_r = 2*r*dfRdP/fR * ( 1/pow(clight,2) - (3*(rho+press/pow(clight,2))/4.0)*dfRdP/fR );




    /////////
    //dPdr
    /////////
    double dPdr= -(P_r_0/(1-alpha_r))*2/(1+p_m*pow(1-beta_r*P_r_0,0.5));


/*
    ///////////////
    // c_1 + c_2*pr + c_3*pr2^tactics
    ///////////////

    double c_1  = 2*(rho+press/pow(clight,2) )/(r*(r-2.0*ggrav*m/pow(clight,2))) * ( m*ggrav/pow(clight,2) - ((tau_rr+(Omega/S)*tau_tt) )*pow(r,3)*0.25  );
    double c_2  = r*2.0/r*(1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*(dOmegadP/Omega + dSdP/S)   );
    double c_3  = r*dOmegadP/Omega*( 1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*(1.5*dOmegadP/Omega - (dSdP/S - dOmegadP/Omega)  )  );

    //c_3= dfRdP/fR * ( 1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*1.5*dfRdP/fR);
    //c_1  = 2*(rho+press/pow(clight,2) )/(r*(r-2.0*ggrav*m/pow(clight,2))) * ( m*ggrav/pow(clight,2) - ((f+kappa_2*(press/pow(clight,2) - rho) ) )*pow(r,3)*0.25  );
    //c_2 = r*2.0/r*(1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*(2*dfRdP/fR)   );
    //c_3 = r*dOmegadP/Omega*( 1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*(1.5*dfRdP/fR  )  );
    // c
    // b
    // a
    if(c_3>pow(10,-30)){dPdr = (-c_2 + p_m*sqrt( pow(c_2,2) - 4*c_3*c_1) )/2*c_3; }
    else{dPdr = -c_1/c_2; }


*/



    ///////////////////
    //dMdr PREPARATIoN
    ///////////////////
    double dOmegadr = dOmegadP*dPdr;
    double dSdr = dSdP*dPdr;


        //New
        double ddSigma1drdP;
        double ddSigma2drdP;

        double ddLambda_2dPP = ddfdPP*0.5 + (1/(8*fQ))*ddfR_2dPP;
        double ddLambda_PP = -0.25*pow(lambda_2,-3./2.)*dLambda_2dP + 0.5*pow(lambda_2,-0.5)*ddLambda_2dPP;

        ddSigma1drdP = pow(2*fQ,0.5)*(-0.25)*pow(lambda_2-kappa_2*(rho+press/pow(clight,2)),-1.5)*(dLambda_2dP*dPdr-kappa_2*dPdr) +
                        pow(2*fQ,0.5)*pow(lambda_2-kappa_2*(rho+press/pow(clight,2)),-0.5) * ddLambda_2dPP*dPdr;

        ddSigma2drdP = pow(2*fQ,0.5)*0.5*(-0.5*pow(lambda_2,-1.5)*dLambda_2dP*dPdr*dLambda_2dP  + 0.5*pow(lambda_2,-0.5)*(pow(kappa_2,2)/pow(clight,2))*dPdr    );



        double ddOmegadrdP = -0.25*pow(sigma_1*sigma_2,-1.5)*(sigma_1*dSigma2dP*dPdr + sigma_2*dSigma1dP*dPdr) * (sigma_1*dSigma2dP + sigma_2*dSigma1dP)
                    + 0.5*pow(sigma_1*sigma_2,-0.5)*( dSigma1dP*dPdr*dSigma2dP + sigma_1*ddSigma2drdP + dSigma2dP*dPdr*dSigma1dP + sigma_2*ddSigma1drdP  );



        double hh = pow(lambda_2 - kappa_2*(rho + press/pow(clight,2)),-0.5);
        double gg = dLambda_2dP - kappa_2*(drhodP + 1/pow(clight,2));
        double dhhdP = -0.5*pow(lambda_2 - kappa_2*(rho + press/pow(clight,2)),-1.5)* (dLambda_2dP - kappa_2*(drhodP + 1/pow(clight,2)) );   //pow(hh,-1)
        double dggdP = ddLambda_2dPP - kappa_2*ddrhodPP;
        double ddSigma1dPP = ddfRdPP*0.5 + pow(2*fQ,0.5)*(dhhdP*gg+hh*dggdP); //dhhdP*gg+
        double ddSigma2dPP = ddfRdPP*0.5 + pow(2*fQ,0.5)*ddLambda_PP;

        double dOmegadSigma1 = 0.5*pow(sigma_1*sigma_2,-0.5)*sigma_2;
        double dOmegadSigma2 = 0.5*pow(sigma_1*sigma_2,-0.5)*sigma_1;
        double ddOmegadSigma1_2 = -0.25*pow(sigma_1*sigma_2, 3./2.)*sigma_2*sigma_2;
        double ddOmegadSigma2_2 = -0.25*pow(sigma_1*sigma_2, 3./2.)*sigma_1*sigma_1;

        double ddOmegadPP = -0.25*pow(sigma_1*sigma_2,-1.5)*pow( sigma_1*dSigma2dP + sigma_2*dSigma1dP,2)    +  0.5*pow(sigma_1*sigma_2,-0.5)*(sigma_1*ddSigma2dPP + 2*dSigma1dP*dSigma2dP + sigma_2*ddSigma1dPP);




//________________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________________


//FROM HERE TO NEXT DOUBLE LINE WAS ATTEMPT TO DERIVE omega_rr and P_rr myself: but did mistakes by not including all derivatives: e.g. of rho

        double ddSdSigma1dSigma2 = 3*sigma_2*sigma_2*(-0.5)*pow(sigma_1*sigma_2,-3./2.) + pow(sigma_2,3)*(-0.5)*(-1.5)*pow(sigma_1*sigma_2,-5./2)*sigma_1;
        double ddSdSigma1Sigma1 = pow(sigma_2,4)*(-0.5)*(-1.5)*pow(sigma_1*sigma_2,-5./2);
        double ddSdSigma2Sigma2 = -0.25*pow(sigma_2*sigma_2*sigma_2/sigma_1,-1.5)*3*sigma_2*sigma_2/sigma_1 + 0.5*pow(sigma_2*sigma_2*sigma_2/sigma_1,-0.5)*6*sigma_2/sigma_1;


        //S_PP
        double ddSdPP = dSdSigma1*ddSigma1dPP + dSdSigma2*ddSigma2dPP + ddSdSigma1Sigma1*pow(dSigma1dP,2) + 2*ddSdSigma1dSigma2*dSigma1dP*dSigma2dP + ddSdSigma2Sigma2*pow(dSigma2dP,2);

        //ddSdPP directly from Mathematica
        ddSdPP = (2*pow(dSigma2dP,2) / ( pow(sigma_1*sigma_2, 0.5) )   )
                - (2*sigma_2*dSigma2dP*(sigma_2*dSigma1dP + sigma_1*dSigma2dP) ) / (pow(sigma_1*sigma_2,1.5))
                + 3*sigma_2*sigma_2 * pow(sigma_2*dSigma1dP - sigma_1*dSigma2dP, 2 ) / (pow(4*sigma_1*sigma_2,2.5))
                + 2*sigma_2*ddSigma2dPP / (pow(sigma_1*sigma_2,0.5))
                - ( sigma_2*sigma_2* ( 2*dSigma1dP*dSigma2dP + sigma_2*ddSigma1dPP + sigma_1*ddSigma2dPP )  ) / (2*pow(sigma_1*sigma_2,1.5));




        //ALPHA(r) - BETA(r)
        //often needed derivative of Omega_P/Omega and of S_P/S
        double dSPDSdP = dSdP*((-1)/pow(S,2))*dSdP + ddSdPP/S;
        double dOmegaPDOmegadP = dOmegadP*((-1)/pow(Omega,2))*dOmegadP + ddOmegadPP/Omega;
        double dalpha_rdr = 0.5*(drho_dp[num_an](press)+ 1/pow(clight,2) )*(dOmegadP/Omega + dSdP/S) + (rho + press/pow(clight,2))*0.5*(dOmegaPDOmegadP + dSPDSdP);
        dalpha_rdr = dalpha_rdr*dPdr;
        double beta_r1 = 2*r*dOmegadP/Omega;
        double beta_r2 = 1- (rho + press/pow(clight,2))*0.5*( 1.5* dOmegadP/Omega - (dOmegadP/Omega - dSdP));
        double dbeta_r1dP = 2.0/dPdr*dOmegadP/Omega + 2*r*dOmegaPDOmegadP;
        double dbeta_r2dP = - (drhodP + 1.0/pow(clight,2) )/2.0*( 1.5* dOmegadP/Omega - (dOmegadP/Omega - dSdP))   - (rho + press/pow(clight,2))*0.5* ( 1.5* dOmegaPDOmegadP - (dOmegaPDOmegadP - dSPDSdP)) ;
        double dbeta_rdP = dbeta_rdP + 2*dOmegadP/Omega*(rho+press/pow(clight,2))*0.5*(0.5*dOmegadP/Omega-dSPDSdP/S)
                              + 2*r*dOmegaPDOmegadP*(rho+press/pow(clight,2))*0.5*(0.5*dOmegadP/Omega-dSPDSdP/S)
                              + 2*r*dOmegadP/Omega*0.5*(0.5*dOmegadP/Omega-dSPDSdP/S)
                              + 2*r*dOmegadP/Omega*(rho+press/pow(clight,2))*0.5*(0.5*dOmegaPDOmegadP - dSPDSdP);
        dbeta_rdP = dbeta_r1dP*beta_r2 + dbeta_r2dP*beta_r1;
        double dbeta_rdr = dbeta_rdP*dPdr;


        double dSPdr = ddSdPP*dPdr;
        double dOmegaPdr = ddOmegadPP*dPdr;
        double dSPDSdr = dSPdr/S + dSdP*(-1)*pow(S,-2)*dSdr;
        double dOmegaPDOmegadr = dOmegaPdr/Omega + dOmegadP*(-1)*pow(Omega,-2)*dOmegadr;
        double v2 = 1.5*dOmegadP/Omega - (dOmegadP/Omega + dSdP/S);
        double dv2dr = 1.5*dOmegaPDOmegadr - (dOmegaPDOmegadr + dSPDSdr);
        double v = 1- 0.5*(rho + press/pow(clight,2))*v2;
        double dvdr = - 0.5*(drho_dp[num_an](press)+ dPdr/pow(clight,2)  )*dv2dr - 0.5*(rho + press/pow(clight,2))*v2;


        double d1Dsqrt_det_sigmadP = -0.5*pow(det_sigma,-3./2.)*(pow(sigma_2,3)*dSigma1dP + 3*sigma_1*pow(sigma_2,2)*dSigma2dP);
        double dtau_rrdP = d1Dsqrt_det_sigmadP*(f*0.5+kappa_2*press/pow(clight,2))+ 1/sqrt_det_sigma*(dfdP*0.5+kappa_2/pow(clight,2));
        double dtau_ttdP = d1Dsqrt_det_sigmadP*(f*0.5-kappa_2*rho)+ 1/sqrt_det_sigma*(dfdP*0.5 - kappa_2*drhodP);



        double P_rr_0 = dPdr*(r*(r-2*m*ggrav/pow(clight,2)))*(m*ggrav/pow(clight,2)-(tau_rr+Omega*(1/S)*tau_tt)*pow(r,3)*0.25)
                       +(rho+press/pow(clight,2))*(2*r-2*m*ggrav/pow(clight,2))*(m*ggrav/pow(clight,2)-(tau_rr+Omega*(1/S)*tau_tt)*pow(r,3)*0.25)
                       +(rho+press/pow(clight,2))*(r*(r-2*m*ggrav/pow(clight,2)))*(  (3*r*r/4.)*m*ggrav/pow(clight,2) - dtau_rrdP*dPdr*pow(r,3)/4 - tau_rr*3+pow(r,2)/4
                                                    - dOmegadP*dPdr*(1/S)*tau_tt*pow(r,3)*0.25
                                                    - Omega*(-1)*pow(S,-2)*dSdP*dPdr*tau_tt*pow(r,3)*0.25
                                                    - Omega*(1/S)*dtau_ttdP*dPdr*pow(r,3)*0.25
                                                    - Omega*(1/S)*tau_tt*3*pow(r,2)*0.25);

        double P_rr= -P_rr_0*(1/(1-alpha_r))*2*1/(1+pow(1-beta_r*P_r_0,0.5))
                     -P_r_0*(-1)*pow(1-alpha_r,-2)*dalpha_rdr*2*1/(1+pow(1-beta_r*P_r_0,0.5))
                     -P_rr_0*(1/(1-alpha_r))*2*(-1)*pow(1+pow(1-beta_r*P_r_0,0.5),-2)
                                                   *0.5*pow(1-beta_r*P_r_0,-0.5)
                                                   *(-1)*(dbeta_rdr*P_r_0 + beta_r*P_rr_0);

        double ddOmegadrr = ddOmegadPP*dPdr*dPdr + dOmegadP*P_rr;


//________________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________________



    //Preparing l,j,k ,  a,b,c,d
    double ll=((dOmegadr/Omega) + (2.0/r))/r;
    double jj=0.5*(3*tau_rr - (Omega/S)*tau_tt);
    double kk= (dOmegadr/Omega)*(    ((2*r-3*m*ggrav/pow(clight,2))/(r*(r-2*m*ggrav/pow(clight,2))))-0.75*(dOmegadr/Omega)      );


    //The following variables are used in Helios paper to calculate Prr. But i calculated Prr manually above, but wrong: did not take derivative with respect to r all variables e.g. rho
    double ss=1;
    double aa= 1+                        0.5*ss* ( (beta_r*P_r_0)/(pow(1-beta_r*P_r_0,0.5)*(1+ss*pow(1-beta_r*P_r_0,0.5)))  );
    double bb=(dalpha_rdr/(1-alpha_r)) + 0.5*ss* ( (dbeta_rdr*P_r_0)/(pow(1-beta_r*P_r_0,0.5)*(1+ss*pow(1-beta_r*P_r_0,0.5)))  );

    double big_phi = tau_rr + (Omega/S)*tau_tt;
    double drhodp = drho_dp[num_an](press);
    double dPhidP = dtau_rrdP + dtau_ttdP*(Omega/S) + tau_tt*(dOmegadP/S - Omega/pow(S,2)*dSdP);
    /*
    double cc = ( (1+drhodp)/(rho+press/pow(clight,2)) - (dPhidP*pow(r,3)*0.25)/(m-big_phi * pow(r,3)*0.25)  )*dPdr -
    (  (2*(r-m*ggrav/pow(clight,2)))/(r*(r-2*m*ggrav/pow(clight,2))) + (3*big_phi*r*r*0.25 )/(m-big_phi*pow(r,3)*0.25)   );
    */
    double cc =  ( (1/pow(clight,2)+drhodp)/(rho+press/pow(clight,2)))*dPdr
            - ((dPhidP*pow(r,3)*0.25)/(m*ggrav/pow(clight,2)-big_phi * pow(r,3)*0.25)  )*dPdr
            - (  (2*(r-m*ggrav/pow(clight,2)))/(r*(r-2*m*ggrav/pow(clight,2)))   + (3*big_phi*r*r*0.25 )/(m*ggrav/pow(clight,2)-big_phi*pow(r,3)*0.25)   );

//
     //aa=0;
     double dd =  (2/(r-2*m*ggrav/pow(clight,2)) + 1/(m*ggrav/pow(clight,2)-big_phi*r*r*r*0.25));





    ///////////////////
    //dMdr
    ///////////////////

    //New dMdr
    double new_dmdr = jj + Ar*( (ddOmegadPP*dPdr*dPdr + dOmegadP*   dPdr*(cc*aa+bb))/Omega + kk   ) ;
    new_dmdr = new_dmdr/(ll- (Ar*dOmegadP*dPdr*dd*aa)/Omega);
    double dm_dr = 4*pi*rho*pow(r,2);

    if(r>1*pow(10,-6))
    {
        dm_dr = new_dmdr*pow(clight,2)/ggrav;
        dm_dr = get_gradients2[1](m, press, r).first;
    }


    //HERE THE GRADIENTS ARE FINALLY CALCULATED WITH ALL THE PRECALCULATED TERMS
    double dp_dr;

    if(r>1*pow(10,-6))
    {
        dp_dr = dPdr;
        //dp_dr = get_gradients2[1](m, press, r).second;
    }
    else
    {
        dp_dr = -ggrav*( rho+press/pow(clight,2) )*(4.0*pi*r*press/pow(clight,2));
    }



    return make_pair(dm_dr, dp_dr);
}
