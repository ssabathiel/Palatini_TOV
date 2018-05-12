
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
            double Q=(3*pow(Rq,2)/8)*(1-(2*kappa_2*(rho+press/pow(clight,2))/Rq)+(2*pow(kappa_2,2)*pow((rho-3*press/pow(clight,2)),2)/(3*pow(Rq,2)))-pow(1-(4*kappa_2*(rho+press/pow(clight,2)))/Rq,0.5));     //make Energy-dependent??

            Q = -(kappa_2*press/pow(clight,2) + f_tilde*0.5 + (Rp/(8*b)) )
            + Rp/(32*b) * pow(    3*((b*R)/(Rp) + f_tildeR) - pow( pow( (b*R)/(Rp) + f_tildeR,2 ) -(4*b*kappa_2*(rho+press/pow(clight,2) ) )/(Rp), 0.5  ),2); //From Helios--> Olmo --> p/m-paper


            double f=R+a*pow(R,2)/Rp + Q/Rq;

            double dRdP = -kappa_2*(   -drhodP + (3/pow(clight,2) )    );
            double dR_2dP = 2*R*dRdP;
            double ddRdPP = kappa_2*ddrhodPP;
            double ddR_2dPP = 2*R*dRdP; //18* kappa_2;
            double fR=1+2*a*R/Rp;
            double fQ=1/Rq;
            //double dfdRP = 1 + 2*a*R/Rp;
            double dQdP = (3*pow(Rq,2)/8)*(-2*kappa_2*(drhodP + 1/pow(clight,2))/Rq + (2*pow(kappa_2,2)*2*(rho-3*press/pow(clight,2))*(drhodP-3/pow(clight,2)))/(3*pow(Rq,2))-0.5*pow(1-4*kappa_2*(rho+press/pow(clight,2))/Rq,-0.5)*(-4*kappa_2*(drhodP + 1/pow(clight,2))/Rq));
            //ddQdPP = (3*pow(Rq,2)/8) * ( 36*pow(kappa_2,2)/(3*Rq*Rq) +0.25*pow(1-4*kappa_2*(rho+press/pow(clight,2))/Rq,-3./2.)*pow(-4*kappa_2/Rq,2) );
            double ddQdPP = 3*Rq*Rq*(1./8.)* (  - (2*kappa_2*(ddrhodPP))/(Rq)
                                        + 2*kappa_2*2*( (drhodP - 3/pow(clight,2))*(drhodP - 3/pow(clight,2))* (rho-3*press/pow(clight,2))*ddrhodPP   )/(3*Rq*Rq)
                                        - 0.5*(  -0.5*pow(1-(4*kappa_2*(rho+press/pow(clight,2) ))/(Rq), -1.5 ) * (-4*kappa_2*(drhodP + 1/pow(clight,2)/Rq))* (-4*kappa_2*(drhodP + 1/pow(clight,2)/Rq))
                                                  + (pow(1-(4*kappa_2*(rho+press/pow(clight,2) ))/(Rq), -0.5 ))*(4*kappa_2*(ddrhodPP)/Rq)   )
                                    );

            double dfdP = dRdP + a*dR_2dP/Rp + dQdP/Rq;
            double ddfdPP = ddRdPP + a*ddR_2dPP/Rp + ddQdPP/Rq;
            double dfRdP = (2*a/Rp)*dRdP;					// multiplied Rp artificially!!!!!!!!!!!!!!!!
            double ddfRdPP = (2*a/Rp)*(-kappa_2*(-ddrhodPP));
            double dfR_2dP = 2*fR*dfRdP; //4*a*dRdP/Rp + 4*pow(a,2)*dR_2dP/pow(Rp,2);
            double ddfR_2dPP = 2*(dfRdP*dfRdP + fR*ddfRdPP); //4*a*a/pow(Rp,2) * kappa_2*kappa_2*18;

            double lambda=pow(kappa_2*press/pow(clight,2) + f/2 + pow(fR,2)/(8*fQ),0.5);
            double lambda_2=pow(lambda,2);

            double sigma_1=fR*0.5+ p_m*pow(2*fQ,0.5)*pow(lambda_2-kappa_2*(rho+press/pow(clight,2)),0.5);
            double sigma_2=fR/2+pow(2*fQ,0.5)*lambda;

            double S=pow(sigma_2,2)/pow(sigma_1*sigma_2,0.5);
            double Omega=pow(sigma_1*sigma_2,0.5);

            double Ar=1-2*m*ggrav/(r*pow(clight,2));

            double det_sigma=sigma_1*pow(sigma_2,3);           //1/math.sqrt(sigma1 * sigma2 * sigma2 * sigma2)
            double sqrt_det_sigma=pow(det_sigma,0.5);

            double tau_rr=(1/sqrt_det_sigma)*(f*0.5+kappa_2*press/pow(clight,2));    //T_rr=pressure for a perfect fluid
            double tau_tt=(1/sqrt_det_sigma)*(f*0.5+kappa_2*(-rho));     //press cancel each other

            double dLambda_2dP = kappa_2/pow(clight,2) + dfdP*0.5 + (1/(8*fQ))*dfR_2dP;
            double dLambdadP = 0.5*pow(lambda_2,-0.5)*dLambda_2dP;


            double dSigma1dP = dfRdP*0.5 + pow(2*fQ,0.5)*pow(lambda_2-kappa_2*(rho+press/pow(clight,2)),-0.5)*0.5*(dLambda_2dP-kappa_2*(drho_dp[num_an](press)    + 1/pow(clight,2)));
            double dSigma2dP = 0.5*dfRdP+pow(2*fQ,0.5)*dLambdadP;

            double dOmegadP = 0.5*pow(sigma_1*sigma_2,-0.5)*sigma_2*dSigma1dP + 0.5*pow(sigma_1*sigma_2,-0.5)*sigma_1*dSigma2dP;
            dOmegadP = pow(sigma_1*sigma_2,-0.5) * 0.5 * (dSigma1dP*sigma_2 + sigma_1* dSigma2dP);

            double dSdSigma1 = -0.5*pow(sigma_2,2)*pow(sigma_1*sigma_2,-3./2.)*sigma_2;
            //double dSdSigma2 = 3*sigma_2*sigma_2/(2*sigma_1)*pow(sigma_2*sigma_2*sigma_2/sigma_1,-0.5);
            double dSdSigma2 = 2*sigma_2*pow(sigma_1*sigma_2,-0.5) + pow(sigma_2,2)*(-0.5)*pow(sigma_1*sigma_2,-1.5)*sigma_1;
            double dSdP = dSdSigma1*dSigma1dP + dSdSigma2*dSigma2dP;
            dSdP = (2*sigma_2*dSigma2dP/(pow(sigma_1*sigma_2,0.5) ) )   - 0.5 * sigma_2*sigma_2*pow(sigma_1*sigma_2,-1.5)* (dSigma1dP*sigma_2 + sigma_1* dSigma2dP);     //skip dSdSigmaversion
            //for dOmegadR=dOmegadP * dPdR    need P_r first

            //ALPHA_BETA
            double alpha_r=(rho+press/pow(clight,2))*0.5*(dOmegadP/Omega+dSdP/S);
            double beta_r= (2*r)*dOmegadP/Omega*  (1-(rho+press/pow(clight,2))*0.5*((3/2)*(dOmegadP/Omega)-(dOmegadP/Omega-dSdP/S) ));



            //P_r_0
            double P_r_0 = ggrav*((rho+press/pow(clight,2))/(r*(r-2*ggrav*m/pow(clight,2))))*(m-(tau_rr+(Omega/S)*tau_tt)*pow(r,3)*0.25*pow(clight,2)/ggrav);




            /////////
            //dPdr
            /////////
            double dPdr= -(P_r_0/(1-alpha_r))*2/(1+p_m*pow(1-beta_r*P_r_0,0.5));







            ///////////////////
            //dMdr PREPARATIoN
            ///////////////////
            double dOmegadr = dOmegadP*dPdr;

            //cout << dOmegadr/Omega << endl;
            double dSdr = dSdP*dPdr;

            //Needs some preparation: dMdr ~ Omega_rr ~ Omega_PP, P_rr
                //Omega_PP
                //needs many second derivatives with resp. to P

                //New
                double ddSigma1drdP;
                double ddSigma2drdP;

                double ddLambda_2dPP = ddfdPP*0.5 + (1/(8*fQ))*ddfR_2dPP;
                double ddLambda_PP = -0.25*pow(lambda_2,-3./2.)*dLambda_2dP + 0.5*pow(lambda_2,-0.5)*ddLambda_2dPP;

                double dddSigma1drdP = pow(2*fQ,0.5)*(-0.25)*pow(lambda_2-kappa_2*(rho+press/pow(clight,2)),-1.5)*(dLambda_2dP*dPdr-kappa_2*dPdr) +
                                pow(2*fQ,0.5)*pow(lambda_2-kappa_2*(rho+press/pow(clight,2)),-0.5) * ddLambda_2dPP*dPdr;

                ddSigma2drdP = pow(2*fQ,0.5)*0.5*(-0.5*pow(lambda_2,-1.5)*dLambda_2dP*dPdr*dLambda_2dP  + 0.5*pow(lambda_2,-0.5)*(pow(kappa_2,2)/pow(clight,2))*dPdr    );



                double ddOmegadrdP = -0.25*pow(sigma_1*sigma_2,-1.5)*(sigma_1*dSigma2dP*dPdr + sigma_2*dSigma1dP*dPdr) * (sigma_1*dSigma2dP + sigma_2*dSigma1dP)
                            + 0.5*pow(sigma_1*sigma_2,-0.5)*( dSigma1dP*dPdr*dSigma2dP + sigma_1*ddSigma2drdP + dSigma2dP*dPdr*dSigma1dP + sigma_2*ddSigma1drdP  );








                               //double dLambdadP = 0.5*pow(lambda_2,-0.5)*dLambda_2dP;
                //mistaken?
                //ddSigma1dPP = pow(2*fQ,0.5)*(-0.5*pow(lambda_2-kappa_2*(rho + press/pow(clight,2)),-3./2.)*pow(dLambda_2dP-kappa_2,2.)+(1/(lambda_2-kappa_2*(rho+press/pow(clight,2))))*ddLambda_2dPP);
                //ddSigma2dPP = pow(2*fQ,0.5)*ddLambda_PP;

                double hh = pow(lambda_2 - kappa_2*(rho + press/pow(clight,2)),-0.5);
                double gg = dLambda_2dP - kappa_2*(drhodP + 1/pow(clight,2));
                double dhhdP = -0.5*pow(lambda_2 - kappa_2*(rho + press/pow(clight,2)),-1.5)* (dLambda_2dP - kappa_2*(drhodP + 1/pow(clight,2)) );   //pow(hh,-1)
                double dggdP = ddLambda_2dPP - kappa_2*ddrhodPP;
                double ddSigma1dPP = ddfRdPP*0.5 + pow(2*fQ,0.5)*(dhhdP*gg+hh*dggdP); //dhhdP*gg+
                double ddSigma2dPP = ddfRdPP*0.5 + pow(2*fQ,0.5)*ddLambda_PP;

                //cout << "ddfRdPP= " << ddfRdPP << endl;



                double dOmegadSigma1 = 0.5*pow(sigma_1*sigma_2,-0.5)*sigma_2;
                double dOmegadSigma2 = 0.5*pow(sigma_1*sigma_2,-0.5)*sigma_1;
                double ddOmegadSigma1_2 = -0.25*pow(sigma_1*sigma_2, 3./2.)*sigma_2*sigma_2;
                double ddOmegadSigma2_2 = -0.25*pow(sigma_1*sigma_2, 3./2.)*sigma_1*sigma_1;

                //ddSigma1dPP = pow(2*fQ,0.5)*(-0.5*pow(lambda_2-kappa_2*(rho + press/pow(clight,2)),-3./2.)*pow(dLambda_2dP-kappa_2,2.)+(1/(lambda_2-kappa_2*(rho+press/pow(clight,2))))*ddLambda_2dPP);
                //ddSigma2dPP = pow(2*fQ,0.5)*ddLambda_PP;


                //OMEGA_PP
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
                double dbeta_rdr = 2*dOmegadP/Omega + 2*r*dOmegaPDOmegadP;
                dbeta_rdr = dbeta_rdr + 2*dOmegadP/Omega*(rho+press/pow(clight,2))*0.5*(0.5*dOmegadP/Omega-dSPDSdP/S)
                                      + 2*r*dOmegaPDOmegadP*(rho+press/pow(clight,2))*0.5*(0.5*dOmegadP/Omega-dSPDSdP/S)
                                      + 2*r*dOmegadP/Omega*0.5*(0.5*dOmegadP/Omega-dSPDSdP/S)
                                      + 2*r*dOmegadP/Omega*(rho+press/pow(clight,2))*0.5*(0.5*dOmegaPDOmegadP - dSPDSdP);
                dbeta_rdr = dbeta_rdr*dPdr;




                //now directly dbeta_rdr: since problem which was not considered here: drdP

                double dSPdr = ddSdPP*dPdr;
                double dOmegaPdr = ddOmegadPP*dPdr;
                double dSPDSdr = dSPdr/S + dSdP*(-1)*pow(S,-2)*dSdr;
                double dOmegaPDOmegadr = dOmegaPdr/Omega + dOmegadP*(-1)*pow(Omega,-2)*dOmegadr;
                double v2 = 1.5*dOmegadP/Omega - (dOmegadP/Omega + dSdP/S);
                double dv2dr = 1.5*dOmegaPDOmegadr - (dOmegaPDOmegadr + dSPDSdr);
                double v = 1- 0.5*(rho + press/pow(clight,2))*v2;
                double dvdr = - 0.5*(drho_dp[num_an](press)+ dPdr/pow(clight,2)  )*dv2dr - 0.5*(rho + press/pow(clight,2))*v2;

                //dbeta_rdr = 2*(dOmegadP/Omega)*v + 2*r*dOmegaPDOmegadr*v + 2*r*(dOmegadP/Omega)*dvdr;


                //dtau_rrdr, dtau_ttdr
                double d1Dsqrt_det_sigmadP = -0.5*pow(det_sigma,-3./2.)*(pow(sigma_2,3)*dSigma1dP + 3*sigma_1*pow(sigma_2,2)*dSigma2dP);
                double dtau_rrdP = d1Dsqrt_det_sigmadP*(f*0.5+kappa_2*press)+ 1/sqrt_det_sigma*(dfdP*0.5+kappa_2);
                double dtau_ttdP = d1Dsqrt_det_sigmadP*(f*0.5-kappa_2*rho)+ 1/sqrt_det_sigma*(dfdP*0.5);



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

                //ddOmegadrr =  ddOmegadrdP*dPdr + dOmegadP*P_rr; //ddOmegadPP*dPdr*dPdr + dOmegadP*P_rr;
                //cout << ddOmegadrr/Omega << endl;
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
            double bb=(dalpha_rdr/(1-alpha_r)) + 0.5*ss* ( (dbeta_rdr*P_r_0)/(pow(1-beta_r*P_r_0,0.5)*(1+ss*pow(1-beta_r*P_r_0,0.5)))  );   //same 2nd term as in a

            double big_phi = (tau_rr + (Omega/S)*tau_tt);
            // Prepare drhodp: therefore need parameter for polytrope. later put this into function
            double gamma_now;
            double K_now;
            if(rho<rho_1)
            {
                gamma_now=gamma_0;
                K_now = K_0;
            }
            else
            {
                gamma_now=gamma_1;
                K_now = K_1;
            }
            double drhodp = drho_dp[num_an](press);//(1.0/gamma_now)*pow(press/K_now, (1.0/gamma_now)-1 )*(1.0/K_now);
            double dPhidP = dtau_rrdP + dtau_ttdP*(Omega/S) + tau_tt*(dOmegadP/S - Omega/pow(S,2)*dSdP);
            double cc = ( (1+drhodp)/(rho+press/pow(clight,2)) - (dPhidP*pow(r,3)*0.25)/(m-big_phi * pow(r,3)*0.25)  )*dPdr -
            (  (2*(r-m*ggrav/pow(clight,2)))/(r*(r-2*m*ggrav/pow(clight,2))) + (3*big_phi*r*r*0.25 )/(m-big_phi*pow(r,3)*0.25)   );
            double dd = (2/(r-2*m*ggrav/pow(clight,2)) + 1/(m-big_phi*r*r*r*0.25));







            ///////////////////
            //dMdr
            ///////////////////

            //OLD dMdr
            //double dmdr= ((3*tau_rr-(Omega/S)*tau_tt)*0.5 + Ar*( ddOmegadrr/Omega + (dOmegadr/Omega)*( (2*r-3*m)/(r*(r-2*m)) - 3*0.25*dOmegadr/Omega    ) )) * pow((dOmegadr/Omega) + 2/r,-1) * r;
            double dmdr= (3*tau_rr-(Omega/S)*tau_tt)*0.5;// + Ar*( ddOmegadrr/Omega + (dOmegadr/Omega)*( (2*r-3*m*ggrav/pow(clight,2))/(r*(r-2*m*ggrav/pow(clight,2))) - 3*0.25*dOmegadr/Omega    ) );

            dmdr=dmdr*pow((dOmegadr/Omega) + 2/r,-1) * r;

            //New dMdr
            double new_dmdr = jj + Ar*( (ddOmegadPP*dPdr*dPdr + dOmegadP*   dPdr*(cc*aa+bb))/Omega + kk   ) ;
            new_dmdr = new_dmdr/(ll- (Ar*dOmegadP*dPdr*dd*aa)/Omega);
            double dm_dr = 4*pi*rho*pow(r,2);

            if(r>1*pow(10,-6))
            {

            //dm_dr = dmdr*pow(clight,2)/ggrav;
            dm_dr = new_dmdr*pow(clight,2)/ggrav;

            }


            double dp_dr;

            //Get dp/dr = ...TOV...
            if(r>1*pow(10,-6))
            {
                dp_dr = dPdr;
            }
            else
            {
                dp_dr = (4.0*pi*r*press/pow(clight,2));
                dp_dr = -ggrav*( rho+press/pow(clight,2) )*dp_dr;
            }



        //////////////////////////////
        // fR-Code-old
        ////////////////////////////

            f=R+a*pow(R,2)/Rp;
            double capF = 1 + 2*a*R/Rp;
            double dcapFdr = (2*a/Rp)*(-kappa_2*(-drhodP));
            double ddcapFdrr = 0;                               //Assumed that 2nd deriv. of rho is 0. can be done numerically as well.
            double R_rho = kappa_2;
            double f_rho = R_rho + 2*a*R*R_rho/Rp;
            double capF_rho = 2*a*R_rho/Rp;
            double lambda_rho = 0.5*(R_rho - (f_rho/capF_rho));
            double e_B = 1/(1-2*m*ggrav/(r*pow(clight,2)) - lambda_rho*pow(r,2)*ggrav/(3.0*pow(clight,2)));
            double alpha = pow(r,2)*( (3.0/4.0)* pow(dcapFdr/capF,2) + (2*dcapFdr)/(r*capF) + (e_B)/(2)*(R-f/capF)    );
            double beta = pow(r,2)* ( (ddcapFdrr/capF) - (3.0/2.0)*pow(dcapFdr/capF,2)  );
            double gamma = r*dcapFdr/(2*capF);
            double dpdr = -(1.0)/(1.0+gamma) * ggrav*((rho+press/pow(clight,2))/(r*(r-2*ggrav*m/pow(clight,2))))
                    *( ggrav*m/pow(clight,2) + (4*pi*ggrav*pow(r,3)*press/pow(clight,4))/(capF) - 0.5*alpha*(r-2*ggrav*m/pow(clight,2)))*(pow(clight,2)/ggrav);
            dmdr = (1.0)/(1.0+gamma) * ( (4*pi*pow(r,2)*rho)/(capF) + (alpha+beta)*pow(clight,2)/(2*ggrav)  - (m)/(r)*(alpha+beta-gamma)   );




            double dp_dr_fR;
            double dm_dr_fR;

            if(r>1*pow(10,-6))     // does not matter how I cange this value: always at this value dm_dr gets negative
            {
                dp_dr_fR = dpdr;
                dm_dr_fR = dmdr;
                //dp_dr = (m+4.0*pi*r*r*r*press/pow(clight,2))/(r*(r-2.0*ggrav*m/pow(clight,2)));
            }
            else
            {
                dp_dr_fR = (4.0*pi*r*press/pow(clight,2));
                dp_dr_fR = -ggrav*( rho+press/pow(clight,2) )*dp_dr_fR;
                dm_dr_fR = 4*pi*rho*pow(r,2);
            }






        //////////////////////////////
        // fR-limit from Helios Paper
        ////////////////////////////


        //dp_dr_limt
        P_r_0 = ggrav*((rho+press/pow(clight,2))/(r*(r-2*ggrav*m/pow(clight,2))))*(m-(tau_rr+tau_tt)*pow(r,3)*0.25*pow(clight,2)/ggrav);
        alpha_r = (rho+press/pow(clight,2))*dfRdP/fR;
        beta_r = 2*r*(dfRdP/fR)*(1-(rho+press/pow(clight,2))*(3*dfRdP)/(4*fR) );
        double dp_dr_lim = - (P_r_0*2)/( (1-alpha_r)*( 1+ pow(1-beta_r*P_r_0,0.5) )     );


        //dm_dr_limt
        double dfRdr = (2*a)/(Rp)*kappa_2*(drhodP*dp_dr_lim - 3*dp_dr_lim/pow(clight,2) );
        double xx = ((dfRdr/fR) + (2.0/r) )*(1.0/r);
        double yy = (3*tau_rr - tau_tt)/(2.0);
        double zz = (dfRdr/fR)*( (2*r-3*ggrav*m/pow(clight,2))/(r*(r-2*ggrav*m/pow(clight,2))   )  - 0.75*(dfRdr)/(fR)   );

        double ww = Ar*(1/fR)*kappa_2*2*a*(1/Rp)*(3/pow(clight,2)-drhodP)*dd*aa*dp_dr_lim;     //check out sign for ww and vv
        // dOmega
        double ddrhodrr = ddrhodPP*pow(dp_dr_lim,2) + drhodP*P_rr;
        ddrhodrr=0;                                                                                     //MODIFIED 2nd derivative
        double vv = -Ar*( (1/fR)*kappa_2*(2*a/Rp)*(-ddrhodrr  + 3*(cc*aa+bb)*dp_dr_lim/pow(clight,2) )   );
        //double vv = Ar*( (1/fR)*dOmegadP*(-ddrhodrr  + 3*(cc*aa+bb)*dp_dr_lim)/pow(clight,2)   );
        //cout <<<< Ar*( (ddOmegadPP*dPdr*dPdr + dOmegadP*   dPdr*(cc*aa+bb))/Omega) << endl;

        double dm_dr_lim = (yy+vv+Ar*zz)/(xx+ww)*( (pow(clight,2))/(ggrav) );



        //////////////////////////////
        // fR-limit from Helios Paper ENDE
        ////////////////////////////













        //////////////////////////////
        // fR - original
        ////////////////////////////



            double dRdT = -kappa_2;
            double dTdP = -drhodP + 3.0/pow(clight,2);
            dRdP = dRdT*dTdP;
            double fRR = 2*a/Rp;
            double ddfRdRR = 0;
            double ddRdTT = 0;
            double dTdrho = -1;

            double dfRdT = fRR*dRdT;
            dfRdP = dfRdT*dTdP;

            capF = fR;

            R_rho = dRdT*dTdrho;
            f_rho = fR*R_rho;
            capF_rho = fRR*R_rho;
            lambda_rho = 0.5*(R_rho - (f_rho/capF_rho));
            e_B = 1/(1-2*m*ggrav/(r*pow(clight,2)));// - lambda_rho*pow(r,2)*ggrav/(3.0*pow(clight,2)));

            double N_1 = (1.0/fR)*fRR*dRdT*(3/pow(clight,2)-drhodP);
            double N_2 = (1.0/fR)*( ddfRdRR*pow(dRdT,2) + fRR*ddRdTT)*pow(3.0/pow(clight,2)-drhodP,2) - (1.0/fR)*fRR*dRdT*ddrhodPP;

            aa = r*N_1*( 1/pow(clight,2)-0.75*(rho+press/pow(clight,2) )*N_1 );
            bb = 2*(1/pow(clight,2)-(rho+press/pow(clight,2) )*N_1 );
            cc = -(rho+press/pow(clight,2) )*( (1-e_B)/(r) - ((8*pi*r*e_B*press*ggrav/pow(clight,4))/fR) + (r*e_B)/(2*fR)*( fR*R-f ) );

            double dPdr_GR=-ggrav*( rho+press/pow(clight,2) )*(m+4.0*pi*r*r*r*press/pow(clight,2))/(r*(r-2.0*ggrav*m/pow(clight,2)));

            //double dPdr;
            if(aa>pow(10,-20)){dPdr = (-bb + p_m*sqrt( pow(bb,2) - 4*aa*cc) )/2*aa; }
            else{dPdr = -cc/bb; }



            /////////
            /// d/dr - Derivatives (now possible since we have dPdr)
            ///

            dcapFdr = dfRdP*dPdr;
            double dRdr = dRdT*dTdP*dPdr;
            double drhodr = drhodP*dPdr;

            /////////
            /// 2nd Derivatives
            ///
            //double ddTdrr = -ddrhodrr+3*ddPdrr/pow(clight,2);
            double ddTdPP = -ddrhodPP;

            ddRdPP = ddRdTT*pow(dTdP,2) + dRdT*ddTdPP;
            double ddcapFdPP = ddfRdRR*pow(dRdP,2) + fRR*ddRdPP;

            alpha = pow(r,2)*( (3.0/4.0)* pow(dcapFdr/capF,2) + (2*dcapFdr)/(r*capF) + (e_B)/(2)*(R-f/capF)    );

            gamma = r*dcapFdr/(2*capF);
            double A_prime = -(1/(1+gamma))*( (1-e_B)/(r) - e_B*8*pi*r*press*(ggrav/pow(clight,4))/fR   + alpha/r );
            dd = 0.5*r*A_prime + 2*gamma +1;


            double ddPdrr = (N_2*pow(dPdr,2) - 0.5*A_prime*( (drhodr + dPdr/pow(clight,2))/(rho+press/pow(clight,2)) - dd/r)
                    - dd/(2*(1+gamma)) * ( N_2*pow(dPdr,2) + (1-e_B)/(r*r) + e_B/fR*(fR*R-f) + 8*pi*e_B*rho*(ggrav/pow(clight,2))/fR + (gamma*(4-3*gamma))/(r*r) )
                    +(e_B)/(2*fR)*( fR*R-f) -  8*pi*e_B*press*(ggrav/pow(clight,4))/fR + (gamma*(2*3*gamma))/(r*r) )
                    * (rho+press/pow(clight,2))/(1+N_1*(rho + press/pow(clight,2) )*( (dd)/(2*(1+gamma)) - 1  ) );
            ddrhodrr = ddrhodPP*pow(dPdr,2) + drhodP*ddPdrr;
            ddcapFdrr = ddcapFdPP*pow(dPdr,2) + dfRdP*ddPdrr;
            beta = pow(r,2)* ( (ddcapFdrr/capF) - (3.0/2.0)*pow(dcapFdr/capF,2)  );

            dmdr = (1.0/(1+(r*dcapFdr)/(2*fR) )) * (4*pi*pow(r,2)*rho/fR - m/r
                                                           + pow(r,2)/4.0*(R-f/fR)*(pow(clight,2)/ggrav)
                                                           + (pow(r,2)/2.0)*(pow(clight,2)/ggrav)*(1-2*m*ggrav/(r*pow(clight,2)))*( ddcapFdrr/fR - 0.75*pow(dcapFdr/fR,2) + (2*dcapFdr)/(r*fR) )   ) + m/r;






        //////////////////////////////
        // fR - original ENDE
        ////////////////////////////


            if(r>1*pow(10,-6))     // does not matter how I cange this value: always at this value dm_dr gets negative
            {
                dp_dr_lim = dp_dr_lim;
                dm_dr_lim = dmdr;//dm_dr_lim;
                //dp_dr = (m+4.0*pi*r*r*r*press/pow(clight,2))/(r*(r-2.0*ggrav*m/pow(clight,2)));
            }
            else
            {
                dp_dr_lim = (4.0*pi*r*press/pow(clight,2));
                dp_dr_lim = -ggrav*( rho+press/pow(clight,2) )*dp_dr_fR;
                dm_dr_lim = 4*pi*rho*pow(r,2);
            }

            dm_dr = dm_dr_lim;
            dp_dr= dp_dr_lim;

    return make_pair(dm_dr, dp_dr);
}






