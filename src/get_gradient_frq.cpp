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
    double dpdrho = dp_drho(rho,press);
    double drhodP = drho_dp[num_an](press);//(1.0/dpdrho); //pow(clight,2); //drho_dp[num_an](press);

    double ddrhodPP = ddrho_dPP[num_an](press);


    rho *= Gdc2;
    press *= Gdc4;
    m = m *= Gdc2;
    dpdrho *= Gdc4/Gdc2;
    drhodP *= Gdc2/Gdc4;
    ddrhodPP *= Gdc2/pow(Gdc4,2);

    ////////////////////////////////////////
    //PRECALCULATE TERMS NECESSARY FOR THE GRADIENTS
    /////////////////////////////////////////


    ////////////////////
    //dPdr PREPARATION
    //////////////////////
    double T=(-rho+3*press/pow(clight,2));
    double R=-kappa_2*T;

    double dRdT = -kappa_2;
    double dTdP = -drhodP + 3.0/pow(clight,2);
    double dRdP = dRdT*dTdP;
    double dR_2dP = 2*R*dRdP;
    double ddRdPP = kappa_2*ddrhodPP;
    double ddR_2dPP = 2*R*dRdP;

    double f_tilde = R+a*pow(R,2)/Rp;
    double f_tildeR = 1 + 2*a*R/Rp;
    double df_tildedP = dRdP + 2*R*dRdP/Rp;
    double df_tildeRdP = 2*a*dRdP/Rp;
    double ddf_tildeRdPP = 2*a*ddRdPP/Rp;
    double ddf_tildedPP = ddRdPP + 2*(dRdP*dRdP + R*ddRdPP)/Rp;



    double b=1.0;
    if(Rq/(32.0) * pow(     -  pow( (R)/(Rq) + f_tildeR,2 ) -(4*kappa_2*(rho+press/pow(clight,2) ) )/(Rq), 1  )<0){b=-1.0;}
    double Q=(3*pow(Rq,2)/8)*(1-(2*kappa_2*(rho+press/pow(clight,2))/Rq)+(2*pow(kappa_2,2)*pow((rho-3*press/pow(clight,2)),2)/(3*pow(Rq,2)))-pow(1-(4*kappa_2*(rho+press/pow(clight,2)))/Rq,0.5));     //make Energy-dependent??

    Q = 2*Rq/b* (-(kappa_2*press/pow(clight,2) + f_tilde*0.5 + (Rq/(8.0*b))*pow(f_tildeR,2) )
    + Rq/(32.0*b) * pow(    3*((R*b)/(Rq) + f_tildeR) - pow( pow( (R*b)/(Rq) + f_tildeR,2 ) -(4*b*kappa_2*(rho+press/pow(clight,2) ) )/(Rq), 0.5  ),2) ); //From Helios--> Olmo --> p/m-paper

    double Q_2 = 3*((R*b)/(Rq) + f_tildeR) - pow( pow( (R*b)/(Rq) + f_tildeR,2 ) -(4*b*kappa_2*(rho+press/pow(clight,2) ) )/(Rq), 0.5  );
    double Q_sqr = pow( (R*b)/(Rq) + f_tildeR,2 ) -(4*b*kappa_2*(rho+press/pow(clight,2) ) )/(Rq);
    double dQ_sqrdP = 2*(b*R/Rq + f_tildeR)*(b*dRdP/Rq + df_tildeRdP) - 4*b*kappa_2*(drhodP+1)/Rq;
    double dQ_2dP = 3*(b*dRdP/Rq + df_tildeRdP) - 0.5*pow(Q_sqr,-0.5)*dQ_sqrdP;
    double dQdP = - (kappa_2 + df_tildedP*0.5 + Rq/(8.0*b)*2*f_tildeR*df_tildeRdP  ) + Rq/(32.0*b)*2*Q_2*dQ_2dP;
    dQdP = dQdP*2*Rq/b;


    double f=R+a*pow(R,2)/Rp + Q/Rq;
    double fR=1+2*a*R/Rp;
    double fRR = 2*a/Rp;
    double fQ=1/Rq;



    //double dQdP = (3*pow(Rq,2)/8)*(-2*kappa_2*(drhodP + 1/pow(clight,2))/Rq + (2*pow(kappa_2,2)*2*(rho-3*press/pow(clight,2))*(drhodP-3/pow(clight,2)))/(3*pow(Rq,2))-0.5*pow(1-4*kappa_2*(rho+press/pow(clight,2))/Rq,-0.5)*(-4*kappa_2*(drhodP + 1/pow(clight,2))/Rq));
    //ddQdPP = (3*pow(Rq,2)/8) * ( 36*pow(kappa_2,2)/(3*Rq*Rq) +0.25*pow(1-4*kappa_2*(rho+press/pow(clight,2))/Rq,-3./2.)*pow(-4*kappa_2/Rq,2) );
    double ddQdPP = 3*Rq*Rq*(1./8.)* (  - (2*kappa_2*(ddrhodPP))/(Rq)
                                + 2*kappa_2*2*( (drhodP - 3/pow(clight,2))*(drhodP - 3/pow(clight,2))* (rho-3*press/pow(clight,2))*ddrhodPP   )/(3*Rq*Rq)
                                - 0.5*(  -0.5*pow(1-(4*kappa_2*(rho+press/pow(clight,2) ))/(Rq), -1.5 ) * (-4*kappa_2*(drhodP + 1/pow(clight,2)/Rq))* (-4*kappa_2*(drhodP + 1/pow(clight,2)/Rq))
                                          + (pow(1-(4*kappa_2*(rho+press/pow(clight,2) ))/(Rq), -0.5 ))*(4*kappa_2*(ddrhodPP)/Rq)   )
                            );
    //ddQdPP = 0;

    double ddQsqrdPP = 2*( pow(b*dRdP/Rq + df_tildeRdP,2) + (b*R/Rq + f_tildeR)*(b*ddRdPP/Rq + ddf_tildeRdPP)) - 4*b*kappa_2*ddrhodPP/Rq;
    double ddQ_2dPP = 3*(b*ddRdPP/Rq  + ddf_tildeRdPP) -  0.5*(-0.5*pow(Q_sqr,-1.5)*dQ_sqrdP*dQ_sqrdP + pow(Q_sqr,-0.5)*ddQsqrdPP   );
    ddQdPP = -(0.5*ddf_tildedPP + Rq/(8.0*b)*2*( pow(df_tildeRdP,2) + f_tildeR*ddf_tildeRdPP   ) )  + Rq/(32*b)*2*(dQ_2dP*dQ_2dP + Q_2*ddQ_2dPP);
    ddQdPP = ddQdPP*2*Rq/b;


    double dfdP = fR*dRdP + fQ*dQdP;
    double ddfdPP = ddRdPP*fR + a*ddR_2dPP/Rp + ddQdPP/Rq;
    double dfRdP = fRR*dRdP;
    double ddfRdPP = (2*a/Rp)*(-kappa_2*(-ddrhodPP));
    double dfR_2dP = 2*fR*dfRdP;
    double ddfR_2dPP = 2*(dfRdP*dfRdP + fR*ddfRdPP);

    double lambda=pow(kappa_2*press/pow(clight,2) + f/2.0 + pow(fR,2)/(8*fQ),0.5);
    double lambda_2=pow(lambda,2);


    //if(lambda_2-kappa_2*(rho+press/pow(clight,2))<0){cout << "hey there " << endl;}
    //if(press/rho>a){p_m=-1;}

    double sigma_1=fR*0.5+ p_m*pow(2*fQ,0.5)*pow(lambda_2-kappa_2*(rho+press/pow(clight,2)),0.5);
    double sigma_2=fR*0.5+pow(2*fQ,0.5)*lambda;

    double S=pow(sigma_2,2)/pow(sigma_1*sigma_2,0.5);
    double Omega=pow(sigma_1*sigma_2,0.5);

    double Ar=1-2*m*ggrav/(r*pow(clight,2));

    double det_sigma=sigma_1*pow(sigma_2,3);
    double sqrt_det_sigma=pow(det_sigma,0.5);

    double tau_rr=(1/sqrt_det_sigma)*(f*0.5+kappa_2*press/pow(clight,2))*fR;
    double tau_tt=(1/sqrt_det_sigma)*(f*0.5+kappa_2*(-rho))*fR;

    tau_rr=(1/pow(fR,2))*(f_tilde*0.5+kappa_2*press/pow(clight,2))*fR;
    tau_tt=(1/pow(fR,2))*(f_tilde*0.5+kappa_2*(-rho))*fR;

    double dLambda_2dP = kappa_2/pow(clight,2) + dfdP*0.5 + (1/(8*fQ))*dfR_2dP;
    double dLambdadP = 0.5*pow(lambda_2,-0.5)*dLambda_2dP;

    double dSigma1dP = dfRdP*0.5 + pow(2*fQ,0.5)*pow(lambda_2-kappa_2*(rho+press/pow(clight,2)),-0.5)*0.5*(dLambda_2dP-kappa_2*(drhodP    + 1/pow(clight,2)));
    double dSigma2dP = 0.5*dfRdP+pow(2*fQ,0.5)*dLambdadP;

    //dSigma1dP = dfRdP;
    //dSigma2dP = dfRdP;

    double  dOmegadP = 0.5*pow(sigma_1*sigma_2,-0.5)*sigma_2*dSigma1dP + 0.5*pow(sigma_1*sigma_2,-0.5)*sigma_1*dSigma2dP;
    dOmegadP = pow(sigma_1*sigma_2,-0.5) * 0.5 * (dSigma1dP*sigma_2 + sigma_1* dSigma2dP);

    double dSdSigma1 = -0.5*pow(sigma_2,2)*pow(sigma_1*sigma_2,-3./2.)*sigma_2;
    //double dSdSigma2 = 3*sigma_2*sigma_2/(2*sigma_1)*pow(sigma_2*sigma_2*sigma_2/sigma_1,-0.5);
    double dSdSigma2 = 2*sigma_2*pow(sigma_1*sigma_2,-0.5) + pow(sigma_2,2)*(-0.5)*pow(sigma_1*sigma_2,-1.5)*sigma_1;
    double dSdP = dSdSigma1*dSigma1dP + dSdSigma2*dSigma2dP;
    dSdP = (2*sigma_2*dSigma2dP/(pow(sigma_1*sigma_2,0.5) ) )   - 0.5 * sigma_2*sigma_2*pow(sigma_1*sigma_2,-1.5)* (dSigma1dP*sigma_2 + sigma_1* dSigma2dP);     //skip dSdSigmaversion

    double alpha_r= (rho+press/pow(clight,2))*(dOmegadP/Omega+dSdP/S)/2.0;
    double beta_r=  (2*r)*dOmegadP/Omega*  (1/pow(clight,2)-(rho+press/pow(clight,2))*((3.0/4.0)*(dOmegadP/Omega)-(dOmegadP/Omega-dSdP/S) ));

    double P_r_0 = ggrav*((rho+press/pow(clight,2))/(r*(r-2*ggrav*m/pow(clight,2))))*(m -(tau_rr+(Omega/S)*tau_tt)*pow(r,3)*0.25*pow(clight,2)/ggrav);
    //+4*pi*r*r*r*press);//
    //double P_r_0 = ((rho+press/pow(clight,2))/(r*(r-2*ggrav*m/pow(clight,2))))*(m*ggrav/pow(clight,2)-(tau_rr+(Omega/S)*tau_tt)*pow(r,3)*0.25)*pow(clight,2);

    //P_r_0 =  ((rho+press/pow(clight,2)) )/(r*(r-2*ggrav*m/pow(clight,2)))*(m*ggrav/pow(clight,2) -  (f+kappa_2*(press/pow(clight,2) - rho))/(fR)*pow(r,3)*0.25   )/pow(clight,2);
    //alpha_r = (rho+press/pow(clight,2) ) * dfRdP/fR;
    //beta_r = 2*r*dfRdP/fR * ( 1/pow(clight,2) - (3*(rho+press/pow(clight,2))/4.0)*dfRdP/fR );

    //P_r_0 = (rho + press)/(r*(r-2*m))*(m- ( (f+kappa_2*(press-rho))/(fR)   )*pow(r,3)/4.0   );
    //alpha_r = (rho+press)*dfRdP;
    //beta_r = 2*r*dfRdP/fR * ( 1- 3.0/4.0*(rho + press)*dfRdP/fR  );


    /////////
    //dPdr
    /////////
    double dPdr= -(P_r_0/(1-alpha_r))*2/(1+p_m*pow(1-beta_r*P_r_0,0.5));



    ///////////////
    // c_1 + c_2*pr + c_3*pr2^tactics
    ///////////////

    double c_1  = 2*(rho+press/pow(clight,2) )/(r*(r-2.0*ggrav*m/pow(clight,2))) * ( m*ggrav/pow(clight,2) - ((tau_rr+(Omega/S)*tau_tt) )*pow(r,3)*0.25  );
    double c_2  = r*2.0/r*(1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*(dOmegadP/Omega + dSdP/S)   );
    double c_3  = r*dOmegadP/Omega*( 1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*(1.5*dOmegadP/Omega - (dSdP/S - dOmegadP/Omega)  )  );

    //c_3= dfRdP/fR * ( 1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*1.5*dfRdP/fR);
    //c_1  = 2*(rho+press/pow(clight,2) )/(r*(r-2.0*ggrav*m/pow(clight,2))) * ( m*ggrav/pow(clight,2) - ((f+kappa_2*(press/pow(clight,2) - rho) )/fR )*pow(r,3)*0.25  );
    //c_2 = r*2.0/r*(1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*(2*dfRdP/fR)   );
    /*
    double aa = r*N_1*( 1/pow(clight,2)-0.75*(rho+press/pow(clight,2) )*N_1 );
    double bb = 2*(1/pow(clight,2)-(rho+press/pow(clight,2) )*N_1 );
    double cc = -(rho+press/pow(clight,2) )*( (1-e_B)/(r) - ((8*pi*r*e_B*press*ggrav/pow(clight,4))/fR) + (r*e_B)/(2*fR)*( fR*R-f ) );
    */
    //c_3 = r*dOmegadP/Omega*( 1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*(1.5*dfRdP/fR  )  );
    // c
    // b
    // a
    //if(c_3<pow(10,-10) ) {cout << "c3 = a =" << c_3 << endl;}
/*
    cout << "c_1= " << c_1 << endl;
    cout << "c_2= " << c_2 << endl;
    cout << "c_3= " << c_3 << endl;
    */

    //if(pow(c_2,2) - 4*c_3*c_1<0 ) {cout << "pow(c_2,2) - 4*c_3*c_1) = " << pow(c_2,2) - 4*c_3*c_1 << endl;}

    if(abs(c_3)>pow(10,-20)){dPdr = dPdr;}// (-c_2 + p_m*sqrt( pow(c_2,2) - 4*c_3*c_1) )/(2*c_3); } //
    else{
        //dPdr = -c_1/c_2;
    }

    if(dPdr!=dPdr && ccount==5)
    {
        //cout << "ccount= " << endl;
        cout << "r_count= " << r_count << endl;
        cout << "dPdr= " << dPdr << endl;
        cout << "P_r_0= " << P_r_0 << endl;
        cout << "alpha_r= " << P_r_0 << endl;
        cout << "beta_r= " << beta_r<< endl;
        cout << "Omega= " << Omega << endl;
        cout << "S= " << S << endl;
        cout << "tau_rr= " << tau_rr << endl;
        cout << "tau_tt= " << tau_tt << endl;
        cout << "sigma_1= " << sigma_1 << endl;
        cout << "sigma_2= " << sigma_1 << endl;
        cout << "dQdP= " << dQdP << endl;
        cout << "1-beta_r*P_r_0= " << 1-beta_r*P_r_0 << endl;



        cout << endl;
    }




    ///////////////////
    //dMdr PREPARATIoN
    ///////////////////
    double dOmegadr = dOmegadP*dPdr/pow(clight,2);
    double dSdr = dSdP*dPdr/pow(clight,2);


        //New
        double ddSigma1drdP;
        double ddSigma2drdP;

        double ddLambda_2dPP = ddfdPP*0.5 + (1/(8*fQ))*ddfR_2dPP;
        double ddLambda_PP = -0.25*pow(lambda_2,-3./2.)*dLambda_2dP + 0.5*pow(lambda_2,-0.5)*ddLambda_2dPP;

        ddSigma1drdP = pow(2*fQ,0.5)*(-0.25)*pow(lambda_2-kappa_2*(rho+press/pow(clight,2)),-1.5)*(dLambda_2dP*dPdr/pow(clight,2)-kappa_2*dPdr/pow(clight,2)) +
                        pow(2*fQ,0.5)*pow(lambda_2-kappa_2*(rho+press/pow(clight,2)),-0.5) * ddLambda_2dPP*dPdr/pow(clight,2);

        ddSigma2drdP = pow(2*fQ,0.5)*0.5*(-0.5*pow(lambda_2,-1.5)*dLambda_2dP*dPdr/pow(clight,2)*dLambda_2dP  + 0.5*pow(lambda_2,-0.5)*(pow(kappa_2,2)/pow(clight,2))*dPdr/pow(clight,2)    );



        double ddOmegadrdP = -0.25*pow(sigma_1*sigma_2,-1.5)*(sigma_1*dSigma2dP*dPdr/pow(clight,2)+ sigma_2*dSigma1dP*dPdr/pow(clight,2)) * (sigma_1*dSigma2dP + sigma_2*dSigma1dP)
                    + 0.5*pow(sigma_1*sigma_2,-0.5)*( dSigma1dP*dPdr/pow(clight,2)*dSigma2dP + sigma_1*ddSigma2drdP + dSigma2dP*dPdr/pow(clight,2)*dSigma1dP + sigma_2*ddSigma1drdP  );



        double hh = 0.5*pow(lambda_2 - kappa_2*(rho + press/pow(clight,2)),-0.5);
        double gg = dLambda_2dP - kappa_2*(drhodP + 1/pow(clight,2));
        double dhhdP = -0.5*pow(lambda_2 - kappa_2*(rho + press/pow(clight,2)),-1.5)* (dLambda_2dP - kappa_2*(drhodP + 1/pow(clight,2)) );   //pow(hh,-1)
        double dggdP = ddLambda_2dPP - kappa_2*ddrhodPP;
        double ddSigma1dPP = ddfRdPP*0.5 + pow(2*fQ,0.5)*(dhhdP*gg+hh*dggdP); //dhhdP*gg+
        double ddSigma2dPP = ddfRdPP*0.5 + pow(2*fQ,0.5)*ddLambda_PP;

        double dOmegadSigma1 = 0.5*pow(sigma_1*sigma_2,-0.5)*sigma_2;
        double dOmegadSigma2 = 0.5*pow(sigma_1*sigma_2,-0.5)*sigma_1;
        double ddOmegadSigma1_2 = -0.25*pow(sigma_1*sigma_2, 3./2.)*sigma_2*sigma_2;
        double ddOmegadSigma2_2 = -0.25*pow(sigma_1*sigma_2, 3./2.)*sigma_1*sigma_1;
/*
        dSigma1dP = dfRdP;
        dSigma2dP = dfRdP;
        sigma_1 = fR;
        sigma_2 = fR;

*/
        //ddSigma1dPP = ddfRdPP;
        //ddSigma2dPP = ddfRdPP;

        double ddOmegadPP = -0.25*pow(sigma_1*sigma_2,-1.5)*pow( sigma_1*dSigma2dP + sigma_2*dSigma1dP,2)
                            +  0.5*pow(sigma_1*sigma_2,-0.5)*(sigma_1*ddSigma2dPP + 2*dSigma1dP*dSigma2dP + sigma_2*ddSigma1dPP);



        //dOmegadP = pow(sigma_1*sigma_2,-0.5) * 0.5 * (dSigma1dP*sigma_2 + sigma_1* dSigma2dP);


//________________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________________


//FROM HERE TO NEXT DOUBLE LINE WAS ATTEMPT TO DERIVE omega_rr and P_rr myself: but did mistakes by not including all derivatives: e.g. of rho

        double ddSdSigma1dSigma2 = 3*sigma_2*sigma_2*(-0.5)*pow(sigma_1*sigma_2,-3./2.) + pow(sigma_2,3)*(-0.5)*(-1.5)*pow(sigma_1*sigma_2,-5./2)*sigma_1;
        double ddSdSigma1Sigma1 = pow(sigma_2,4)*(-0.5)*(-1.5)*pow(sigma_1*sigma_2,-5./2);
        double ddSdSigma2Sigma2 = -0.25*pow(sigma_2*sigma_2*sigma_2/sigma_1,-1.5)*3*sigma_2*sigma_2/sigma_1 + 0.5*pow(sigma_2*sigma_2*sigma_2/sigma_1,-0.5)*6*sigma_2/sigma_1;


        //S_PP
        double S_1 = pow(sigma_2,2);
        double S_2 = pow(sigma_1*sigma_2,-0.5);
        double dS_1dP = 2*sigma_2*dSigma2dP;
        double dS_2dP1 = -0.5*pow(sigma_1*sigma_2,-1.5);
        double dS_2dP2 = dSigma1dP*sigma_2 + dSigma2dP*sigma_1;
        double dS_2dP = dS_2dP1*dS_2dP2;

        double dS_21dP = 0.75*pow(sigma_1*sigma_2,-2.5)*dS_2dP2;
        double dS_22dP = ddSigma1dPP*sigma_2 + 2*dSigma1dP*dSigma2dP + ddSigma2dPP*sigma_1;
        double ddS_1dPP = 2* ( dSigma2dP*dSigma2dP + sigma_2*ddSigma2dPP );
        double ddS_2dPP = dS_21dP*dS_2dP2 + dS_2dP1*dS_22dP;
       //double ddSdPP = dSdSigma1*ddSigma1dPP + dSdSigma2*ddSigma2dPP + ddSdSigma1Sigma1*pow(dSigma1dP,2) + 2*ddSdSigma1dSigma2*dSigma1dP*dSigma2dP + ddSdSigma2Sigma2*pow(dSigma2dP,2);
        double ddSdPP = ddS_1dPP*S_2 + 2*dS_1dP*dS_2dP + ddS_2dPP*S_1;
        //ddSdPP = ddfRdPP;
/*
        cout << "=============" << endl;
        cout << "ddOmegadPP = " << ddOmegadPP << endl;
        cout << "ddSdPP_1= " << ddS_1dPP*S_2 << endl;
        cout << "ddSdPP_2= " << dS_1dP*dS_2dP << endl;
        cout << "ddSdPP_3= " << ddS_2dPP*S_1 << endl;
        */

        /*
        //ddSdPP directly from Mathematica
        ddSdPP = (2*pow(dSigma2dP,2) / ( pow(sigma_1*sigma_2, 0.5) )   )
                - (2*sigma_2*dSigma2dP*(sigma_2*dSigma1dP + sigma_1*dSigma2dP) ) / (pow(sigma_1*sigma_2,1.5))
                + 3*sigma_2*sigma_2 * pow(sigma_2*dSigma1dP - sigma_1*dSigma2dP, 2 ) / (pow(4*sigma_1*sigma_2,2.5))
                + 2*sigma_2*ddSigma2dPP / (pow(sigma_1*sigma_2,0.5))
                - ( sigma_2*sigma_2* ( 2*dSigma1dP*dSigma2dP + sigma_2*ddSigma1dPP + sigma_1*ddSigma2dPP )  ) / (2*pow(sigma_1*sigma_2,1.5));
*/



        //ALPHA(r) - BETA(r)
        //often needed derivative of Omega_P/Omega and of S_P/S
        double dSPDSdP = dSdP            *((-1)/pow(S,2))    *dSdP + ddSdPP/S;
        double dOmegaPDOmegadP = dOmegadP*((-1)/pow(Omega,2))*dOmegadP + ddOmegadPP/Omega;
        double dalpha_rdr = 0.5*(drhodP+ 1/pow(clight,2) )*(dOmegadP/Omega + dSdP/S) + (rho + press/pow(clight,2))*0.5*(dOmegaPDOmegadP + dSPDSdP);
        dalpha_rdr = dalpha_rdr*dPdr/pow(clight,2);
        double beta_r1 = 2*r/dPdr*dOmegadP/Omega;
        double beta_r2 = 1- (rho + press/pow(clight,2))*0.5*( 1.5* dOmegadP/Omega - (dOmegadP/Omega - dSdP));
        double dbeta_r1dP = 2.0/pow(clight,2)*dOmegadP/Omega + 2*r*dOmegaPDOmegadP;
        double dbeta_r2dP = - (drhodP + 1.0/pow(clight,2) )/2.0*( 1.5* dOmegadP/Omega - (dOmegadP/Omega - dSdP))   - (rho + press/pow(clight,2))*0.5* ( 1.5* dOmegaPDOmegadP - (dOmegaPDOmegadP - dSPDSdP)) ;

        double dbeta_rdP = dbeta_r1dP*beta_r2 + dbeta_r2dP*beta_r1;
        double dbeta_rdr = dbeta_rdP*dPdr/pow(clight,2);







        dbeta_rdr = beta_r/r
                    + (2*r)*(ddOmegadPP*dPdr/Omega - dOmegadP/pow(fR,2)*dOmegadP*dPdr)*(1- (rho + press/pow(clight,2))*0.5*( 1.5* dOmegadP/Omega - (dOmegadP/Omega - dSdP/S)) )
                    + (2*r)*dOmegadP/Omega*(-0.75*(drhodP + 1)*dPdr*( dOmegadP/Omega - (dOmegadP/Omega - dSdP/S)/1.5) - 0.75*(rho+press)*(  (ddOmegadPP*dPdr/fR - dfRdP/pow(fR,2)*dfRdP*dPdr) - (dOmegaPDOmegadP - dSPDSdP)/1.5)*dPdr );
/*
        ddOmegadPP = ddfRdPP;
        dalpha_rdr = (drhodP +1)*dPdr*dfRdP/fR + (rho + press)* (ddfRdPP*dPdr/fR - dfRdP/pow(fR,2)*dfRdP*dPdr);
        dbeta_rdr = beta_r/r
                    + (2*r)*(ddfRdPP*dPdr/fR - dfRdP/pow(fR,2)*dfRdP*dPdr)*(1-0.75*(rho+press)*dfRdP/fR)
                    + (2*r)*dfRdP/fR*      (-0.75*(drhodP + 1)*dPdr *dfRdP/fR                                        -0.75*(rho+press)*(ddfRdPP*dPdr/fR - dfRdP/pow(fR,2)*dfRdP*dPdr) );
*/

        double dSPdr = ddSdPP*dPdr/pow(clight,2);
        double dOmegaPdr = ddOmegadPP*dPdr/pow(clight,2);
        double dSPDSdr = dSPdr/S + dSdP*(-1)*pow(S,-2)*dSdr;
        double dOmegaPDOmegadr = dOmegaPdr/Omega + dOmegadP*(-1)*pow(Omega,-2)*dOmegadr;
        double v2 = 1.5*dOmegadP/Omega - (dOmegadP/Omega + dSdP/S);
        double dv2dr = 1.5*dOmegaPDOmegadr - (dOmegaPDOmegadr + dSPDSdr);
        double v = 1- 0.5*(rho + press/pow(clight,2))*v2;
        double dvdr = - 0.5*(drhodP+ dPdr/pow(clight,2)  )*dv2dr - 0.5*(rho + press/pow(clight,2))*v2;


        double d1Dsqrt_det_sigmadP = -0.5*pow(det_sigma,-3./2.)*(pow(sigma_2,3)*dSigma1dP + 3*sigma_1*pow(sigma_2,2)*dSigma2dP);
        //d1Dsqrt_det_sigmadP = -1.0/pow(fR,2)*dfRdP;
        double dtau_rrdP = d1Dsqrt_det_sigmadP*(f*0.5+kappa_2*press/pow(clight,2))+ 1/sqrt_det_sigma*(dfdP*0.5+kappa_2/pow(clight,2));
        double dtau_ttdP = d1Dsqrt_det_sigmadP*(f*0.5-kappa_2*rho)+ 1/sqrt_det_sigma*(dfdP*0.5 - kappa_2*drhodP);



        double P_rr_0 = dPdr/pow(clight,2)*(r*(r-2*m*ggrav/pow(clight,2)))*(m*ggrav/pow(clight,2)-(tau_rr+Omega*(1/S)*tau_tt)*pow(r,3)*0.25)
                       +(rho+press/pow(clight,2))*(2*r-2*m*ggrav/pow(clight,2))*(m*ggrav/pow(clight,2)-(tau_rr+Omega*(1/S)*tau_tt)*pow(r,3)*0.25)
                       +(rho+press/pow(clight,2))*(r*(r-2*m*ggrav/pow(clight,2)))*(  (3*r*r/4.)*m*ggrav/pow(clight,2) - dtau_rrdP*dPdr*pow(r,3)/4 - tau_rr*3+pow(r,2)/4
                                                    - dOmegadP*dPdr/pow(clight,2)*(1/S)*tau_tt*pow(r,3)*0.25
                                                    - Omega*(-1)*pow(S,-2)*dSdP*dPdr/pow(clight,2)*tau_tt*pow(r,3)*0.25
                                                    - Omega*(1/S)*dtau_ttdP*dPdr/pow(clight,2)*pow(r,3)*0.25
                                                    - Omega*(1/S)*tau_tt*3*pow(r,2)*0.25);

        double P_rr= -P_rr_0*(1/(1-alpha_r))*2*1/(1+pow(1-beta_r*P_r_0,0.5))
                     -P_r_0*(-1)*pow(1-alpha_r,-2)*dalpha_rdr*2*1/(1+pow(1-beta_r*P_r_0,0.5))
                     -P_rr_0*(1/(1-alpha_r))*2*(-1)*pow(1+pow(1-beta_r*P_r_0,0.5),-2)
                                                   *0.5*pow(1-beta_r*P_r_0,-0.5)
                                                   *(-1)*(dbeta_rdr*P_r_0 + beta_r*P_rr_0);

        double ddOmegadrr = ddOmegadPP*dPdr/pow(clight,2)*dPdr/pow(clight,2) + dOmegadP*P_rr;


//________________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________________

    //P_r_0 = P_r_0/pow(clight,2);

    //Preparing l,j,k ,  a,b,c,d
    double ll=((dOmegadr/Omega) + (2.0/r))/r;
    double jj=0.5*(3*tau_rr - (Omega/S)*tau_tt);

    double dfRdr = dfRdP*dPdr;
    //jj = (f+ kappa_2*(rho + 3*press) )/(2*fR);
    //ll = ((dfRdr/fR) + (2.0/r))/r;
    double kk= (dOmegadr/Omega)*(    ((2*r-3*m*ggrav/pow(clight,2))/(r*(r-2*m*ggrav/pow(clight,2))))-0.75*(dOmegadr/Omega)      );

    //dalpha_rdr = 0.0;
    //alpha_r = 0.0;
    //dbeta_rdr = 0.0;
    //beta_r = 0.0;

    //The following variables are used in Helios paper to calculate Prr. But i calculated Prr manually above, but wrong: did not take derivative with respect to r all variables e.g. rho
    double ss=1;
    double aa= 1+                        0.5*ss* ( (beta_r*P_r_0)/(pow(1-beta_r*P_r_0,0.5)*(1+ss*pow(1-beta_r*P_r_0,0.5)))  );
    double bb=(dalpha_rdr/(1-alpha_r)) + 0.5*ss* ( (dbeta_rdr*P_r_0)/(pow(1-beta_r*P_r_0,0.5)*(1+ss*pow(1-beta_r*P_r_0,0.5)))  );


    double big_phi = (tau_rr + (Omega/S)*tau_tt);
    //double drhodp = drho_dp[num_an](press);
    double dPhidP = dtau_rrdP + dtau_ttdP*(Omega/S) + tau_tt*(dOmegadP/S - Omega/pow(S,2)*dSdP);
    /*
    double cc = ( (1+drhodp)/(rho+press/pow(clight,2)) - (dPhidP*pow(r,3)*0.25)/(m-big_phi * pow(r,3)*0.25)  )*dPdr -
    (  (2*(r-m*ggrav/pow(clight,2)))/(r*(r-2*m*ggrav/pow(clight,2))) + (3*big_phi*r*r*0.25 )/(m-big_phi*pow(r,3)*0.25)   );
    */
    double cc =  ( (1/pow(clight,2)+drhodP)/(rho+press/pow(clight,2)))*dPdr/pow(clight,2)
            - ((dPhidP*pow(r,3)*0.25)/(m*ggrav/pow(clight,2)-big_phi * pow(r,3)*0.25)  )*dPdr/pow(clight,2)
            - (  (2*(r-m*ggrav/pow(clight,2)))/(r*(r-2*m*ggrav/pow(clight,2)))   + (3*big_phi*r*r*0.25 )/(m*ggrav/pow(clight,2)-big_phi*pow(r,3)*0.25)   );

//
     //aa=0;
     double dd =  (2/(r-2*m*ggrav/pow(clight,2)) + 1/(m*ggrav/pow(clight,2)-big_phi*r*r*r*0.25));





     //aa= 1+                        0.5*ss* ( (beta_r*P_r_0)/(pow(1-beta_r*P_r_0,0.5)*(1+ss*pow(1-beta_r*P_r_0,0.5)))  );
     //bb=(dalpha_rdr/(1-alpha_r)) + 0.5*ss* ( (dbeta_rdr*P_r_0)/(pow(1-beta_r*P_r_0,0.5)*(1+ss*pow(1-beta_r*P_r_0,0.5)))  );

     big_phi = (f+ kappa_2*(press - rho) )/fR;
     dPhidP = dfdP + kappa_2*(1-drhodP)/fR - (f+ kappa_2*(press - rho) )/pow(fR,2)*dfRdP;
/*
     cc =  ( (1/pow(clight,2)+drhodP)/(rho+press/pow(clight,2)))*dPdr/pow(clight,2)
             - ((dPhidP*pow(r,3)*0.25)/(m*ggrav/pow(clight,2)-big_phi * pow(r,3)*0.25)  )*dPdr/pow(clight,2)
             - (  (2*(r-m*ggrav/pow(clight,2)))/(r*(r-2*m*ggrav/pow(clight,2)))   + (3*big_phi*r*r*0.25 )/(m*ggrav/pow(clight,2)-big_phi*pow(r,3)*0.25)   );
*/



    ///////////////////
    //dMdr
    ///////////////////

    //New dMdr
    double new_dmdr = jj + Ar*( (ddOmegadPP*dPdr/pow(clight,2)*dPdr/pow(clight,2) + dOmegadP*dPdr/pow(clight,2)*(cc*aa+bb))/Omega + kk   ) ;
    new_dmdr = new_dmdr/(ll- (Ar*dOmegadP*dPdr/pow(clight,2)*dd*aa)/Omega);

    //if(c_3>pow(10,-30)){new_dmdr = 4*pi*rho*pow(r,2);}

/*
    //////////////////////
    /// dMdr with the role Model of f(R) equations: see Reexamination of polytropic spheres in Palatini f(R) gravity, Olmo
    /////

    new_dmdr = jj*pow(r,2) + (r*(r-2*ggrav*m/pow(clight,2)))*0.5
            * ( (ddOmegadPP/Omega - 0.75*pow(dOmegadP/Omega,2))*dPdr*dPdr + dOmegadP/Omega*P_rr  +  ((r-m*ggrav/pow(clight,2))/(r*(r-2*m*ggrav/pow(clight,2))))*dOmegadP/Omega*dPdr    );

    new_dmdr = new_dmdr/(1+r*0.5*(dOmegadP/Omega)*dPdr  );
*/
    //double dm_dr = 4*pi*rho*pow(r,2);

    double dm_dr = new_dmdr*pow(clight,2)/ggrav;
    double dp_dr = dPdr;

    dm_dr = get_gradients2[1](m/Gdc2, press/Gdc4, r).first*Gdc2;
    //dp_dr = get_gradients2[1](m /Gdc2, press /Gdc4, r).second*Gdc4;


    //double dm_dr_fR = get_gradients2[1](m, press, r).first;
    //double dp_dr_fR = get_gradients2[1](m, press, r).second;


    //if((ccount==21 && dm_dr_fR/dm_dr > 1.5 && r>1*pow(10,-6)) || (ccount==15 &&r_count%100==0 && r>1*pow(10,-6)) )
    if(ccount==21 && ccount == 22 && (r_count%100==0 || Ar*( (dOmegadP*dPdr/pow(clight,2)*(cc*aa+bb))/Omega    )/jj>0.1) )
    {
        cout << endl;
        cout << "===========" << endl;
        cout << "f(R,Q)" << endl;
        cout << "r_count: " << r_count << endl;
        cout << "r: " << r << endl;
        //cout << "dmdr fR:  " << dm_dr_fR << endl;
        cout << "dmdr fRQ: " << dm_dr << endl;
        cout << "dmdr fRQ: " << (jj + Ar*( (ddOmegadPP*dPdr/pow(clight,2)*dPdr/pow(clight,2) + dOmegadP*dPdr/pow(clight,2)*(cc*aa+bb))/Omega + kk   ) )/(ll- (Ar*dOmegadP*dPdr/pow(clight,2)*dd*aa)/Omega) *pow(clight,2)/ggrav << endl;
        cout << "1st dmdr: " << ( jj )                                                                                                                 /(ll- (Ar*dOmegadP*dPdr/pow(clight,2)*dd*aa)/Omega)*pow(clight,2)/ggrav<< endl;
        cout << "2nd dmdr:: " << (Ar*( (ddOmegadPP*dPdr*dPdr/pow(clight,4))/Omega    ))                                                                /(ll- (Ar*dOmegadP*dPdr/pow(clight,2)*dd*aa)/Omega)*pow(clight,2)/ggrav<< endl;
        cout << "3rd dmdr:: " << (Ar*( (dOmegadP*dPdr/pow(clight,2)*(cc*aa+bb))/Omega    ))                                                            /(ll- (Ar*dOmegadP*dPdr/pow(clight,2)*dd*aa)/Omega)*pow(clight,2)/ggrav<< endl;
        cout << "4th dmdr: " << (  Ar*kk )                                                                                                             /(ll- (Ar*dOmegadP*dPdr/pow(clight,2)*dd*aa)/Omega)*pow(clight,2)/ggrav<< endl;

    }


    dm_dr /= Gdc2;
    m /= Gdc2;
    dp_dr /= Gdc4;



    return make_pair(dm_dr, dp_dr);
}
