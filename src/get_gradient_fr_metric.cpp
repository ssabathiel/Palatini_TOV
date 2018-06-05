
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
#include "header/ply_analytical.h"
#include "header/meta_functions.h"
#include "header/get_gradient_fr.h"

#include "header/get_gradient_frq.h"
#include "header/get_gradient_gr.h"
#include "header/get_gradient_fr_lim.h"
#include "header/get_gradient_fr_metric.h"

using namespace constants;
using namespace std;



get_gradients_functions get_gradients4[]=
    {
        get_gradients_GR,
        get_gradients_fR,
        get_gradients_fRQ,
        get_gradients_fR_lim
    };


// Links
// https://arxiv.org/pdf/1307.7977.pdf  // That's where you get the equations from.
//https://arxiv.org/pdf/0910.5480.pdf
// https://arxiv.org/pdf/1003.3179.pdf  //Good one: metric f(R) for many EOSs
// https://arxiv.org/pdf/1301.5189.pdf  // Also good to compare
// https://arxiv.org/pdf/1112.4154.pdf  // f(R,Q)
// https://arxiv.org/pdf/1408.3856.pdf
// https://ac.els-cdn.com/S0370269315000404/1-s2.0-S0370269315000404-main.pdf?_tid=de192040-c4b8-4fee-af3a-bbb59ea00d91&acdnat=1527331338_49adf0649a4a7ef587a1f5c71de2db73
//
pair<double, double> get_gradients_fR_metric(double m, double press, double r)
{


    double rho = rho_of_p[num_an](press);
    ////////////////////////////////////////
    //PRECALCULATE TERMS NECESSARY FOR THE GRADIENTS
    /////////////////////////////////////////



    ////////////////////
    //dPdr PREPARATION
    //////////////////////

    bool fR_Olmo = 1;
    double dpdrho = dp_drho(rho,press);
    double drhodP = drho_dp[num_an](press);

    double ddrhodPP = ddrho_dPP[num_an](press);


    rho *= Gdc2;
    press *= Gdc4;
    m = m *= Gdc2;
    dpdrho *= Gdc4/Gdc2;
    drhodP *= Gdc2/Gdc4;
    ddrhodPP *= Gdc2/pow(Gdc4,2);
//cout << "ddrhodPP= " << ddrhodPP << endl;


    // dPdr
    double T=(-rho+3*press/pow(clight,2));
    double R=-kappa_2*T;
    double h = pow(R,2);
    double hR = 2*R;
    double hRR = 2.0;
    double ddhRdRR = 0.0;
    double dTdP = -drhodP + 3;
    double dRdT = -kappa_2;
    double dRdP = dRdT*dTdP;
    double alpha = 1.0/Rp;
    double r_g = 147473.0;
    double Ar = pow(1-2*m/r,-1);

    double aa = -Ar/pow(r,2)*(rho+press);
    double bb = m + 4*pi*pow(r,3)*press;
    double cc = 0.25*(hR*R-h);
    double dd = 4*pi*press*hR;


    double C1 = -aa*(bb-alpha*pow(r,3)*(cc+dd));
    double C2 = 1+aa*alpha*pow(r,3)*1.0/(Ar*r)*hRR*dRdP;
    double C3 = aa*alpha*pow(r,3)*1.0/(2*Ar*(rho+press))*hRR*dRdP;


    double dPdr;
    if(abs(C3)>pow(10,-20)){ dPdr = (-C2 + (sqrt( pow(C2,2) - 4*C3*C1) ))/(2*C3); } //
    else                   { dPdr = -C1/C2;                                    }



    //Prepare 2nd derivatives for Prr for dMdr (taken from )
    double f=R+a*pow(R,2)/Rp;
    double fR=1+2*a*R/Rp;
    double fRR = 2*a/Rp;
    double ddfRdRR = 0;
    double ddRdTT = 0;
    double ddTdPP = -ddrhodPP;
    //double dRdP = dRdT*dTdP;
    double ddRdPP = ddRdTT*pow(dTdP,2) + dRdT*ddTdPP;

    double dfdP = fR*dRdP;
    double dfRdT = fRR*dRdT;
    double dfRdP = dfRdT*dTdP;
    dfRdP = fRR*dRdP;
    double ddfRdPP = ddfRdRR*pow(dRdP,2)  + fRR*ddRdPP;
    double capF = fR;

    double N_1 = (1.0/fR)*fRR*dRdT*(3.0/pow(clight,2) - drhodP);
    double N_2 = (1.0/fR)*( ddfRdRR*pow(dRdT,2) + fRR*ddRdTT)*pow(3.0/pow(clight,2)-drhodP,2) - (1.0/fR)*fRR*dRdT*ddrhodPP;

    double dcapFdr = dfRdP*dPdr/pow(clight,2);
    double dRdr = dRdT*dTdP*dPdr/pow(clight,2);
    double drhodr = drhodP*dPdr/pow(clight,2);// /pow(clight,2);

    /////////
    /// 2nd Derivatives
    ///

    double ddcapFdPP = ddfRdRR*pow(dRdP,2) + fRR*ddRdPP;   //= 0 in R^2 theory
    double e_B = 1.0/(1-2*m*ggrav/(r*pow(clight,2)));
    double alpha2 = pow(r,2)*( (3.0/4.0)* pow(dcapFdr/capF,2) + (2*dcapFdr)/(r*capF) + (e_B)/(2.0)*(R-f/capF)    );
    double gamma = r*dcapFdr/(2*capF);
    double A_prime = -(1.0/(1+gamma))*( (1.0-e_B)/r - e_B*8*pi*r*press*(ggrav/pow(clight,4))/fR   + alpha2/r );
    //A_prime = -2*dPdr/(rho+press);
    double dd2 = 0.5*r*A_prime + 2*gamma +1;

    double ddPdrr = (N_2*pow(dPdr/pow(clight,2),2) - 0.5*A_prime*( (drhodr + dPdr/pow(clight,2))/(rho+press/pow(clight,2)) - dd2/r)
            - dd2/(2*(1+gamma)) * ( N_2*pow(dPdr/pow(clight,2),2) + (1-e_B)/(r*r) + e_B/(2*fR)*(fR*R-f) + 8*pi*e_B*rho*(ggrav/pow(clight,2))/fR + (gamma*(4-3*gamma))/(r*r) )
            +(e_B)/(2*fR)*( fR*R-f) -  8*pi*e_B*press*(ggrav/pow(clight,4))/fR + (gamma*(2-3*gamma))/(r*r) )
            * (rho+press/pow(clight,2))/(1/pow(clight,2)+N_1*(rho + press/pow(clight,2) )*( (dd2)/(2*(1+gamma)) - 1  ) );

    double dhRdP = hRR*dRdP;
    double dhRdr = dhRdP*dPdr;

    double ddhRdPP = ddhRdRR*dRdP + hRR*ddRdPP;
    double ddhRdrr = ddhRdPP*dPdr*dPdr + dhRdP*ddPdrr;

    double dmdr = 4*pi*rho*pow(r,2) - alpha*pow(r,2)*( 4*pi*rho*hR - 0.25*(hR*R - h) -  1.0/(2*Ar)*( (2.0/r - A_prime/(2*Ar) )*dhRdr +  ddhRdrr        )  );





    ////////
    /// f(R)

    double dp_dr;
    double dm_dr;
    dp_dr = dPdr;//get_gradients4[1](m/Gdc2, press/Gdc4, r).second*Gdc4; //dPdr;
    dm_dr = dmdr;//get_gradients4[1](m/Gdc2, press/Gdc4, r).first*Gdc2;   //dmdr;

    dm_dr /= Gdc2;
    m /= Gdc2;
    dp_dr /= Gdc4;

    return make_pair(dm_dr, dp_dr);
}






