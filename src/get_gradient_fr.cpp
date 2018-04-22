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

double p_m = +1;

pair<double, double> get_gradients_fR(double m, double press, double r)
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

    double fR=1+2*a*R/Rp;
    double dRdT = -kappa_2;
    double dTdP = -drhodP + 3.0/pow(clight,2);
    double dRdP = dRdT*dTdP;
    double fRR = 2*a/Rp;
    double ddfRdRR = 0;
    double ddRdTT = 0;
    double dTdrho = -1;

    double dfRdT = fRR*dRdT;
    double dfRdP = dfRdT*dTdP;

    double capF = fR;

    double R_rho = dRdT*dTdrho;
    double f_rho = fR*R_rho;
    double capF_rho = fRR*R_rho;
    double lambda_rho = 0.5*(R_rho - (f_rho/capF_rho));
    double e_B = 1/(1-2*m*ggrav/(r*pow(clight,2)));// - lambda_rho*pow(r,2)*ggrav/(3.0*pow(clight,2)));

    double N_1 = (1.0/fR)*fRR*dRdT*(3/pow(clight,2)-drhodP);
    double N_2 = (1.0/fR)*( ddfRdRR*pow(dRdT,2) + fRR*ddRdTT)*pow(3.0/pow(clight,2)-drhodP,2) - (1.0/fR)*fRR*dRdT*ddrhodPP;

    double aa = r*N_1*( 1/pow(clight,2)-0.75*(rho+press/pow(clight,2) )*N_1 );
    double bb = 2*(1/pow(clight,2)-(rho+press/pow(clight,2) )*N_1 );
    double cc = -(rho+press/pow(clight,2) )*( (1-e_B)/(r) - ((8*pi*r*e_B*press*ggrav/pow(clight,4))/fR) + (r*e_B)/(2*fR)*( fR*R-f ) );

    double dPdr_GR=-ggrav*( rho+press/pow(clight,2) )*(m+4.0*pi*r*r*r*press/pow(clight,2))/(r*(r-2.0*ggrav*m/pow(clight,2)));

    double dPdr;
    if(aa>pow(10,-20)){dPdr = (-bb + p_m*sqrt( pow(bb,2) - 4*aa*cc) )/2*aa; }
    else{dPdr = -cc/bb; }



    /////////
    /// d/dr - Derivatives (now possible since we have dPdr)
    ///

    double dcapFdr = dfRdP*dPdr;
    double dRdr = dRdT*dTdP*dPdr;
    double drhodr = drhodP*dPdr;

    /////////
    /// 2nd Derivatives
    ///
    //double ddTdrr = -ddrhodrr+3*ddPdrr/pow(clight,2);
    double ddTdPP = -ddrhodPP;

    double ddRdPP = ddRdTT*pow(dTdP,2) + dRdT*ddTdPP;
    double ddcapFdPP = ddfRdRR*pow(dRdP,2) + fRR*ddRdPP;

    double alpha = pow(r,2)*( (3.0/4.0)* pow(dcapFdr/capF,2) + (2*dcapFdr)/(r*capF) + (e_B)/(2)*(R-f/capF)    );

    double gamma = r*dcapFdr/(2*capF);
    double A_prime = -(1/(1+gamma))*( (1-e_B)/(r) - e_B*8*pi*r*press*(ggrav/pow(clight,4))/fR   + alpha/r );
    double dd = 0.5*r*A_prime + 2*gamma +1;


    double ddPdrr = (N_2*pow(dPdr,2) - 0.5*A_prime*( (drhodr + dPdr/pow(clight,2))/(rho+press/pow(clight,2)) - dd/r)
            - dd/(2*(1+gamma)) * ( N_2*pow(dPdr,2) + (1-e_B)/(r*r) + e_B/fR*(fR*R-f) + 8*pi*e_B*rho*(ggrav/pow(clight,2))/fR + (gamma*(4-3*gamma))/(r*r) )
            +(e_B)/(2*fR)*( fR*R-f) -  8*pi*e_B*press*(ggrav/pow(clight,4))/fR + (gamma*(2*3*gamma))/(r*r) )
            * (rho+press/pow(clight,2))/(1+N_1*(rho + press/pow(clight,2) )*( (dd)/(2*(1+gamma)) - 1  ) );
    double ddrhodrr = ddrhodPP*pow(dPdr,2) + drhodP*ddPdrr;
    double ddcapFdrr = ddcapFdPP*pow(dPdr,2) + dfRdP*ddPdrr;
    double beta = pow(r,2)* ( (ddcapFdrr/capF) - (3.0/2.0)*pow(dcapFdr/capF,2)  );

    double dmdr = (1.0/(1+(r*dcapFdr)/(2*fR) )) * (4*pi*pow(r,2)*rho/fR - m/r
                                                   + pow(r,2)/4.0*(R-f/fR)*(pow(clight,2)/ggrav)
                                                   + (pow(r,2)/2.0)*(pow(clight,2)/ggrav)*(1-2*m*ggrav/(r*pow(clight,2)))*( ddcapFdrr/fR - 0.75*pow(dcapFdr/fR,2) + (2*dcapFdr)/(r*fR) )   ) + m/r;

    //double dmdr = (1.0)/(1.0+gamma) * ( (4*pi*pow(r,2)*rho)/(capF) + (alpha+beta)*pow(clight,2)/(2*ggrav)  - (m)/(r)*(alpha+beta-gamma)   );




//if(dmdr<0){neg_dmdr++;}
//if(r>1*pow(10,3) && ccount==1 && neg_dmdr==1)
if(r>1*pow(10,3) && ccount==1 && (r_count== 1018 || r_count ==1019))
{
cout << "r_count: " << r_count << endl;
cout << "fR-dmdr= " << dmdr << endl;
cout << "GR-dmdr= " << 4*pi*rho*pow(r,2) << endl;
//cout << "1/(1+gamma) " << (1.0/(1+(r*dcapFdr)/(2*fR) )) << endl;
//cout << "(1-2*m*ggrav/(r*pow(clight,2))) " << (1-2*m*ggrav/(r*pow(clight,2))) << endl;
cout << "1st dmdr " << (1.0/(1+(r*dcapFdr)/(2*fR) )) * (4*pi*pow(r,2)*rho/fR    ) << endl;
cout << "2nd dmdR " <<  (1.0/(1+(r*dcapFdr)/(2*fR) ))*(-m/r) << endl;
cout << "2nd dmdr " << (1.0/(1+(r*dcapFdr)/(2*fR) )) * pow(r,2)/4.0*(R-f/fR)*(pow(clight,2)/ggrav) << endl;
cout << "3rd dmdr " << (1.0/(1+(r*dcapFdr)/(2*fR) )) * (pow(r,2)/2.0)*(pow(clight,2)/ggrav)
        *(1-2*m*ggrav/(r*pow(clight,2)))*
        ( ddcapFdrr/fR) << endl;
cout << "4th dmdr " << (1.0/(1+(r*dcapFdr)/(2*fR) )) * (pow(r,2)/2.0)*(pow(clight,2)/ggrav)
        * (1-2*m*ggrav/(r*pow(clight,2)))*
        (  0.75*pow(dcapFdr/fR,2)) << endl;
cout << "5th dmdr " << (1.0/(1+(r*dcapFdr)/(2*fR) )) * (pow(r,2)/2.0)*(pow(clight,2)/ggrav)
        * (1-2*m*ggrav/(r*pow(clight,2)))*
        (  (2*dcapFdr)/(r*fR) )   << endl;
cout << "6th dmdr " << m/r << endl;
cout << "rho= " << rho << endl;
cout << "press= " << press << endl;
cout << "ddcapFdrr= " << ddcapFdrr << endl;
/*
cout << " 1st dmdr" << (1.0/(1+(r*dcapFdr)/(2*fR) )) * (4*pi*pow(r,2)*rho/fR - m/r
                                                        + pow(r,2)/4.0*(R-f/fR)*(pow(clight,2)/ggrav)
                                                        + (pow(r,2)/2.0)*(pow(clight,2)/ggrav)
                                                        *(1-2*m*ggrav/(r*pow(clight,2)))*
                                                        ( ddcapFdrr/fR - 0.75*pow(dcapFdr/fR,2) + (2*dcapFdr)/(r*fR) )   )  << endl;
*/
cout << "--------" << endl;
}



    double dpdr = dPdr;


    ////////
    /// f(R)
    //double dmdr = 4*pi*rho*pow(r,2);


    //dpdr=dp_dr;
    double dp_dr;
    double dm_dr;

    if(r>1*pow(10,3))     // does not matter how I cange this value: always at this value dm_dr gets negative
        {
            dp_dr = dpdr;
            dm_dr = dmdr;
            //dp_dr = (m+4.0*pi*r*r*r*press/pow(clight,2))/(r*(r-2.0*ggrav*m/pow(clight,2)));
        }
        else
        {
            dp_dr = (4.0*pi*r*press/pow(clight,2));
            dp_dr = -ggrav*( rho+press/pow(clight,2) )*dp_dr;
            dm_dr = 4*pi*rho*pow(r,2);
        }


    return make_pair(dm_dr, dp_dr);
}






