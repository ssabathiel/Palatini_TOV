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

using namespace constants;
using namespace std;

double p_m = +1;

get_gradients_functions get_gradients3[]=
    {
        get_gradients_GR,
        get_gradients_fR,
        get_gradients_fRQ,
        get_gradients_fR_lim
    };

void build_f(vector <double> coeff, double R, double f, double fR, double fRR, double fRRR)
{
    f=0;
    for(int n=1;n<=coeff.size();n++)
    {
        f += coeff[n]*pow(R,n);
        if(n-1>=0) { fR += n*coeff[n]*pow(R,n-1); }
        if(n-2>=0) { fRR += n*(n-1)*coeff[n]*pow(R,n-2); }
        if(n-3>=0) { fRRR += n*(n-1)*(n-2)*coeff[n]*pow(R,n-3);}
    }

}

pair<double, double> get_gradients_fR(double m, double press, double r)
{

    double rho = rho_of_p[num_an](press);
    //double rho_ply = rho_of_p[2](press);



    ////////////////////////////////////////
    //PRECALCULATE TERMS NECESSARY FOR THE GRADIENTS
    /////////////////////////////////////////



    ////////////////////
    //dPdr PREPARATION
    //////////////////////

    bool fR_Olmo = 1;
    double dpdrho = dp_drho(rho,press);
    double drhodP = drho_dp[num_an](press);
    double drhodP_sly = drho_dp[0](press);
    double drhodP_ply = drho_dp[2](press);
    double ddrhodPP = ddrho_dPP[num_an](press);

    rho *= Gdc2;
    press *= Gdc4;
    m = m *= Gdc2;
    dpdrho *= Gdc4/Gdc2;
    drhodP *= Gdc2/Gdc4;
    drhodP_sly *= Gdc2/Gdc4;
    drhodP_ply *= Gdc2/Gdc4;
    ddrhodPP *= Gdc2/pow(Gdc4,2);



    double T=(-rho+3*press/pow(clight,2));
    double R=-kappa_2*T;


    double f=R+a*pow(R,2)/Rp;
    double fR=1+2*a*R/Rp;
    double fRR = 2*a/Rp;
    double ddfRdRR = 0;
    double ddRdTT = 0;
    double ddTdPP = -ddrhodPP;



/*
    double f=R+a*pow(R,4)/Rp;
    double fR=1+4*a*pow(R,3)/Rp;
    double fRR = 12*a*pow(R,2)/Rp;
    double ddfRdRR = 24*a*pow(R,1)/Rp;
*/

/*
    double f=R+a*pow(R,3)/Rp;
    double fR=1+3*a*pow(R,2)/Rp;
    double fRR = 6*a*pow(R,1)/Rp;
    double ddfRdRR = 6*a*pow(R,0)/Rp;
*/

    double dRdT = -kappa_2;
    double dTdP = -drhodP + 3.0/pow(clight,2);
    double dRdP = dRdT*dTdP;
    double ddRdPP = ddRdTT*pow(dTdP,2) + dRdT*ddTdPP;

    double dTdrho = -1+3*dpdrho/pow(clight,2);
    double dfdP = fR*dRdP;
    double dfRdT = fRR*dRdT;
    double dfRdP = dfRdT*dTdP;
    dfRdP = fRR*dRdP;
    double ddfRdPP = ddfRdRR*pow(dRdP,2)  + fRR*ddRdPP;
    double capF = fR;

    double R_rho = dRdT*dTdrho;
    double f_rho = fR*R_rho;
    double capF_rho = fRR*R_rho;
    double lambda_rho = 0.5*(R_rho - (f_rho/capF_rho));
    double e_B = 1/(1-2*m*ggrav/(r*pow(clight,2)));// - lambda_rho*pow(r,2)*ggrav/(3.0*pow(clight,2)));

    double N_1 = (1.0/fR)*fRR*dRdT*(3.0/pow(clight,2) - drhodP);
    double N_2 = (1.0/fR)*( ddfRdRR*pow(dRdT,2) + fRR*ddRdTT)*pow(3.0/pow(clight,2)-drhodP,2) - (1.0/fR)*fRR*dRdT*ddrhodPP;
             //N_2   = ddfRdPP/fR;
    double aa = r*N_1*( 1/pow(clight,2)-0.75*(rho+press/pow(clight,2) )*N_1 );
    double bb = 2*(1/pow(clight,2)-(rho+press/pow(clight,2) )*N_1 );
    double cc = -(rho+press/pow(clight,2) )*( (1-e_B)/(r) - ((8*pi*r*e_B*press*ggrav/pow(clight,4))/fR) + (r*e_B)/(2*fR)*( fR*R-f ) );

    //aa= r*dfRdP/fR * ( 1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*1.5*dfRdP/fR);
    //aa= r*N_1 * ( 1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*1.5*N_1);
    //cc  = 2*(rho+press/pow(clight,2) )/(r*(r-2.0*ggrav*m/pow(clight,2))) * ( m*ggrav/pow(clight,2) - ((f+kappa_2*(press/pow(clight,2) - rho) ) )*pow(r,3)*0.25  );
    //bb = r*2.0/r*(1.0/pow(clight,2) - (rho+press/pow(clight,2) )*0.5*(2*dfRdP/fR)   );
    //cc = 2*(rho+press/pow(clight,2) )/(r*(r-2.0*ggrav*m/pow(clight,2))) * ( m*ggrav/pow(clight,2) - (f+kappa_2*(press/pow(clight,2) - rho) )*pow(r,3)*0.25  );


    double dPdr;
    if(aa>pow(10,-20) && ( pow(bb,2) - 4*aa*cc)>0 )
    {
        dPdr = (-bb + sqrt( pow(bb,2) - 4*aa*cc) )/(2*aa);
    }
    else{dPdr = -cc/bb; }  //cout << "aa=" << aa << endl;

    double dPdr_abc = dPdr;




    ////////////////////////////
    // f(R) OLMO - dPdr
    ///////////////////////////
    //dfRdP = 2*a/Rp*(-kappa_2*(-drhodP+3.0));
    double P_r_0 = (rho + press)/(r*(r-2*m))*(m- ( (f+kappa_2*(press-rho))/(fR)   )*pow(r,3)/4.0   );
    double alpha_r = (rho+press)*dfRdP/fR;
    double beta_r = 2*r*dfRdP/fR * ( 1- 3.0/4.0*(rho + press)*dfRdP/fR  );

    if(fR_Olmo==1){dPdr= -(P_r_0/(1-alpha_r))*2/(1.0+pow(1-beta_r*P_r_0,0.5));}





    ////////////////////////////
    // f(R) OLMO - dPdr  - END
    ///////////////////////////


    /////////
    /// d/dr - Derivatives (now possible since we have dPdr)
    ///

    double dcapFdr = dfRdP*dPdr/pow(clight,2);
    double dRdr = dRdT*dTdP*dPdr/pow(clight,2);
    double drhodr = drhodP*dPdr/pow(clight,2);// /pow(clight,2);

    /////////
    /// 2nd Derivatives
    ///

    double ddcapFdPP = ddfRdRR*pow(dRdP,2) + fRR*ddRdPP;   //= 0 in R^2 theory
    double alpha = pow(r,2)*( (3.0/4.0)* pow(dcapFdr/capF,2) + (2*dcapFdr)/(r*capF) + (e_B)/(2)*(R-f/capF)    );
    double gamma = r*dcapFdr/(2*capF);
    double A_prime = -(1/(1+gamma))*( (1.0-e_B)/r - e_B*8*pi*r*press*(ggrav/pow(clight,4))/fR   + alpha/r );
    //A_prime = -2*dPdr/(rho+press);
    double dd = 0.5*r*A_prime + 2*gamma +1;

    double ddPdrr = (N_2*pow(dPdr/pow(clight,2),2) - 0.5*A_prime*( (drhodr + dPdr/pow(clight,2))/(rho+press/pow(clight,2)) - dd/r)
            - dd/(2*(1+gamma)) * ( N_2*pow(dPdr/pow(clight,2),2) + (1-e_B)/(r*r) + e_B/(2*fR)*(fR*R-f) + 8*pi*e_B*rho*(ggrav/pow(clight,2))/fR + (gamma*(4-3*gamma))/(r*r) )
            +(e_B)/(2*fR)*( fR*R-f) -  8*pi*e_B*press*(ggrav/pow(clight,4))/fR + (gamma*(2-3*gamma))/(r*r) )
            * (rho+press/pow(clight,2))/(1/pow(clight,2)+N_1*(rho + press/pow(clight,2) )*( (dd)/(2*(1+gamma)) - 1  ) );
    //double ddrhodrr = ddrhodPP*pow(dPdr/pow(clight,2),2) + drhodP*ddPdrr;  //extra pow(clight,2)
    double ddcapFdrr = (ddcapFdPP*pow(dPdr/pow(clight,2),2) + dfRdP*ddPdrr);    // quiet 0
    double beta = pow(r,2)* ( (ddcapFdrr/capF) - (3.0/2.0)*pow(dcapFdr/capF,2)  );

    double dmdr = (1.0/(1+(r*dcapFdr)/(2*fR) )) * (4*pi*pow(r,2)*rho/fR*(ggrav/pow(clight,2)) - m*(ggrav/pow(clight,2))/r
                                                   + pow(r,2)/4.0*(R-f/fR)
                                                   + (pow(r,2)/2.0)*(1-2*m*ggrav/(r*pow(clight,2)))*( ddcapFdrr/fR - 0.75*pow(dcapFdr/fR,2) + (2*dcapFdr)/(r*fR) )   ) + m*(ggrav/pow(clight,2))/r;

    //dmdr = (1.0)/(1.0+gamma) * ( (4*pi*pow(r,2)*rho)/(capF)*(ggrav/pow(clight,2)) + (alpha+beta) - (m*(ggrav/pow(clight,2)))/(r)*(alpha+beta-gamma)   );




    ////////////////////////////
    // f(R) OLMO - dMdr
    ///////////////////////////

    double Ar = 1.0/(e_B);
    double dfRdr = dfRdP*dPdr;
    double ll=((dfRdr/fR) + (2.0/r))/r;
    double jj= (f+ kappa_2*(rho + 3*press) )/(2*fR);
    double kk= (dfRdr/fR)*(    ((2*r-3*m)/(r*(r-2*m)))-0.75*(dfRdr/fR)      );


    double dalpha_rdr = (drhodP +1)*dPdr*dfRdP/fR + (rho + press)* (ddfRdPP*dPdr/fR - dfRdP/pow(fR,2)*dfRdP*dPdr);


    double dbeta_rdr = beta_r/r
            + (2*r)*(ddfRdPP*dPdr/fR - dfRdP/pow(fR,2)*dfRdP*dPdr)*(1-0.75*(rho+press)*dfRdP/fR)
            + (2*r)*dfRdP/fR*(-0.75*(drhodP + 1)*dPdr*dfRdP/fR -0.75*(rho+press)*(ddfRdPP*dPdr/fR - dfRdP/pow(fR,2)*dfRdP*dPdr) );


    double ss=1;
    aa= 1+                        0.5*ss* ( (beta_r*P_r_0)/(pow(1-beta_r*P_r_0,0.5)*(1+ss*pow(1-beta_r*P_r_0,0.5)))  );
    bb=(dalpha_rdr/(1-alpha_r)) + 0.5*ss* ( (dbeta_rdr*P_r_0)/(pow(1-beta_r*P_r_0,0.5)*(1+ss*pow(1-beta_r*P_r_0,0.5)))  );

    double big_phi = (f+ kappa_2*(press - rho) )/fR;
    double dPhidP = dfdP + kappa_2*(1-drhodP)/fR - (f+ kappa_2*(press - rho) )/pow(fR,2)*dfRdP;

    cc =  ( (1/pow(clight,2)+drhodP)/(rho+press/pow(clight,2)))*dPdr/pow(clight,2)
            - ((dPhidP*pow(r,3)*0.25)/(m*ggrav/pow(clight,2)-big_phi * pow(r,3)*0.25)  )*dPdr/pow(clight,2)
            - (  (2*(r-m*ggrav/pow(clight,2)))/(r*(r-2*m*ggrav/pow(clight,2)))   + (3*big_phi*r*r*0.25 )/(m*ggrav/pow(clight,2)-big_phi*pow(r,3)*0.25)   );

    dd =  (2/(r-2*m) + 1/(m-big_phi*r*r*r*0.25));


    double new_dmdr = jj + Ar*( (ddfRdPP*dPdr/pow(clight,2)*dPdr/pow(clight,2) + dfRdP*dPdr/pow(clight,2)*(cc*aa+bb))/fR + kk   ) ;

    new_dmdr = new_dmdr/(ll- (Ar*dfRdP*dPdr/pow(clight,2)*dd*aa)/fR);

    if(fR_Olmo==1){dmdr = new_dmdr;}

    ////////////////////////////
    // f(R) OLMO - dMdr  - END
    ///////////////////////////


    double dpdr = dPdr;


    ////////
    /// f(R)

    double dp_dr;
    double dm_dr;
    dp_dr = dpdr;
    dm_dr = dmdr*pow(clight,2)/ggrav;

    dm_dr /= Gdc2;
    m /= Gdc2;
    dp_dr /= Gdc4;




    if(dmdr!=dmdr || dmdr<0)
    {
        /*
        cout << "dpdr= " << dPdr << endl;
        cout << "dMdr= " << dmdr << endl;
        cout <<"rho= " << rho << endl;
        cout << "rho_ply= " << rho_ply << endl;
        */
    }



    return make_pair(dm_dr, dp_dr);
}






