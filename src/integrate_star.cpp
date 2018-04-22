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
#include "header/glob_variables.h"
#include "header/eos_functions.h"
#include "header/eos_analytical.h"
#include "header/integrate_star.h"

#include "header/meta_functions.h"
#include "header/get_gradient_fr.h"
#include "header/get_gradient_gr.h"
#include "header/get_gradient_frq.h"

using namespace constants;
using namespace std;

bool dr_adapted = 1;
bool single_star=0;
bool euler_method = 0;
int runge_step=0;


get_gradients_functions get_gradients[]=
    {
        get_gradients_GR,
        get_gradients_fR,
        get_gradients_fRQ,
    };




///////////////////////////////////////////////
//INTEGRATE MASS AND RADIUS, from r=0 to R(p=0)
///////////////////////////////////////////////

//important, right order of initialization and first steps of parameters

pair<double, double> tov_integrate(double rho_C)
{
    double press = p_of_rho_num(rho_C);

    double press_now = press;

    double m_now = 0.0;
    double r=0.;
    int count = 0;
    double rho = rho_C;
    double dm_dr;
    double dp_dr;
    r_count = 0;

    double dr = 1.0*pow(10.0,(-20));      //=0.001km = 1m = 100cm


    while(press_now > min_press && count < max_r_iterations)
    {

        count +=1;
        r_count+=1;

        //Calculate/Take p,m- step to even next radius step

        //RK4
        //if(ccount==10 && r_count<440){cout << "Here: dm_dr, dp_dr calculated with r and old m/p" << endl;}
        if(euler_method == false){
        double m_old = m_now;
        double p_old = press_now;
        pair<double, double> mass_press_now = RK4_step(m_old, p_old,r, dr);
        m_now = mass_press_now.first;
        press_now = mass_press_now.second;}

        //EulerMethod
        if(euler_method == true)
        {
            pair<double, double> dp_dm_dr = get_gradients[th](m_now, press_now, r);
            dm_dr = dp_dm_dr.first;
            dp_dr = dp_dm_dr.second;

            m_now = m_now + dm_dr*dr;
            press_now = press_now + dp_dr*dr;
        }
        double rho_now= rho_of_p_num(press_now);



        for(int l=0;l<p_rho_profiles.size();l++)
        {
            if(ccount == l )
            {
                fprintf(p_rho_profiles[l],"%.10lf %.10lf %.10lf\n",press_now,rho_now,r);
            }
        }

        r=r+dr;

        //Adapt dr-step if wanted
        if(dr_adapted==true)
        {   //cout << "next come beast for dr_adapted " << endl;
            pair<double, double> dp_dm_dr = get_gradients[th](m_now, press_now, r);
            runge_step=0;
            dm_dr = dp_dm_dr.first;
            dp_dr = dp_dm_dr.second;

            dr = abs(delta*1./((1.0/m_now)*dm_dr - (1.0/press_now)*dp_dr));
        }


    }

    return make_pair(m_now/msun, r*pow(10,0));

}







//////////////////////////////////////////
//NEXT MASS and PRESS in next RADIUS STEP: not by Newtonian p_i = p_(i-1) + dp/dr*dr, but by 4 intermediate steps --> RK4
/////////////////////////////////////////

pair<double, double> RK4_step(double mass_old, double press_old, double r, double dr)
{
    double mass_temp;
    double press_temp;
    pair<double, double> dm_dp_dr;

    //k1
    mass_temp = mass_old;
    press_temp = press_old;
    dm_dp_dr = get_gradients[th](mass_temp, press_temp, r);
    double m_k1 = dm_dp_dr.first*dr;
    double p_k1 = dm_dp_dr.second*dr;
    runge_step++;


    //k2
    mass_temp = mass_old + 0.5*m_k1;
    press_temp = press_old + 0.5*p_k1;
    dm_dp_dr = get_gradients[th](mass_temp, press_temp, r+0.5*dr);
    double m_k2 = dm_dp_dr.first*dr;
    double p_k2 = dm_dp_dr.second*dr;
    runge_step++;

    //k3
    mass_temp = mass_old + 0.5*m_k2;
    press_temp = press_old + 0.5*p_k2;
    dm_dp_dr = get_gradients[th](mass_temp, press_temp, r+0.5*dr);
    double m_k3 = dm_dp_dr.first*dr;
    double p_k3 = dm_dp_dr.second*dr;
    runge_step++;

    //k4
    mass_temp = mass_old + m_k3;
    press_temp = press_old + p_k3;
    dm_dp_dr = get_gradients[th](mass_temp, press_temp, r+dr);
    double m_k4 = dm_dp_dr.first*dr;
    double p_k4 = dm_dp_dr.second*dr;

    double next_mass = mass_old + (m_k1 + m_k2 + m_k2 + m_k3 + m_k3 + m_k4) / 6;
    double next_press = press_old + (p_k1 + p_k2 + p_k2 + p_k3 + p_k3 + p_k4) / 6;
    runge_step++;

    return make_pair(next_mass, next_press);
}






