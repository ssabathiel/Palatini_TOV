#include "header/eos_functions.h"
#include "header/options.h"
#include "linalg_library/spline.h"
#include <iomanip>
#include "linalg_library/interpolation.h"
#include "header/glob_variables.h"
#include "header/eos_filetoarray.h"

using namespace constants;
using namespace std;

///////////////////////////////////////////////////////
//Equation of states: created about another subfunction (not directly the splines pf, rhof) in order to be easily modified to a 1) CONSTANT rho, 2) EXPONENTIAL function 3) or SPLINE
//////////////////////////////////////////////////////

tk::spline pf;       //pf = p(rho)-function
tk::spline rhof;     //rhof = rho(p)-function



double p_of_rho_num(double rho)
{



    double p;
    if(EOS_mode == "Spline"){p=5;}//pf(rho);}      //in Arnetts paper: p=[dynes/cm²], 1 dynes = 10^-5 Newton, but newton includes also Metres --> cm, still to do
    else if(EOS_mode == "Const"){p=5;}//pf(rho);}
    else if(EOS_mode == "Polytrope")
    {   if(rho<rho_1)
        {
            p = K_0*pow(rho,gamma_0);
        }
        else
        {
             p = K_1*pow(rho,gamma_1);
        }
    }
    //p = 1*pow(10,28);
    //return p;
    return alglib::spline1dcalc(press_of_rho_alg, rho);

}

double rho_of_p_num(double p)
{

    double rho=-1;
    double p_copy=p;
    int counter =0;
    if(EOS_mode == "Spline")
    {
        //rho = rhof(p_copy);

    }

    else if(EOS_mode == "Const"){rho = rho_center;}
    else if(EOS_mode == "Polytrope")
    {
        if(p<p_1)
        {
            rho = pow((p/(K_0)),1.0/gamma_0);
        }
        else
        {
           rho = pow((p/(K_1)),1.0/gamma_1);
        }
    }
    //return rho;
    return alglib::spline1dcalc(rho_of_press_alg, p);

}



/*
//Numerical Derivative (=difference quotient)
double drho_dp_function(double p)
{

    //rho1
    double rho1=-1;
    double p_copy=p;
    int counter = 0;
    while(rho1<0 && counter < 5)
    {   counter +=1;
        rho1 = rhof(p_copy);
        if(p_copy<pow(10,12)){rho1=7.86;}
        p_copy=p_copy/2;
    }


    //rho2
    double rho2=-1;
    p_copy=p+1;
    counter = 0;
    while(rho2<0 && counter < 5)
    {   counter +=1;
        rho2 = rhof(p_copy);
        if(p_copy<pow(10,12)){rho2=7.86;}
        p_copy=p_copy/2;
    }


    return rho2-rho1;

}

*/


double drho_dp_num(double p)
{
    int N=10;
    vector<double> line(N,0);
    double x = p;
    double h = p/1.;//0000000;
    double err;
    const int ntab=10; //Sets maximum size of tableau.
    const double con=1.4, con2=(con*con); //Stepsize decreased by CON at each iteration.
    const double big=pow(10,20);
    const double safe=2.0; //Return when error is SAFE worse than the
    int i,j; //best so far.
    double errt,fac,hh,ans;
    vector< vector<double> > aa(N,line);
    // Set up sizes. (HEIGHT x WIDTH)
    aa.resize(ntab);
    for (int i = 0; i < ntab; ++i)
    {
        aa.resize(ntab);
    }

aa[0][0] = 5;
    if (h == 0.0)
    {
        throw("h must be nonzero in dfridr.");
    }
    hh=h;


    aa[0][0]=(rho_of_p_num(x+hh)-rho_of_p_num(x-hh))/(2.0*hh);

    err=big;

    for (i=1;i<ntab;i++)
    {
        hh /= con;
        aa[0][i]=(rho_of_p_num(x+hh)-rho_of_p_num(x-hh))/(2.0*hh); //Try new, smaller stepsize.
        fac=con2;
        for (j=1;j<=i;j++)
        { //Compute extrapolations of various orders, requiring no new function evaluations.
            aa[j][i]=(aa[j-1][i]*fac-aa[j-1][i-1])/(fac-1.0);
            fac=con2*fac;
            errt=max(abs(aa[j][i]-aa[j-1][i]),abs(aa[j][i]-aa[j-1][i-1]));
            //The error strategy is to compare each new extrapolation to one order lower, both
            //at the present stepsize and the previous one.

            if (errt <= err)
            { //If error is decreased, save the improved answer.
                err=errt;
                ans=aa[j][i];
            }
        }
        //If higher order is worse by a significant factor SAFE, then quit early
        if (abs(aa[i][i]-aa[i-1][i-1]) >= safe*err)
        {
        break;
        }


    }

    return ans;
}

void create_spline(char * file_name)
{
    //Get .csv data (from given file) to two vectors --> rhos[], press[]

    const char* write_dat = "other.dat";
    vector< vector<string> > csvData_main;
    csvData_main = csv_to_array(file_name, write_dat);

    //Now have EOS data in STRING MATRIX (size: data_size *2) = (number of shells * (rho & press))    //since Polytropic
    //Now want to have 2*double vectors in data_size in same size
    int data_size = csvData_main.size();
    cout << "Data size= " << data_size << endl;
    vector<double> rhos(data_size);
    vector<double> press(data_size);

    rho_of_press_alg = string_to_double_EOS_rop(csvData_main,rhos,press);
    press_of_rho_alg = string_to_double_EOS_por(csvData_main,rhos,press);

    //Create interpolation-functions of the discrete EOS data-set
    pf.set_points(rhos,press);
    rhof.set_points(press,rhos);

}

double ddrho_dPP_num(double pressy)
{
    double ddrhodPP = 0;
    return ddrhodPP;
}


void create_points_from_analytical_EOS(char *file_name, double (*analytical_fct)(double rho))
{
    FILE *ifp=fopen(file_name,"w");
    rewind(ifp);

    FILE *ifp2=fopen("/home/silvester/fRQ_Projectt/build-fRQ_TOV2-Desktop-Debug/Plotting/EOS_Plotting/Sly_csv","w");
    rewind(ifp2);

    double press;
    double rho = pow(10,6);

    while(rho<7.0*pow(10,15))
    {
        press = analytical_fct(rho);
        fprintf(ifp,"%.10lf,%.10lf\n",rho,press);
        fprintf(ifp2,"%.10lf,%.10lf\n",rho,press);
        rho *= 1.8;
    }
    fclose(ifp);
}
