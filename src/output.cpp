#include "header/eos_filetoarray.h"

using namespace std;

////////////////////////
// WRITES RESULTS ON FILES
////////////////////////

//Get output (collumns not seperated by "," as in CSV, but by space " ") to array
void output_to_arrays(const char* i_file_name,  vector<double> &rhos, vector<double> &masses, vector<double> &radii )
{
        // removed all the windowing stuff.
    ifstream in;
    in.open(i_file_name);

    double tmp_rho;
    double tmp_mass;
    double tmp_radius;

    if (in.is_open())
    {
        // read each column into an array
        for (int i = 0; i < 364; i++)
        { // what do you do if the file ends early or contains bad information?
            in >> tmp_rho >>
                  tmp_mass >>
                  tmp_radius;
            rhos.push_back(tmp_rho);
            masses.push_back(tmp_mass);
            radii.push_back(tmp_radius);
        }
        in.close();
    }

    else
    {
        cout << "error reading file" << endl;
        return;
    }
}

string DoublePowToString(double Rpp)
{

    char Rp_value[4];
    char Rp_pot[4];

    memset(Rp_value, 0, sizeof(Rp_value));
    snprintf(Rp_value, sizeof(Rp_value), "%g", Rpp);
    std::string strObj4(Rp_value);

    memset(Rp_pot, 0, sizeof(Rp_pot));
    snprintf(Rp_pot, sizeof(Rp_pot), "%g", abs(log10(Rpp) ));
    std::string strObj2(Rp_pot);

    string pow_sign;
    if(log10(Rpp)<0)
    {
        pow_sign="-";
    }
    else
    {
        pow_sign="+";
    }

    cout << Rp_value[0] << endl;
    cout << Rp_value[1] << endl;
    cout << Rp_value[2] << endl;

    string Rp_string=Rp_value;
    Rp_string.append("*10^");
    Rp_string.append(pow_sign);


    Rp_string.append(Rp_pot);



    return Rp_string;

}



string create_theory_ID()
{
    string theory_ID;
    string Rp_string;
    string Rq_string;

    Rp_string = DoublePowToString(Rp);
    Rq_string = DoublePowToString(Rq);


    if(GR_theory==1){theory_ID  = "GR";}
    if(fR_theory==1)
    {
        theory_ID  = "fR_Rp_";
        theory_ID.append(Rp_string);

    }
    if(fRQ_theory==1)
    {
        theory_ID  = "fRQ_Rp_";
        theory_ID.append(Rp_string);
        theory_ID.append("_Rq_");
        theory_ID.append(Rq_string);
    }

    return theory_ID;
}




