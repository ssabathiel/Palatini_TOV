#include "header/eos_filetoarray.h"

using namespace std;

////////////////////////
//GET EOS DATA FROM FILE
////////////////////////


//Get string matrix out of CSV table
void readCSV(istream &input, vector< vector<string> > &output)
{
   string csvLine;
    // read every line from the stream
    while( getline(input, csvLine) )
    {
            istringstream csvStream(csvLine);
            vector<string> csvColumn;
            string csvElement;

            // read every element from the line that is seperated by commas
            // and put it into the vector or strings
            while( getline(csvStream, csvElement, ',') )
            {
                    csvColumn.push_back(csvElement);
            }
            output.push_back(csvColumn);
    }
}



//Write vector on file AND RETURN string MATRIX which, could have been done in readCSV directly
vector< vector<string> > csv_to_array(const char* i_file_name, const char* o_file_name)
{

    ofstream myfile;
    string a;
    fstream file(i_file_name, ios::in);
    myfile.open (o_file_name);
    if(!file.is_open())
    {
           cout << "File not found!\n";
           //return 1;
    }

    // typedef to save typing for the following objecc
    typedef vector< vector<string> > csvVector;
    csvVector csvData;

    readCSV(file, csvData);
    // print out read data to prove reading worked
    for(csvVector::iterator i = csvData.begin(); i != csvData.end(); ++i)
    {
            for(vector<string>::iterator j = i->begin(); j != i->end(); ++j)
            {
              a=*j;
              myfile <<a<<",";
            }
            myfile <<"\n";
    }

    myfile.close();
    return csvData;
}





alglib::spline1dinterpolant string_to_double_EOS_rop(vector< vector<string> > csvData_main,vector<double> &rhos, vector<double> &press)
{

    //First string matrix to double matrix and then to two double vectors
    vector< vector<double> > csvData_double;
    int data_size = csvData_main.size();
    //cout << "Data size= " << data_size << endl;

    // Set up sizes. (HEIGHT x WIDTH)
    csvData_double.resize(data_size);
    for (int i = 0; i < data_size; ++i)
    csvData_double[i].resize(2);

    alglib::real_1d_array rhos_alglib;
    alglib::real_1d_array press_alglib;
    rhos_alglib.setlength(data_size);
    press_alglib.setlength(data_size);

    for(int i=0; i<data_size; i++)
    {
        for(int j=0; j<2; j++)
        {
            string s = csvData_main[i][j];
            csvData_double[i][j] = atof(s.c_str());
            if(j==0)
            {
                rhos[i] = atof(s.c_str());
                rhos_alglib(i) = atof(s.c_str());
                //rhos[i] = log(rhos[i]);
            }
            else
            {
                press[i] = atof(s.c_str());
                press_alglib(i) = atof(s.c_str());
                //press[i] = log(press[i]);
            }
        }
    }


    alglib::spline1dbuildlinear(press_alglib, rhos_alglib, rho_of_press_alg);

    return rho_of_press_alg;

}



alglib::spline1dinterpolant string_to_double_EOS_por(vector< vector<string> > csvData_main,vector<double> &rhos, vector<double> &press)
{

    //First string matrix to double matrix and then to two double vectors
    vector< vector<double> > csvData_double;
    int data_size = csvData_main.size();
    //cout << "Data size= " << data_size << endl;

    // Set up sizes. (HEIGHT x WIDTH)
    csvData_double.resize(data_size);
    for (int i = 0; i < data_size; ++i)
    csvData_double[i].resize(2);

    alglib::real_1d_array rhos_alglib;
    alglib::real_1d_array press_alglib;
    rhos_alglib.setlength(data_size);
    press_alglib.setlength(data_size);

    for(int i=0; i<data_size; i++)
    {
        for(int j=0; j<2; j++)
        {
            string s = csvData_main[i][j];
            csvData_double[i][j] = atof(s.c_str());
            if(j==0)
            {
                //rhos[i] = atof(s.c_str());
                rhos_alglib(i) = atof(s.c_str());
                //rhos[i] = log(rhos[i]);
            }
            else
            {
                //press[i] = atof(s.c_str());
                press_alglib(i) = atof(s.c_str());
                //press[i] = log(press[i]);
            }
        }
    }


    alglib::spline1dinterpolant press_of_rho_alg;
    alglib::spline1dbuildlinear(rhos_alglib, press_alglib, press_of_rho_alg);

    return press_of_rho_alg;

}






