#include <iostream>
#include <string>

// Definition of the MG cosmology strcture
struct mg_cosmology
{

  // Vector of background values to be filled with mg_import function
  std::vector<double> a_vec;
  std::vector<double> H_vec;
  std::vector<double> mg_field_vec;
  std::vector<double> mg_field_p_vec;
  std::vector<double> particleHorizon_vec;

  // Pointers to associated double * arrays to the above vectors (necessary as inputs of GSL interpolation)
  double * a;
  double * H;
  double * mg_field;
  double * mg_field_p;
  double * particleHorizon;

  // Value of the size of the vectors
  int last_int;

  // Interpolation structures (allocated in main)
  gsl_interp_accel * acc_H;
  gsl_spline * spline_H;
  gsl_interp_accel * acc_mg_field;
  gsl_spline * spline_mg_field;
  gsl_interp_accel * acc_mg_field_p;
  gsl_spline * spline_mg_field_p;
  gsl_interp_accel * acc_particleHorizon;
  gsl_spline * spline_particleHorizon;

};

// Function to import the background data from file_name
//  The function will import only the necessary values (from z_in)
bool mg_import(const double a_in, const double fourpiG, mg_cosmology * mg_cosmo, const string file_name) {

    ifstream fin;

    fin.open("mg_bk.csv", ios::in);

    // Counter set to -2 to eliminate header and fist row (bad results)
    int i = -2;

    vector<string> row;
    string line, value;

    do {

        row.clear();

        // Read an entire row and store it in a string variable 'line'
        getline(fin, line);

        // Split the row in words (separated by space)
        stringstream s(line);

        // Populate row vector with different values (eliminating ,)
        getline(s, value, ',');

        while (s >> value) {

            row.push_back(value);
        }

        if(i >= 0){

          if ( stod(row[0]) >= 9. * a_in / 10. ){

            // Populating vectors with file values
            mg_cosmo->a_vec.push_back( stod(row[0]) );
            mg_cosmo->H_vec.push_back( stod(row[1]) * sqrt(fourpiG) );
            mg_cosmo->particleHorizon_vec.push_back( stod(row[8]) / sqrt(fourpiG) );
            mg_cosmo->mg_field_vec.push_back( stod(row[9]) );
            mg_cosmo->mg_field_p_vec.push_back( stod(row[10]) );

          }

        }

        if (i < 0) i++;

    } while (line != "");

    // Set length of vectors for GSL use
    mg_cosmo->last_int = mg_cosmo->a_vec.size()-1;

    // Defining the pointers to standard arrays
    mg_cosmo->a = &mg_cosmo->a_vec[0];
    mg_cosmo->H = &mg_cosmo->H_vec[0];
    mg_cosmo->particleHorizon = &mg_cosmo->particleHorizon_vec[0];
    mg_cosmo->mg_field = &mg_cosmo->mg_field_vec[0];
    mg_cosmo->mg_field_p = &mg_cosmo->mg_field_p_vec[0];

    // Check if the import is good
    if ( mg_cosmo->a_vec.size() > 0 && mg_cosmo->a_vec.size() == mg_cosmo->H_vec.size() && mg_cosmo->a_vec.size() == mg_cosmo->particleHorizon_vec.size() && mg_cosmo->a_vec.size() == mg_cosmo->mg_field_vec.size() && mg_cosmo->a_vec.size() == mg_cosmo->mg_field_p_vec.size()) {
      return true;
    }
    else {
      return false;
    }

}
