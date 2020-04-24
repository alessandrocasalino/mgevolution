// Alessandro Casalino
//
// Compile with (Mac with default homebrew gsl 2.6) gcc-9 -O2 quintessence_evolve.c -o quintessence_evolve.exe -L/usr/local/Cellar/gsl/2.6/lib -I/usr/local/Cellar/gsl/2.6/include -lgsl
// Run with ./quintessence_evolve.exe

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


// PHYSICAL PARAMETERS VALUES

// INITIAL CONDITIONS
// Initial value of the scale factor
double a_init = 1e-15;

// Values of fractional density for the cosmological matter today
double Omega_rad_0 = 8e-5;
double Omega_b_0 = 0.0486;
double Omega_Lambda_0 = 0.;//0.6911;
double Omega_cdm_0 = 0.2589;

// Physical constants for output conversions
double _G_ = 6.67428e-11;                 /**< Newton constant in m^3/Kg/s^2 */
double _MPc_over_m_ = 3.085677581282e22;  // Conversion factor from Mpc to m
double _Gyr_over_Mpc_ = 3.06601394e2;     // Conversion factor from Gyr to Mpc
double _c_ = 2.99792458e8;                // Speed of light in m/s
double _H0_ = 67.74;                      // H0 of LCDM in Km/s/Mpc

double fourpiG = 1.;

double T_CONF = 1.;

// COMPUTATION PARAMETERS VALUES

// Number of points used in the computation
int points = (int) 1e6;

// Raise this value to make the csv file smaller, but decreasing resolution
// The value inserted is the ratio between the number of values written in a full resolution file / decreased resolution file
int csv_resolution = 10;

double delta_bisection = 1e-2;
double mg_field_init_min = 1e-20;
double mg_field_init_max = 2e1; //1e16 max for LCDM (model 1)


// QUINTESSENCE PARAMETERS
// Initial value of the quintessence velocity (phi derivated with respect to tau)
double mg_field_p_bk_0 = 0.;

// TEST MODE
// Provides some informations on the terminal and the .csv file during the computation (1: on, others: off)
int TEST_MODE = 1;


// Definition of the POTENTIAL
// For a list of the models see above
double mg_pot(const double mg_field_bk) {

  return mg_field_bk * mg_field_bk / 2.;

}

// Definitions of DIFFERENTIAL EQUATION system
// Need to have first order equations to use Runge Kutta
//
double a_p_rk4(const double a, const double a_p, const double mg_field_bk, const double mg_field_p_bk)
{
	return a_p;
}
double a_pp_rk4(const double a, const double a_p, const double mg_field_bk, const double mg_field_p_bk)
{
  double a3 = a * a * a;
  double rhoa3 = (Omega_cdm_0 + Omega_b_0) + (Omega_Lambda_0 * a3) + (Omega_rad_0 / a);
  double Pa3 = - (Omega_Lambda_0 * a3) + 1./3. * (Omega_rad_0 / a );

  return fourpiG / 3. * ( (rhoa3 - 3. * Pa3) - a * mg_field_p_bk * mg_field_p_bk + 4. * a3 * mg_pot(mg_field_bk) );
}
double mg_field_bk_p_rk4(const double a, const double a_p, const double mg_field_bk, const double mg_field_p_bk)
{
	return mg_field_p_bk;
}
double mg_field_bk_pp_rk4(const double a, const double a_p, const double mg_field_bk, const double mg_field_p_bk)
{
	return - 2. * a_p/a * mg_field_p_bk - mg_pot(mg_field_bk) * a * a;
}
double Hconf(const double a, const double mg_field_bk, const double mg_field_p_bk)
{
	return sqrt((2. * fourpiG / 3.) * ( ((Omega_cdm_0 + Omega_b_0) / a) + (Omega_Lambda_0 * a * a) + (Omega_rad_0 / a / a) + (mg_field_p_bk * mg_field_p_bk / 2.) + (a * a * mg_pot(mg_field_bk)) ));
}

// Integrand for the particle horizon integral
double particleHorizonIntegrand(double a, double mg_field_bk, double mg_field_p_bk)
{
  //return 2. / (sqrt(a) * Hconf(a, mg_field_bk, mg_field_p_bk));
	return 1. / ( a * Hconf(a, mg_field_bk, mg_field_p_bk) );
}
// Particle horizon integral step
double particleHorizon(const int i, double * a, double * mg_field_bk, double * mg_field_p_bk) {
  double h = a[i]-a[i-1];
  double fa = particleHorizonIntegrand(a[i-1],mg_field_bk[i-1],mg_field_p_bk[i-1]);
  double fb = particleHorizonIntegrand(a[i],mg_field_bk[i],mg_field_p_bk[i]);

  return h*(fa+fb)/2.;
}

// Function used to print results stored in vectors as a csv file
void csv(double * t, double * a, double * a_p, double * mg_field_bk, double * mg_field_p_bk, double * particleHorizonVec, char * filename) {

    FILE *fp;
    fp = fopen (filename, "w+");
    fprintf(fp, "%s, %s, %s, %s, %s, %s, %s, %s", "t", "a(t)", "H(t)/H0", "H_prime(t)/H0^2", "Omega_df", "Omega_r", "Omega_b", "Omega_cdm");
    fprintf(fp, ", %s", "omega_df");
    fprintf(fp, ", %s", "PH"); // Particle horizon
    fprintf(fp, ", %s", "mg_field");
    fprintf(fp, ", %s", "mg_field_p");
    if(TEST_MODE==1) fprintf(fp, ", %s", "H_check(t)/H0");
    fprintf(fp, "\n");

    // Time conversion factor
    // double tcf = 1./_H0_/(60.*60.*24.*365.*1e9)*_MPc_over_m_/1000.;

    int i = 0;
    double particleHorizonRes = 0.;

    for(i=1;i<points;i++){
      particleHorizonRes += particleHorizon(i,a,mg_field_bk,mg_field_p_bk);
      particleHorizonVec[i] = particleHorizonRes;
    }

    i = 0;

    double H_test, H_cosmo = 0.;

    while( a[i] <= 1.0){

      double H = Hconf(a[i],mg_field_bk[i],mg_field_p_bk[i]); // this is the Hubble constant with conformal time
      double H_cosmo = H / a[i];
      double H_prime = 1. / a[i] * a_pp_rk4(a[i], a_p[i], mg_field_bk[i], mg_field_p_bk[i]) - H * H; // this is Hubble prime with conformal time
      if(T_CONF==0) H_prime = ( H_prime - H * H )/a[i]/a[i]; // this is Hubble prime with cosmological time
      if(T_CONF==0) H = H_cosmo; // this it the Hubble constant with cosmological time

      fprintf(fp, "%e, %e, %e, %e, %e, %e, %e, %e", t[i], a[i], H / sqrt(2. * fourpiG / 3.), H_prime / sqrt(2. * fourpiG / 3.) / sqrt(2. * fourpiG / 3.), 2. * fourpiG / 3. * ( (mg_field_p_bk[i] * mg_field_p_bk[i] / 2.) + (a[i] * a[i] * mg_pot(mg_field_bk[i])) ) /H_cosmo /H_cosmo , 2. * fourpiG / 3. * Omega_rad_0 /a[i] /a[i] /a[i] /a[i] /H_cosmo  / H_cosmo , 2. * fourpiG / 3. * Omega_b_0 /a[i] /a[i] /a[i] /H_cosmo  / H_cosmo , 2. * fourpiG / 3. * Omega_cdm_0 /a[i] /a[i] /a[i] /H_cosmo  / H_cosmo );
      fprintf(fp, ", %e", (-mg_pot(mg_field_bk[i])+mg_field_p_bk[i]*mg_field_p_bk[i]/2./a[i]/a[i])/(mg_pot(mg_field_bk[i])+mg_field_p_bk[i]*mg_field_p_bk[i]/2./a[i]/a[i]));

      fprintf(fp, ", %e", particleHorizonVec[i]);

      fprintf(fp, ", %e", mg_field_bk[i]);
      fprintf(fp, ", %e", mg_field_p_bk[i]);

      if(TEST_MODE==1){

        H_test = a_p[i]/a[i];
        if(T_CONF==0) H_test = H_test / a[i];

        fprintf(fp, ", %e", H_test / sqrt(2. * fourpiG / 3.));

      }

      fprintf(fp, "\n");

      if(a[i+csv_resolution]<= 1.0){
        i=i+csv_resolution;
      }
      else{
        i++;
      }

    }

    fclose(fp);

}

int scan_for_a0 (double * a) {

  int i=0;

  for(i=0;a[i] <= 1.0;i++);

  return i;

}

// This divides an interval in a logarithmic scale of basis 10, and store results in input vector v
void logscale10 (double * v, double A, double B, int points){

    int i;

    double a = log10(A);
    double b = log10(B);

    double h = (b - a) / (points-1.0);

    for(i=0;i<points;i++){
        v[i] = pow(10.0, a + i * h);
    }

}

// This function is a temporal step of the Runge-Kutta 4 method:
// evolves the system for a time equal to dtau (h in the program)
void mg_rungekutta4bg(double * f, const double dtau)
{

  double k_F[4], k_A[4], k_rhof[4], k_rhor[4], k_rhob[4];
  double a = f[1];
  double a_p = f[2];
  double mg_field_bk = f[3];
  double mg_field_p_bk = f[4];

	double k1a, k2a, k3a, k4a;
	double k1ap, k2ap, k3ap, k4ap;
	double k1f, k2f, k3f, k4f;
	double k1fp, k2fp, k3fp, k4fp;

	k1a 	= a_p_rk4(a, a_p, mg_field_bk, mg_field_p_bk);
	k1ap 	= a_pp_rk4(a, a_p, mg_field_bk, mg_field_p_bk);
	k1f 	= mg_field_bk_p_rk4(a, a_p, mg_field_bk, mg_field_p_bk);
	k1fp 	= mg_field_bk_pp_rk4(a, a_p, mg_field_bk, mg_field_p_bk);

	k2a 	= a_p_rk4(a + k1a * dtau / 2., a_p + k1ap * dtau / 2., mg_field_bk + k1f * dtau / 2., mg_field_p_bk + k1fp * dtau / 2.);
	k2ap 	= a_pp_rk4(a + k1a * dtau / 2., a_p + k1ap * dtau / 2., mg_field_bk + k1f * dtau / 2., mg_field_p_bk + k1fp * dtau / 2.);
	k2f 	= mg_field_bk_p_rk4(a + k1a * dtau / 2., a_p + k1ap * dtau / 2., mg_field_bk + k1f * dtau / 2., mg_field_p_bk + k1fp * dtau / 2.);
	k2fp 	= mg_field_bk_pp_rk4(a + k1a * dtau / 2., a_p + k1ap * dtau / 2., mg_field_bk + k1f * dtau / 2., mg_field_p_bk + k1fp * dtau / 2.);

	k3a 	= a_p_rk4(a + k2a * dtau / 2., a_p + k2ap * dtau / 2., mg_field_bk + k2f * dtau / 2., mg_field_p_bk + k2fp * dtau / 2.);
	k3ap 	= a_pp_rk4(a + k2a * dtau / 2., a_p + k2ap * dtau / 2., mg_field_bk + k2f * dtau / 2., mg_field_p_bk + k2fp * dtau / 2.);
	k3f 	= mg_field_bk_p_rk4(a + k2a * dtau / 2., a_p + k2ap * dtau / 2., mg_field_bk + k2f * dtau / 2., mg_field_p_bk + k2fp * dtau / 2.);
	k3fp 	= mg_field_bk_pp_rk4(a + k2a * dtau / 2., a_p + k2ap * dtau / 2., mg_field_bk + k2f * dtau / 2., mg_field_p_bk + k2fp * dtau / 2.);

	k4a 	= a_p_rk4(a + k3a * dtau, a_p + k3ap * dtau, mg_field_bk + k3f * dtau, mg_field_p_bk + k3fp * dtau);
	k4ap 	= a_pp_rk4(a + k3a * dtau, a_p + k3ap * dtau, mg_field_bk + k3f * dtau, mg_field_p_bk + k3fp * dtau);
	k4f 	= mg_field_bk_p_rk4(a + k3a * dtau, a_p + k3ap * dtau, mg_field_bk + k3f * dtau, mg_field_p_bk + k3fp * dtau);
	k4fp 	= mg_field_bk_pp_rk4(a + k3a * dtau, a_p + k3ap * dtau, mg_field_bk + k3f * dtau, mg_field_p_bk + k3fp * dtau);

	f[1] += dtau * (k1a + 2. * k2a + 2. * k3a + k4a) / 6.;
	f[2] += dtau * (k1ap + 2. * k2ap + 2. * k3ap + k4ap) / 6.;
	f[3] += dtau * (k1f + 2. * k2f + 2. * k3f + k4f) / 6.;
	f[4] += dtau * (k1fp + 2. * k2fp + 2. * k3fp + k4fp) / 6.;

}

// This function evolves the system with the Runge-Kutta 4 method until t_stop
double rk4(double * t, double * a, double * a_p, double * mg_field_bk, double * mg_field_p_bk, double Omega_f_0) {

    double f[5];
    int i = 0;

    // Call the function for a logarithmic scale of time
    // We don't start from t = 0 to avoid problems with quintessence potentials
    logscale10(t,1e-9,20.,points);

    while( i < points - 1 ){

      f[0] = t[i];
      f[1] = a[i];
      f[2] = a_p[i];
      f[3] = mg_field_bk[i];
      f[4] = mg_field_p_bk[i];

      mg_rungekutta4bg(f, t[i+1]-t[i]);

      a[i+1] = f[1];
      a_p[i+1] = f[2];
      mg_field_bk[i+1] = f[3];
      mg_field_p_bk[i+1] = f[4];
      i++;

    }

    int j = scan_for_a0(a);

    double H        = Hconf(a[j], mg_field_bk[j], mg_field_p_bk[j]);
    double Omega_f  = (2. * fourpiG / 3.) * ( (mg_field_p_bk[j] * mg_field_p_bk[j] / 2.) + (a[j] * a[j] * mg_pot(mg_field_bk[j])) ) /H /H;
    return Omega_f  - Omega_f_0;

}


double bisection (double min, double max, double * t, double * a, double * a_p, double * mg_field_bk, double * mg_field_p_bk, const double Omega_f_0) {

  double C = (min+max) / 2.;

  while(fabs((max-min)/min)>delta_bisection){

    a_p[0] = a_init * Hconf(a[0], mg_field_init_min, mg_field_p_bk[0]);
    mg_field_bk[0] = min;
    double rk4_min = rk4(t, a, a_p, mg_field_bk, mg_field_p_bk, Omega_f_0);
    a_p[0] = a_init * Hconf(a[0], C, mg_field_p_bk[0]);
    mg_field_bk[0] = C;
    double rk4_C = rk4(t, a, a_p, mg_field_bk, mg_field_p_bk, Omega_f_0);

    if(rk4_min*rk4_C>=0) {
      min=C;
    }
    else {
      max=C;
    }

    C=(max+min)/2.;

    if (TEST_MODE == 1) printf("TEST_MODE ON - min: %e , max: %e, C: %e, rk4_min: %e , rk4_C: %e \n", min, max, C, rk4_min, rk4_C);

  }

  double result = (max+min)/2.;

  if(TEST_MODE==1) printf("\n");
  printf("\t-> Result of bisection method is mg_field_bk: %e (internal units).\n", result);
  //if(TEST_MODE==1) printf("\t--> Confront with LCDM value: %e (internal units).\n", (2. * fourpiG / 3.) * ( Omega_Lambda_0 + 6. * OmegaCDM_0 )/ pow(a_init,3.));

  return result;

}

int main() {

    double Omega_f_0 = 1. - Omega_Lambda_0 - Omega_rad_0 - Omega_b_0 - Omega_cdm_0;

    int i = 0;

    printf("\n\t\t----------------------------------------\n\n");

    // Definition of the vector needed for the evolution functions
    double * t; double * a; double * a_p; double * mg_field_bk; double * mg_field_p_bk; double * particleHorizonVec;
    t                   = (double *) malloc(sizeof(double) * points);
    a                   = (double *) malloc(sizeof(double) * points);
    a_p                 = (double *) malloc(sizeof(double) * points);
    mg_field_bk         = (double *) malloc(sizeof(double) * points);
    mg_field_p_bk       = (double *) malloc(sizeof(double) * points);
    particleHorizonVec  = (double *) malloc(sizeof(double) * points);

    if(!t||!a||!a_p||!mg_field_bk||!mg_field_p_bk||!particleHorizonVec){
      printf("Error! The memory cannot be allocated. The program will be terminated.\n");
      exit(1);
    }

    // Initial conditions (tau=0)
    a[0] = a_init;
    mg_field_p_bk[0] = mg_field_p_bk_0;

    printf(" Searching for best initial value for the dark fluid ... \n \n");
    mg_field_bk[0] = bisection(mg_field_init_min, mg_field_init_max, t, a, a_p, mg_field_bk, mg_field_p_bk, Omega_f_0);
    a_p[0] = a_init * Hconf(a_init, mg_field_bk[0], mg_field_p_bk[0]);

    printf("\n\n Evolving the system ...\n");

    rk4(t, a, a_p, mg_field_bk, mg_field_p_bk, Omega_f_0);

    char filename[50];
    sprintf (filename, "mg_bk.csv");
    int last_int = scan_for_a0(a);

    printf("\n RESULTS:\n");
    printf("\t-> H0: %f \n", a_p[last_int]/a[last_int]/sqrt(2. * fourpiG / 3.));
    printf("\t-> number of points %d \n", last_int);
    //printf("\t-> Age of the Universe: %f Gyr\n", t[last_int] /_H0_/(60.*60.*24.*365.*1e9)*_MPc_over_m_/1000.);

    csv(t, a, a_p, mg_field_bk, mg_field_p_bk, particleHorizonVec, filename);

    printf("\n The results are saved in '.csv' files. The name is labelled with the value of c, and the model (m).\n");

    if(TEST_MODE==1) printf("\n TEST_MODE ON: check the values of H in .csv file. They must be equal!");

    printf("\n\t\t----------------------------------------\n");

    double * a_int; double * a_p_int; double * mg_field_bk_int; double * mg_field_p_bk_int; double * particleHorizonVec_int;
    a_int                   = (double *) malloc(sizeof(double) * last_int);
    a_p_int                 = (double *) malloc(sizeof(double) * last_int);
    mg_field_bk_int         = (double *) malloc(sizeof(double) * last_int);
    mg_field_p_bk_int       = (double *) malloc(sizeof(double) * last_int);
    particleHorizonVec_int  = (double *) malloc(sizeof(double) * last_int);

    if(!a_int||!a_p_int||!mg_field_bk_int||!mg_field_p_bk_int||!particleHorizonVec_int){
      printf("Error! The memory cannot be allocated. The program will be terminated.\n");
      exit(1);
    }

    memcpy(a_int, a, last_int * sizeof(double));
    memcpy(a_p_int, a_p, last_int * sizeof(double));
    memcpy(mg_field_bk_int, mg_field_bk, last_int * sizeof(double));
    memcpy(mg_field_p_bk_int, mg_field_p_bk, last_int * sizeof(double));
    memcpy(particleHorizonVec_int, particleHorizonVec, last_int * sizeof(double));

    free(a);free(a_p);free(mg_field_bk);free(mg_field_p_bk);free(particleHorizonVec);

    // Spline interpolation with gsl
    gsl_interp_accel *acc_mg_field = gsl_interp_accel_alloc();
    gsl_spline * spline_mg_field = gsl_spline_alloc(gsl_interp_cspline,last_int);
    gsl_interp_accel *acc_mg_field_p = gsl_interp_accel_alloc();
    gsl_spline * spline_mg_field_p = gsl_spline_alloc(gsl_interp_cspline,last_int);
    gsl_interp_accel *acc_a_p = gsl_interp_accel_alloc();
    gsl_spline * spline_a_p = gsl_spline_alloc(gsl_interp_cspline,last_int);
    gsl_interp_accel *acc_particleHorizon = gsl_interp_accel_alloc();
    gsl_spline * spline_particleHorizon = gsl_spline_alloc(gsl_interp_cspline,last_int);

    gsl_spline_init(spline_mg_field,a_int,mg_field_bk_int,last_int);
    gsl_spline_init(spline_mg_field_p,a_int,mg_field_p_bk_int,last_int);
    gsl_spline_init(spline_a_p,a_int,a_p_int,last_int);
    gsl_spline_init(spline_particleHorizon,a_int,particleHorizonVec_int,last_int);

    double a_eval = 1e-2;
    printf("spline eval: %e %e \n",a_eval,gsl_spline_eval(spline_mg_field,a_eval,acc_mg_field));
    printf("spline eval: %e %e \n",a_eval,gsl_spline_eval(spline_mg_field_p,a_eval,acc_mg_field_p));
    printf("spline eval: %e %e \n",a_eval,gsl_spline_eval(spline_a_p,a_eval,acc_a_p));
    printf("spline eval: %e %e \n",a_eval,gsl_spline_eval(spline_particleHorizon,a_eval,acc_particleHorizon));
    printf("%d %e \n",last_int, particleHorizonVec_int[gsl_interp_bsearch(a_int,a_eval,0,last_int-1)]);

    gsl_spline_free(spline_mg_field);gsl_interp_accel_free(acc_mg_field);
    gsl_spline_free(spline_mg_field_p);gsl_interp_accel_free(acc_mg_field_p);
    gsl_spline_free(spline_a_p);gsl_interp_accel_free(acc_a_p);
    gsl_spline_free(spline_particleHorizon);gsl_interp_accel_free(acc_particleHorizon);
    free(t);free(a_int);free(a_p_int);free(mg_field_bk_int);free(mg_field_p_bk_int);free(particleHorizonVec_int);

    printf("\n");

    exit(0);

}
