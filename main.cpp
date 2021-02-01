/*
=================================================================================
Authors: Tatsiana Bardachova, Samkeyat Shohan, Pablo Conrat
Date: 01.02.2021
Description: 1D Radiation-Convection Model 
             Representative wavelenght parametrization for Thermal Radiative Transfer
             Solar Radiative Transfer
                 Water Vapor Feedback
=================================================================================
*/

#include <cstdio>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <netcdf>
#include "./repwvl_V1.0_cpp/repwvl_thermal.h"

using namespace std;


/*
=================================================================
Declaration of Constants & Parameters
=================================================================
*/

class Consts {
public: 
    static const double kappa; // adiabatic exponent [/]
    static const double c_air; // specific heat capacity [J/kg K]
    static const double g;     // gravity acceleration [m/s^2]
    static const double E_0; // heating rate from surface [W/m^2]
    static const double h;     // Planck constant [J/s]
    static const double c;     // speed of light in vacuum [m/s]
    static const double kB;    // Boltzmann constant [J/K]
    static const double sigma; // Stefan–Boltzmann constant [W/m^2 K^4]
    static const double M;     // molar mass of dry air [kg/mol]
    static const double R0;    // universal gas constant [J/mol K]
    static const double Tkelvin; // temperature for converting Kelvin into Celsius, [K]
  
    static const int nlayer; // number of layers
    static const int nlevel; // number of levels
    static const int nangle; // number of angles

    static const float max_dT; // maximal temperature change per timestep for stability [K]
    static const int n_steps;  // number of timesteps [/]
    static const int output_steps; // intervall in which the model produces output [/]
    
    static const double tau_s; // optical thickness in the solar spectral range [1/m] ?
    static const double mu_s; // solar zenith angle [°]
    static const int doublings; // number of doublings in the doubling-adding method [/]
    static const double g_asym; // asymmetry factor [/]
    static const double albedo; // surface albedo [/]
    static const float daytime; // proportion of day with sun [/]
};

const double Consts::kappa = 2.0 / 7.0;
const double Consts::c_air = 1004;
const double Consts::g = 9.80665; 
const double Consts::E_0 = 1361.0;
const double Consts::h = 6.62607e-34;
const double Consts::c = 299792458;
const double Consts::kB = 1.380649e-23;
const double Consts::sigma = 5.670373e-8;
const double Consts::M = 0.02896;
const double Consts::R0 = 8.3144;
const double Consts::Tkelvin = 273.15;

const int Consts::nlayer = 20; 
const int Consts::nlevel = Consts::nlayer + 1; 
const int Consts::nangle = 30; 

const float Consts::max_dT = 5;
const int Consts::n_steps = 2000;
const int Consts::output_steps = 10;

const double Consts::tau_s = 2.1;
const double Consts::mu_s = cos( 60 * M_PI / 180.0 );
const int Consts::doublings = 20;
const double Consts::g_asym = 0.85;
const double Consts::albedo = 0.12;
const float Consts::daytime = 0.5;

/*
=================================================================
Output Functions
=================================================================
*/

// simple version of an output function, gets called for one timestep
void output_conv(const float &time, const vector<double> &player,
                 const vector<double> &Tlayer, const vector<double> &theta) {
    
  freopen("output.txt","a",stdout);

  // print the pressure, temperature, and potential temperature of each layer and the time
  for (int i=0; i<Consts::nlayer; ++i) {
    printf("%d,%f,%f,%f,%f\n",
            i, player[i], Tlayer[i], theta[i], time);
  }
    
  return;
}


/*
=================================================================
 Coordinate Change Functions
=================================================================
*/

void t_to_theta(const vector<double> &Tlayer, vector<double> &theta, const vector<double> &conversion_factors) {
  transform(Tlayer.begin(), Tlayer.end(), 
            conversion_factors.begin(), theta.begin(), multiplies<double>());
  return;
}

double t_to_theta(double temperature, double conversion_factor) {
  double theta;
  theta = temperature * conversion_factor;
  return theta;
}

void theta_to_t(const vector<double> &theta, vector<double> &Tlayer, const vector<double> &conversion_factors) {
  transform(theta.begin(), theta.end(), 
            conversion_factors.begin(), Tlayer.begin(), divides<double>());
  return;
}
    

/*
=================================================================
Thermodynamics
=================================================================
*/

void calculate_timestep(const double &dp, vector<double> &dE, double &timestep) {
  timestep = Consts::max_dT / *max_element(dE.begin(), dE.end()) * (Consts::c_air * dp * 100.0) / Consts::g;
  if(timestep > 3600 * 12){
     timestep = 3600 * 12;
  }
  return;
}

void thermodynamics(vector<double> &Tlayer, const double &dp, vector<double> &dE,
                    const double &timestep, double &T_surface, const vector<double> conversion_factors) {
  
  // heating rate
  for (int i=0; i<Consts::nlayer; ++i){
    Tlayer[i] += dE[i] * timestep * Consts::g / (Consts::c_air * dp * 100.0);
  }

  // assume the surface temperature to be the potential temperature of the lowermost layer
  T_surface = t_to_theta(Tlayer[Tlayer.size()-1], conversion_factors[Tlayer.size()-1]);

  return;
}


/*
=================================================================
 Functions for radiative transfer
=================================================================
*/

// Planck function computation, includes wavelength weight implementation
void cplkavg(vector<double> &B, double wvl, const double weight, const vector<double> &Tlayer) {

  wvl *= 1e-9; // from [nm] to [m]
  for (int i=0; i<Consts::nlayer; ++i){
    B[i] = weight * 2 * Consts::h * pow(Consts::c, 2) / (pow(wvl, 5) *
                   (exp(Consts::h * Consts::c / (wvl * Consts::kB * Tlayer[i])) - 1)) / 1e9; 
    // [weight unit] * [W/(m2 nm sterad)]
  }
  return;
}

double cplkavg(double wvl, const double weight, const double &Tlayer) {

  wvl *= 1e-9; // from [nm] to [m]
  double radiance = weight * 2 * Consts::h * pow(Consts::c, 2) / (pow(wvl, 5) * 
                            (exp(Consts::h * Consts::c / (wvl * Consts::kB * Tlayer)) - 1)) / 1e9;
  // [weight unit] * [W/(m2 nm sterad)]    
  return radiance;
}

// absorption coeficient(=emissivity)
void emissivity(vector<double> &alpha, double* tau, const double &mu){
  for (int i=0; i<Consts::nangle; ++i){
    alpha[i] = 1.0 - exp(- tau[i] / mu);
  }
  return;
}

void doubling_adding(double &r_dir, double &s_dir, double &t_dir, double &r, double &t) {
  // assume asymmetry factor g = 0
  double tau = (1 - Consts::g_asym) * Consts::tau_s;
  double dtau = tau / pow(2, Consts::doublings);
  double r_new; double one_minus_rsq;
  double t_new;
  // temporary variables for iterative loop
  double r_dir_new; double s_dir_new; double t_dir_new;

  r = 0.5 * dtau/Consts::mu_s;
  t = 1.0 - r;
  r_dir = dtau/Consts::mu_s * 0.5; // (1-exp(-dtau/mu))
  s_dir = r_dir;
  t_dir = 1 - dtau/Consts::mu_s;
  
  //freopen("output.txt","a",stdout);
  //printf("dtau %f, r %f, t %f, s_dir %f, r_dir %f, t_dir %f \n", dtau, r, t, s_dir, r_dir, t_dir);

  for(int i=0; i < Consts::doublings; ++i){
    one_minus_rsq = (1 - r * r);
    r_new = r + (r * t * t)/one_minus_rsq;
    t_new = (t * t)/one_minus_rsq;

    t_dir_new = pow(t_dir, 2);
    s_dir_new = (t * s_dir + t_dir * r_dir * r * t)/one_minus_rsq + t_dir * s_dir;
    r_dir_new = (t * s_dir * r + t * t_dir * r)/one_minus_rsq + r_dir;

    r = r_new;
    t = t_new;

    t_dir = t_dir_new;
    s_dir = s_dir_new;
    r_dir = r_dir_new;

    dtau = 2 * dtau;
    //printf("iteration: %d, dtau %f, r %f, t %f, r_dir %f, s_dir %f, t_dir %f \n", i, dtau, r, t, r_dir, s_dir, t_dir);
  }

  return;
}

double solar_radiative_transfer_setup(double &r_total, const double &r_dir, const double &s_dir, const double &t_dir, const double &r, const double &t) {
  
  // integrate surface albedo into reflectivity of earth
  r_total = r_dir + (t_dir + s_dir)/(1 - Consts::albedo * r) * t * Consts::albedo;
  
  double solar_irr = Consts::daytime * Consts::E_0 * Consts::mu_s * (1 - r_total);  
  //freopen("output.txt","a",stdout);
  //printf("solar irradiance: %f \n", solar_irr);
  return solar_irr;
}

// Magnus equation for saturated vapour pressure
double magnus(const double &Tlevel) {    
    return 6.1094 * exp (17.625 * (Tlevel - Consts::Tkelvin) / (Tlevel - Consts::Tkelvin + 243.04)); // [hPa];
}

void water_vapor_feedback(vector<double> &e_sat, vector<double> &Tlevel, 
                          const vector<double> &rel_hum, vector<double> &plevel, double* H2O_VMR){
    
     for (int i=0; i<Consts::nlevel; ++i){
         e_sat[i] = magnus(Tlevel[i]);
         H2O_VMR[i] = rel_hum[i] * e_sat[i] / plevel[i];
     }
    return;
}

void monochromatic_radiative_transfer(vector<double> &B, vector<double> &alpha, 
                                      vector<double> &E_down, vector<double> &E_up,
                                      const int &i_rad, double* tau, const double weight, int &nwvl, double* wvl,
                                      vector<double> &mu, const double &dmu,
                                      const vector<double> &Tlayer, const double &T_surface) {

  for (int imu=0; imu<Consts::nangle; ++imu) {
    
    // boundary conditions
    double L_down = 0.0;
    double L_up = cplkavg(wvl[i_rad], weight, T_surface);
    E_up[Consts::nlevel-1] += 2 * M_PI * L_up * mu[imu] * dmu;
    
    emissivity(alpha, tau, mu[imu]);
      
    for (int ilev=1; ilev<Consts::nlevel; ++ilev) {
      L_down = (1 - alpha[ilev-1]) * L_down + alpha[ilev-1] * B[ilev-1];
      E_down[ilev] += 2 * M_PI * L_down * mu[imu] * dmu;
    }
      
    for (int ilev=Consts::nlevel-2; ilev >= 0; --ilev) {
      L_up = (1 - alpha[ilev])*L_up + alpha[ilev] * B[ilev];
      E_up[ilev] += 2 * M_PI * L_up * mu[imu] * dmu;
    }
  
  }
  return;
}

void radiative_transfer(vector<double> &B, vector<double> &alpha,
                        vector<double> &E_down, vector<double> &E_up, vector<double> &dE, const double solar_irr,
                        vector<double> &mu, const double &dmu, 
                        vector<double> &Tlayer, const double &T_surface,
                        double** tau, double* weight, int &nwvl, double* wvl) {
  
  fill(E_down.begin(), E_down.end(), 0.0);
  fill(E_up.begin(), E_up.end(), 0.0);
    
  for (int i_rad=0; i_rad<nwvl; ++i_rad) {
    
    cplkavg(B, wvl[i_rad], weight[i_rad], Tlayer);
      
    monochromatic_radiative_transfer(B, alpha, E_down, E_up, i_rad, 
                                     tau[i_rad], weight[i_rad], nwvl, wvl, mu, dmu, Tlayer, T_surface);
  }
    
  for (int i=0; i<Consts::nlayer; ++i){
    dE[i] = E_down[i] - E_down[i+1] + E_up[i+1] - E_up[i];
  }

  dE[dE.size()-1] += solar_irr + E_down[Consts::nlevel-1] - E_up[Consts::nlevel-1];

  return;
}


/*
=================================================================
Initialization of model
=================================================================
*/

int main() {
      
  double dp = 1000.0 / (double) Consts::nlayer;
  double dmu = 1.0 / (double) Consts::nangle;
  double T_surface = 288.2;   
  double timestep = 0.0;
  float time = 0.0;
  int delete_check = 0;
  double r_dir;
  double s_dir;
  double t_dir;
  double r;
  double t;
  double r_total;

  vector<double> plevel(Consts::nlevel); // vector of pressures between the layers
  vector<double> player(Consts::nlayer); // vector of pressures for each layer
  vector<double> zlevel(Consts::nlevel); // vector of heights between the layers
  vector<double> Tlevel(Consts::nlevel); // vector of temperatures between the layers
  vector<double> Tlayer(Consts::nlayer); // vector of temperatures for each layer
  vector<double> theta(Consts::nlayer);  // vector of pot. temperatures for each layer
  vector<double> conversion_factors(Consts::nlayer); // vector for the conversion factors between t and theta
  vector<double> mu(Consts::nangle); // vector of cosines of zenith angles, characterize direction of radiation
  vector<double> B(Consts::nlayer); // vector of weighted radiance for one wavelength according to Planck's law
  vector<double> alpha(Consts::nangle); // vector of emissivity for one tau value for every angle
  vector<double> dE(Consts::nlayer); // vector of net radiative fluxes after radiative transfer
  vector<double> E_down(Consts::nlevel); // vector of downgoing thermal irradiances for each level
  vector<double> E_up(Consts::nlevel);   // vector of upgoing thermal irradiances for each level  
  vector<double> e_sat(Consts::nlevel);     // vector of saturated vapour pressure for each level
  vector<double> rel_hum(Consts::nlevel);   // vector of relative humidity for each level

  
  /*
  ====================================================================================
  Reading profiles from file (z, p, T and profiles of all trace spacies)
  ====================================================================================
  */
    
  vector<vector<double>> data;   // 2D vector for reading column by column
  ifstream input_file("./repwvl_V1.0_cpp/test.atm");
  string line;
  double value;

  // remove the header (4 lines)
  for (int i = 0; i < 4; i++) {
    getline(input_file, line);
  }
  // allocate all columns (9 columns)
  for (int i = 0; i < 9; i++) {
    data.push_back(vector<double>());
  }
    
  while (getline(input_file, line)) {
    stringstream ss(line);
    int col_index = 0;
      while(ss >> value){
        data[col_index].push_back(value);
        col_index++;
      }
  }
    
  input_file.close();

  zlevel = data[0];
  plevel = data[1];
  Tlevel = data[2];
  
  double* H2O_VMR = data[4].data(); // array of optical thickness profile for HO2
  double* O3_VMR  = data[5].data(); // array of optical thickness profile for O3
  double* CO2_VMR = data[6].data(); // array of optical thickness profile for CO2
  double* CH4_VMR = data[7].data(); // array of optical thickness profile for CH4
  double* N2O_VMR = data[8].data(); // array of optical thickness profile for N2O
  
  // define these species and initialize as zero, need them as arguments in read_tau function
  double* CO_VMR = new double[Consts::nlevel]; // array of optical thickness profile for CO
  double* O2_VMR = new double[Consts::nlevel]; // array of optical thickness profile for O2
  double* HNO3_VMR = new double[Consts::nlevel]; // array of optical thickness profile for HNO3
  double* N2_VMR = new double[Consts::nlevel]; // array of optical thickness profile for N2
  
  // CO2 concentation factor
  double factor = 1;  
    
  for (int i=0; i<Consts::nlevel; ++i) {
    // convert from ppm to absolute concentrations
    H2O_VMR[i] *= 1E-6; O3_VMR[i] *= 1E-6; CO2_VMR[i] *= factor * 1E-6; CH4_VMR[i] *= 1E-6; N2O_VMR[i] *= 1E-6;
    
    // initialize missing species  to 0 
    CO_VMR[i] = 0.0; O2_VMR[i] = 0.0; HNO3_VMR[i] = 0.0; N2_VMR[i] = 0.0;
  }

  for (int i=0; i<Consts::nlevel; ++i) {
    // initialize irradiance vectors to 0
    E_down[i] = 0.0; 
    E_up[i] = 0.0;
      
    e_sat[i] = magnus(Tlevel[i]);
    rel_hum[i] = H2O_VMR[i] * plevel[i] / e_sat[i];
  }

  for (int i=0; i<Consts::nlayer; ++i) {
    player[i] = (plevel[i] + plevel[i+1]) / 2.0; // computation of pressure between the levels
    Tlayer[i] = (Tlevel[i] + Tlevel[i+1]) / 2.0; // computation of temperature between the levels
    conversion_factors[i] = pow(1000.0 / player[i], Consts::kappa); // computation of conversion factors
    theta[i] = Tlayer[i] * conversion_factors[i]; // computation of theta for each layer
    
    B[i]  = 0.0; // initialize radiance to 0
    dE[i] = 0.0; // initialize net radiances to 0
  }

  for (int i=0; i<Consts::nangle; ++i) {
    mu[i] = dmu / 2.0 + dmu * (double) i; // angles are spaced equally between 0 and pi/2
                                          //mu[i] is the center of the i-interval
    alpha[i] = 0.0;  // initialize emissivity to 0
  }

    
  /*
  ==============================================================================================
   Initialization of optical thickness
  ==============================================================================================
  */ 
    
   int nwvl=0;  // number of wavelengths
   double *wvl = NULL;    // array of wavelengths 
   double *weight = NULL; // array of wavelenght weightes
   double **tau = NULL;   // 2D array of optical thicknesses
      
   read_tau("./repwvl_V1.0_cpp/Reduced100Forcing.nc", Consts::nlevel, plevel, Tlevel,
            H2O_VMR, CO2_VMR, O3_VMR, N2O_VMR, CO_VMR, CH4_VMR, O2_VMR, HNO3_VMR, N2_VMR,
            &tau, &wvl, &weight, &nwvl);
    
    
   // computation of plevel for new temperature profile ???
   plevel[Consts::nlevel-1] = plevel[Consts::nlevel-1];
   plevel[Consts::nlevel-2] = (player[Consts::nlayer-1] + plevel[Consts::nlevel-1]) / 2.0;
     for (int ilyr=Consts::nlayer - 2; ilyr>=0; ilyr--) {
       plevel[ilyr] = (player[ilyr+1] + player[ilyr]) / 2.0;
     }
  
  /*
  ====================================================================================
  Setup of solar radiative transfer
  ====================================================================================
  */
  
  freopen("output.txt","a",stdout);
  printf("\n ===================new run===================== \n");
  doubling_adding(r_dir, s_dir, t_dir, r, t);
  double solar_irr = solar_radiative_transfer_setup(r_total, r_dir, s_dir, t_dir, r, t);
  printf("Planetary albedo for an optical thickness of %2.3f and an cos(SZA) of %2.2f: %f \n", Consts::tau_s, Consts::mu_s, r_total);
  printf("rdir %f, sdir %f, and tdir %f \n", r_dir, s_dir, t_dir);
  printf("sum: %f \n", r_dir + t_dir + s_dir);
  printf("solar irradiance: %f \n", solar_irr);
  
  /*
  =================================================================
  Model Run
  =================================================================
  */
  
  // print at what timestep the model is
  printf("layer,player,Tlayer,theta,time\n" );
    
  // loop over time steps
  for (int i=0; i<=Consts::n_steps; i++) {
      
    delete_check = 0; 
      
    // calculate theta values from new T values
    t_to_theta(Tlayer, theta, conversion_factors);

    // sort theta to simulate a stabilizing convection
    sort(theta.begin(), theta.end(), greater<double>());
    theta_to_t(theta, Tlayer, conversion_factors);

    // call output_conv function every n=output_steps times
    if (i % Consts::output_steps == 0) {
      // call output function
      output_conv(time, player, Tlayer, theta);
    }
    
    // recall of optical thickness computations for new temperature profile (every step)
    if (i != 0 and i % 1 == 0) {
        
      for (int iwvl=0; iwvl<nwvl; ++iwvl){
        delete[] tau[iwvl]; 
      }
      delete[] tau; delete[] weight; delete[] wvl;
        
      nwvl=0; 
        
      wvl = NULL;
      weight = NULL; 
      tau = NULL;   
      
      // computation of Tlevel for new temperature profile ???
      Tlevel[Consts::nlevel-1] = T_surface;
      Tlevel[Consts::nlevel-2] = (Tlayer[Consts::nlayer - 1] + T_surface) / 2.0;
        for (int ilyr=Consts::nlayer-2; ilyr>=0; ilyr--) {
          Tlevel[ilyr] = (Tlayer[ilyr+1] + Tlayer[ilyr]) / 2.0;
        }
        
      water_vapor_feedback(e_sat, Tlevel, rel_hum, plevel, H2O_VMR);
      
      read_tau("./repwvl_V1.0_cpp/Reduced100Forcing.nc", Consts::nlevel, plevel, Tlevel,
               H2O_VMR, CO2_VMR, O3_VMR, N2O_VMR, CO_VMR, CH4_VMR, O2_VMR, HNO3_VMR, N2_VMR,
               &tau, &wvl, &weight, &nwvl);
        
      delete_check = 1;
    }
      
    radiative_transfer(B, alpha, E_down, E_up, dE, solar_irr, mu, dmu, Tlayer, T_surface, tau, weight, nwvl, wvl); 
      
    calculate_timestep(dp, dE, timestep);

    thermodynamics(Tlayer, dp, dE, timestep, T_surface, conversion_factors);

    // compute current time in hours from start
    time += (float) timestep / 3600; // [hours]

  } 

  if (delete_check == 0) {
        
      for (int i=0; i<nwvl; ++i){
        delete[] tau[i]; 
      }
      delete[] tau; delete[] weight; delete[] wvl; 
  }
      
  delete[] CO_VMR; delete[] O2_VMR; delete[] HNO3_VMR; delete[] N2_VMR;

  return 0;
}