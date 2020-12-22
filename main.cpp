/*
=================================================================
Authors: Tatsiana Bardachova, Samkeyat Shohan, Pablo Conrat
Date: 17.12.2020
Description: 1D Radiation-Convection Model
=================================================================
*/


#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include "cplkavg.h"
#include "ascii.h"
#include <omp.h>
#include <chrono>
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
    static const double E_abs; // heating rate from surface [W/m^2]
    static const double sigma; // Stefan–Boltzmann constant [W/m^2 K^4]
    static const double M;  // molar mass of dry air [kg/mol]
    static const double R0; // universal gas constant [J/mol K]
  
    static const int nlayer; // number of layers
    static const int nlevel; // number of levels
    static const int nangle; // number of angles

    static const double dt; // time step length [s]
    static const int n_steps; // number of timesteps [/]
    static const int output_steps; // intervall in which the model produces output [/]
};

const double Consts::kappa = 2.0 / 7.0;
const double Consts::c_air = 1004;
const double Consts::g = 9.80665; 
const double Consts::E_abs = 235.0;
const double Consts::sigma = 5.670373e-8;
const double Consts::M = 0.02896;
const double Consts::R0 = 8.3144;

const int Consts::nlayer = 20; 
const int Consts::nlevel = Consts::nlayer + 1; 
const int Consts::nangle = 30; 

const double Consts::dt = 360.0; 
const int Consts::n_steps = 100000;
const int Consts::output_steps = 500;


/*
=================================================================
Output Functions
=================================================================
*/

// first version of an output function, gets called for one timestep
void output_conv(const float &time, const vector<double> &player,
                 const vector<double> &Tlayer, const vector<double> &theta) {
    
  freopen("output.txt","a",stdout);
  // print at what timestep the model is
  printf("output_conv at %2.1f hours \n", time);
  // print the pressure, temperature and potential temperature of each layer
  for (int i=0; i<Consts::nlayer; ++i) {
    printf("level: %3d p: %6.1f T: %5.1f theta: %5.1f\n",
            i, player[i], Tlayer[i], theta[i]);
  }
  fclose(stdout);
  return;
}

// first version of an output function, gets called for one timestep
void output_rad(const vector<double> &zlevel, const vector<double> &plevel,
                const vector<double> &Tlevel, const vector<double> &E_down, const vector<double> &E_up) {
  // print the altitude, pressure and temperature of each level
  for (int i=0; i<Consts::nlevel; ++i) {
    printf("lvl: %3d z: %6.3f p: %5.1f T: %5.2f E_dn: %6.3f E_up: %6.3f\n",
           i, zlevel[i], plevel[i], Tlevel[i], E_down[i], E_up[i]);
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

double t_to_theta(double Temperature, double conversion_factor) {
  double theta;
  theta = Temperature * conversion_factor;
  return(theta);
}

void theta_to_t(const vector<double> &theta, vector<double> &Tlayer, const vector<double> &conversion_factors) {
  transform(theta.begin(), theta.end(),
            conversion_factors.begin(), Tlayer.begin(), divides<double>());
  return;
}

// barometric formula (includes conversion to km)
double p_to_z(const double &plevel,const double &T_surface) {
    return Consts::c_air * T_surface / Consts::g * (1 - pow(plevel / 1000.0, Consts::R0 / (Consts::c_air * Consts::M))) / 1000.0;
}


/*
=================================================================
Thermodynamics
=================================================================
*/

void thermodynamics(vector<double> &Tlayer, const double &dp, vector<double> &dE,
                    double &T_surface, const vector<double> conversion_factors) {
  
  // heating rate
  for (int i=0; i<Consts::nlayer; ++i){
    Tlayer[i] += dE[i] * Consts::dt * Consts::g / (Consts::c_air * dp * 100.0);
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

// absorption coeficient 
double alpha (const double &tau, const double &mu){
  return 1.0 - exp(- tau / mu);
}

void monochromatic_radiative_transfer(vector<double> &E_down, vector<double> &E_up,
                                      const int &i_rad, const vector<double> &tau, int &nwvl, double* wvl,
                                      vector<double> &mu, const double &dmu,
                                      const vector<double> &Tlayer, const double &T_surface) {

  for (int imu=0; imu<Consts::nangle; ++imu) {

    // boundary conditions
    double L_down = 0.0;
    double L_up = cplkavg(wvl[i_rad], wvl[i_rad+1], T_surface);
    E_up[Consts::nlevel-1] += 2 * M_PI * L_up * mu[imu] * dmu;
    
    for (int ilev=1; ilev<Consts::nlevel; ++ilev) {
      L_down = (1 - alpha(tau[ilev-1], mu[imu])) * L_down + alpha(tau[ilev-1], mu[imu]) * cplkavg(wvl[i_rad], wvl[i_rad+1], Tlayer[ilev-1]);
      E_down[ilev] += 2 * M_PI * L_down * mu[imu] * dmu;
    }
      
    for (int ilev=Consts::nlevel-2; ilev >= 0; --ilev) {
      L_up = (1 - alpha(tau[ilev], mu[imu]))*L_up + alpha(tau[ilev], mu[imu]) * cplkavg(wvl[i_rad], wvl[i_rad+1], Tlayer[ilev]);
      E_up[ilev] += 2 * M_PI * L_up * mu[imu] * dmu;
    }
      
  }

  return;
}

void radiative_transfer(vector<double> &Tlayer, vector<double> &E_down, vector<double> &E_up, vector<double> &dE,
                        vector<double> &mu, const double &dmu, const double &T_surface,
                        vector<vector<double>> &tau, int &nwvl, double* wvl) {
  
  fill(E_down.begin(), E_down.end(), 0.0);
  fill(E_up.begin(), E_up.end(), 0.0);

  #pragma omp parallel for schedule(static)
  for (int i_rad=0; i_rad<nwvl-1; ++i_rad) {
    monochromatic_radiative_transfer(E_down, E_up, i_rad, tau[i_rad], nwvl, wvl, mu, dmu, Tlayer, T_surface);
  }

  for (int i=0; i<Consts::nlayer; ++i){
    dE[i] = E_down[i] - E_down[i+1] + E_up[i+1] - E_up[i];
  }

  dE[dE.size()-1] += Consts::E_abs + E_down[Consts::nlevel-1] - E_up[Consts::nlevel-1];

  return;
}

// calculate optical thickness for the chosen atmospheric composition
void optical_thickness(int nwvl, int nlyr, int ngases, double** gases[], double factors[], 
                       vector<vector<double>> &tau) {
    
    for (int i=0; i<ngases; ++i) {
      for (int iwvl=0; iwvl<nwvl; ++iwvl) {
        for (int ilyr=0; ilyr<nlyr; ilyr++) {
          tau[iwvl][ilyr] += gases[i][iwvl][ilyr] * factors[i];
        }
      }
    }
    
  return;
}


/*
=================================================================
Initialization of model
=================================================================
*/


int main() {
    
  omp_set_num_threads(32);
    
  double dp = 1000.0 / (double) Consts::nlayer;
  double dT = 100.0 / (double) Consts::nlayer;
  double dmu = 1.0 / (double) Consts::nangle;
  double T_surface = 288.0;   

  vector<double> player(Consts::nlayer); // vector of pressures for each layer
  vector<double> plevel(Consts::nlevel); // vector of pressures between the layers
  vector<double> zlevel(Consts::nlevel); // vector of heights between the layers
  vector<double> Tlevel(Consts::nlevel); // vector of temperatures between the layers
  vector<double> Tlayer(Consts::nlayer); // vector of temperatures for each layer
  vector<double> theta(Consts::nlayer);  // vector of pot. temperatures for each layer
  vector<double> conversion_factors(Consts::nlayer); // vector for the conversion factors between t and theta
  vector<double> mu(Consts::nangle); // vector of cosines of zenith angles, characterize direction of radiation
  vector<double> dE(Consts::nlayer); // vector of net radiative fluxes after radiative transfer
  vector<double> E_down(Consts::nlevel); // vector of downgoing thermal irradiances for each layer
  vector<double> E_up(Consts::nlevel);   // vector of upgoing thermal irradiances for each layer
  vector<vector<double>> tau; // 2D vector of optical thickness for every layer and every wavelength

    
  for (int i=0; i<Consts::nlevel; ++i) {
    plevel[i] = dp * (double) i; // the pressure levels are spaced equally between 0 and 1000 hPa 
    zlevel[i] = p_to_z(plevel[i], T_surface); // compute height of pressure levels by barometric formula
      
    // initialize irradiance vectors to 0
    E_down[i] = 0.0; 
    E_up[i] = 0.0;
  }


  for (int i=0; i<Consts::nlayer; ++i) {
    Tlayer[i] = 180.0 + dT * (double) i; // just a first guess for the T-profile for each layer
    player[i] = (plevel[i] + plevel[i+1]) / 2.0; // computation of pressure between the levels
    conversion_factors[i] = pow(1000.0 / player[i], Consts::kappa); // computation of conversion factors
    theta[i] = Tlayer[i] * conversion_factors[i]; // computation of theta for each layer

    dE[i] = 0.0; // initialize net radiances to 0
  }

  for (int i=0; i<Consts::nangle; ++i) {
    mu[i] = dmu / 2.0 + dmu * (double) i; // angles are spaced equally between 0 and pi/2
                                          //mu[i] is the center of the i-interval
  }
    
    
   /*
  =========================================================================================================
   Initialization of optical thickness
   Include: reading tau profiles of individual gases and creation of 2D tau vector for combination of gases
  =========================================================================================================
  */ 

  int nwvl=0, nlyr=Consts::nlayer;  // number of wavelength, number of atmospheric layers
    
  double *wvl = NULL;   // array of wavelengths 
  double **tauCO2 = NULL;  // 2D array of optical thickness for CO2
  double **tauH2O = NULL;  // 2D array of optical thickness for H2O
  double **tauN2O = NULL;  // 2D array of optical thickness for N2O
  double **tauCH4 = NULL;  // 2D array of optical thickness for CH4
  double **tauO3  = NULL;  // 2D array of optical thickness for O3
    
  char tauCO2filename[FILENAME_MAX] = "./lbl.arts/lbl.co2.asc";
  char tauH2Ofilename[FILENAME_MAX] = "./lbl.arts/lbl.h2o.asc";
  char tauN2Ofilename[FILENAME_MAX] = "./lbl.arts/lbl.n2o.asc";
  char tauCH4filename[FILENAME_MAX] = "./lbl.arts/lbl.ch4.asc";
  char tauO3filename[FILENAME_MAX]  = "./lbl.arts/lbl.o3.asc";
    
  ASCII_file2xy2D (tauCO2filename, &nwvl, &nlyr, &wvl, &tauCO2);
  ASCII_file2xy2D (tauH2Ofilename, &nwvl, &nlyr, &wvl, &tauH2O);
  ASCII_file2xy2D (tauN2Ofilename, &nwvl, &nlyr, &wvl, &tauN2O);
  ASCII_file2xy2D (tauCH4filename, &nwvl, &nlyr, &wvl, &tauCH4);
  ASCII_file2xy2D (tauO3filename,  &nwvl, &nlyr, &wvl, &tauO3);
    
  // choose composition of atmosphere
  int ngases = 2;  // number of gases 
  double*** gases = new double**[ngases];
  gases[0] = tauH2O; gases[1] = tauCO2;  // specify individual gases
      
  double* factors = new double[ngases];  // array contains ratio of individual gases
  factors[0] = 1.0; factors[1] = 280.0 / 400.0; 
    
  // initialization of tau as vector<vector<double>>  
  for (int i=0; i<nwvl; ++i) {
    tau.push_back(vector<double>());
    for (int ilyr=0; ilyr<nlyr; ++ilyr) {
      tau[i].push_back(0);
    }
  }
  
  // define optical thickness values for every wavelength and every layer
  optical_thickness(nwvl, nlyr, ngases, gases, factors, tau);
    
      
  /*
  =================================================================
  Model Run
  =================================================================
  */

  // loop over time steps
  for (int i=0; i<=Consts::n_steps; i++) {
      
    //chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
      
    // compute current time in hours from start
    float time = (float) i * Consts::dt / 360; // [hours]
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
    radiative_transfer(Tlayer, E_down, E_up, dE, mu, dmu, T_surface, tau, nwvl, wvl);

    thermodynamics(Tlayer, dp, dE, T_surface, conversion_factors);
      
    //chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
    //double runtime = chrono::duration<double>(end - start).count();
    //cout << "Single iteration runtime: " << runtime << endl;
  }
  
  for (int i=0; i<nwvl; ++i){
      delete[] tauCO2[i]; delete[] tauH2O[i]; delete[] tauN2O[i]; delete[] tauCH4[i]; delete[] tauO3[i];
  }
    
  for (int i=0; i<ngases; ++i){
    for (int iwvl=0; iwvl<nwvl; ++iwvl){
      delete[] gases[i][iwvl];
    }
  }
    
  delete[] tauCO2; delete[] tauH2O; delete[] tauN2O; delete[] tauCH4; delete[] tauO3;
  delete[] wvl; delete[] factors;

  return 0;
}