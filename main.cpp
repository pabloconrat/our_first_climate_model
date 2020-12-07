#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>
#include "cplkavg.h"
using namespace std;

// declaration and initialization of physical constants and model parameters
class Consts {
public: 
    static const double kappa; // adiabatic exponent [/]
    static const double c_air; // specific heat capacity [J/kg K]
    static const double g;     // gravity acceleration [m/s^2]
    static const double T0; // sea level standard temperature [K]
    static const double E_abs; // heating rate from surface [W/m^2]
    static const double sigma; // Stefan–Boltzmann constant [W/m^2 K^4]
    static const double eps;   // emissivity [/]
    static const double M;  // molar mass of dry air [kg/mol]
    static const double R0; // universal gas constant [J/mol K]
  
    static const int nlayer; // number of layers
    static const int nlevel; // number of levels
    static const int nangle; // number of angles

    static const double dt; // time step length [s]
    static const int n_steps; // number of timesteps [/]
    static const int output_steps; // intervall in which the model produces output [/]
    static const double cooling_rate; // prescribed cooling rate [K/s]

    static const vector<double> lamdas; // vector of wavelengths dividing the intervals
    static const vector<double> total_taus; // vector containing the total optical thicknesses
    static const int nlamda; // number of wavelengths
};

const double Consts::kappa = 2.0 / 7.0;
const double Consts::c_air = 1004;
const double Consts::g = 9.80665; 
const double Consts::T0 = 288.0;
const double Consts::E_abs = 235;
const double Consts::sigma = 5.670373e-8;
const double Consts::eps = 0.3;
const double Consts::M = 0.02896;
const double Consts::R0 = 8.3144;

const int Consts::nlayer = 25; 
const int Consts::nlevel = Consts::nlayer + 1; 
const int Consts::nangle = 30; 

const double Consts::dt = 360.0; 
const int Consts::n_steps = 10000;
const int Consts::output_steps = 1000;
const double Consts::cooling_rate = - 1.475/ 86400.0; 

const vector<double> Consts::lamdas = {1, 1e6}; // [nm]
const vector<double> Consts::total_taus = {1};
const int Consts::nlamda = Consts::lamdas.size();


// first version of an output function, gets called for one timestep
void output_conv(const float &time, const vector<double> &player,
            const vector<double> &Tlayer, const vector<double> &theta) {
  // print at what timestep the model is
  printf("output_conv at %2.1f hours \n", time);
  // print the pressure, temperature and potential temperature of each layer
  for (int i=0; i<player.size(); i++) {
    printf("%3d %6.1f %5.1f %5.1f\n",
           i, player[i], Tlayer[i], theta[i]);
  }
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

void t_to_theta(const vector<double> &Tlayer, vector<double> &theta, const vector<double> &conversion_factors) {
  transform(Tlayer.begin(), Tlayer.end(),
            conversion_factors.begin(), theta.begin(), multiplies<double>());
  return;
}

void theta_to_t(const vector<double> &theta, vector<double> &Tlayer, const vector<double> &conversion_factors) {
  transform(theta.begin(), theta.end(),
            conversion_factors.begin(), Tlayer.begin(), divides<double>());
  return;
}

// heating function - so far only surface heating
void thermodynamics(vector<double> &Tlayer, const double &dp) {
  // computation of surface heating
  Tlayer[Tlayer.size()-1]+= Consts::E_abs * Consts::dt * Consts::g / (Consts::c_air * dp * 100.0);
    
  // computation of thermal cooling at the ground (according to Stefan-Boltzmann law)
  Tlayer[Tlayer.size()-1]-= Consts::eps * pow(Tlayer[Tlayer.size()-1], 4) * Consts::sigma * Consts::dt * Consts::g / (Consts::c_air * dp * 100.0); 
    
  // prescribed cooling of atmosphere
  double delta_T = Consts::cooling_rate * Consts::dt;
  transform(Tlayer.begin(), Tlayer.end(), Tlayer.begin(), bind2nd(plus<double>(), delta_T));
  return;
}

// barometric formula (includes conversion to km)
double p_to_z(const double &plevel) {
    return Consts::c_air * Consts::T0 / Consts::g * (1 - pow(plevel / 1000.0, Consts::R0 / (Consts::c_air * Consts::M))) / 1000.0;
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
                                      const int &i_rad, const vector<double> &tau,
                                      vector<double> &mu, const vector<double> &Tlayer, const double &dmu) {
  // substitute this for loop by for_each?
  for (int imu=0; imu<Consts::nangle; ++imu) {

    // boundary conditions
    double L_down = 0.0;
    double L_up = cplkavg(Consts::lamdas[i_rad], Consts::lamdas[i_rad+1], Consts::T0);
    E_up[Consts::nlevel-1] += 2 * M_PI * L_up * mu[imu] * dmu;
    
    for (int ilev=1; ilev<Consts::nlevel; ++ilev) {
      L_down = (1 - alpha(tau[ilev - 1], mu[imu])) * L_down + alpha(tau[ilev - 1], mu[imu]) * cplkavg(Consts::lamdas[i_rad], Consts::lamdas[i_rad+1], Tlayer[ilev - 1]);
      E_down[ilev] += 2 * M_PI * L_down * mu[imu] * dmu;
    }
    for (int ilev=Consts::nlevel-2; ilev >= 0; --ilev) {
      L_up = (1 - alpha(tau[ilev], mu[imu]))*L_up + alpha(tau[ilev], mu[imu]) * cplkavg(Consts::lamdas[i_rad], Consts::lamdas[i_rad+1], Tlayer[ilev]);
      E_up[ilev] += 2 * M_PI * L_up * mu[imu] * dmu;
    }
  }

  return;
}

void radiative_transfer(vector<double> &Tlayer, vector<double> &E_down, vector<double> &E_up, 
                        vector<double> &mu, const double &dmu) {
  
  for (int i_rad=0; i_rad<Consts::nlamda-1; i_rad++) {
    double dtau = Consts::total_taus[i_rad] / Consts::nlayer;
    vector<double> tau(Consts::nlayer, dtau);

    monochromatic_radiative_transfer(E_down, E_up, i_rad, tau, mu, Tlayer, dmu);
  }
  return;
}

/*
=================================================================
Initialization of model
=================================================================
*/

int main() {
  double dp = 1000.0 / (double) Consts::nlayer;
  double dT = 100.0 / (double) Consts::nlayer;
  double dmu = 1.0 / (double) Consts::nangle;   

  vector<double> player(Consts::nlayer); // vector of pressures for each layer
  vector<double> plevel(Consts::nlevel); // vector of pressures between the layers
  vector<double> zlevel(Consts::nlevel); // vector of heights between the layers
  vector<double> Tlevel(Consts::nlevel); // vector of temperatures between the layers
  vector<double> Tlayer(Consts::nlayer); // vector of temperatures for each layer
  vector<double> mu(Consts::nangle); // vector of cosines of zenith angles, characterize direction of radiation
  vector<double> E_down(Consts::nlevel); // vector of downgoing thermal irradiances for each layer
  vector<double> E_up(Consts::nlevel); // vector of upgoing thermal irradiances for each layer
  vector<double> Radiances(Consts::nlayer); // initialize vector of radiances
  vector<double> theta(Consts::nlayer); // vector of pot. temperatures for each layer
  vector<double> conversion_factors(Consts::nlayer); // vector for the conversion factors between t and theta

  for (int i=0; i<Consts::nlevel; ++i) {
    plevel[i] = dp * (double) i; // the pressure levels are spaced equally between 0 and 1000 hPa 
    zlevel[i] = p_to_z(plevel[i]); // compute height of pressure levels by barometric formula      
    // Initialize irradiance vectors to 0
    E_down[i] = 0.0; 
    E_up[i] = 0.0;
  }

  for (int i=0; i<Consts::nlayer; i++) {
    Tlayer[i] = 180.0 + dT * (double) i; // just a first guess for the T-profile for each layer
    player[i] = (plevel[i]+plevel[i+1])/2.0; // computation of pressure between the levels
    conversion_factors[i] = pow(1000.0 / player[i], Consts::kappa); // computation of conversion factors
    theta[i] = Tlayer[i] * conversion_factors[i]; // computation of theta for each layer
  }

  for (int i=0; i<Consts::nangle; i++) {
    mu[i] = dmu / 2.0 + dmu * (double) i; // angles are spaced equally between 0 and pi/2
                                          //mu[i] is the center of the i-interval
  }

  /* end of initialization */

  // loop over time steps
  for (int i=0; i<=Consts::n_steps; i++) {
    // compute current time in hours from start
    float time = (float) i * Consts::dt / 360; // [hours]
    // calculate theta values from new T values
    t_to_theta(Tlayer, theta, conversion_factors);

    // sort theta to simulate a stabilizing mixing
    sort(theta.begin(), theta.end(), greater<double>());
    theta_to_t(theta, Tlayer, conversion_factors);

    // call output_conv function every n=output_steps times
    if (i % Consts::output_steps == 0) {
      // call output function
      output_conv(time, player, Tlayer, theta);
    }

    radiative_transfer(Tlayer, E_down, E_up, mu, dmu);

    thermodynamics(Tlayer, dp);
  }

  return 0;
}
