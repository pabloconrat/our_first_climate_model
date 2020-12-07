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
    static const double c_air; // specific heat capacity [J/kg K]
    static const double g;  // gravity acceleration [m/s^2]
    static const double T0; // sea level standard temperature [K]
    static const double M;  // molar mass of dry air [kg/mol]
    static const double R0; // universal gas constant [J/mol K]

    static const int nlayer; // number of layers
    static const int nlevel; // number of levels
    static const int nangle; // number of angles
    static const vector<double> lamdas; // vector of wavelengths dividing the intervals
    static const vector<double> total_taus; // vector containing the total optical thicknesses
    static const int nlamda; // number of wavelengths
};

const double Consts::c_air = 1004.0;
const double Consts::g = 9.80665; 
const double Consts::T0 = 288.0;
const double Consts::M = 0.02896;
const double Consts::R0 = 8.3144;

const int Consts::nlayer = 10; 
const int Consts::nlevel = Consts::nlayer + 1; 
const int Consts::nangle = 30; 
const vector<double> Consts::lamdas = {1, 1e6}; // [nm]
const vector<double> Consts::total_taus = {1};
const int Consts::nlamda = Consts::lamdas.size();

// first version of an output function, gets called for one timestep
void output(const vector<double> &zlevel, const vector<double> &plevel,
            const vector<double> &Tlevel, const vector<double> &E_down, const vector<double> &E_up) {
  // print the altitude, pressure and temperature of each level
  for (int i=0; i<Consts::nlevel; ++i) {
    printf("lvl: %3d z: %6.3f p: %5.1f T: %5.2f E_dn: %6.3f E_up: %6.3f\n",
           i, zlevel[i], plevel[i], Tlevel[i], E_down[i], E_up[i]);
  }
  return;
}

// barometric formula (includes conversion to km)
double p_to_z(const double &plevel) {
    return Consts::c_air * Consts::T0 / Consts::g * (1 - pow(plevel / 1000.0, Consts::R0 / (Consts::c_air * Consts::M))) / 1000.0;
}

// absorption coeficient 
double alpha (const double &tau, const double &mu){
  return 1.0 - exp(- tau / mu);
}

void monochromatic_radiative_transfer(vector<double> &E_down, vector<double> &E_up,
                                      const int &i_rad, const vector<double> &tau,
                                      vector<double> &mu, const vector<double> &T, const double &dmu) {

  // substitute this for loop by for_each?
  for (int imu=0; imu<Consts::nangle; ++imu) {

    // boundary conditions
    double L_down = 0.0;
    double L_up = cplkavg(Consts::lamdas[i_rad], Consts::lamdas[i_rad+1], Consts::T0);
    E_up[Consts::nlevel-1] += 2 * M_PI * L_up * mu[imu] * dmu;
    
    for (int ilev=1; ilev<Consts::nlevel; ++ilev) {
      L_down = (1 - alpha(tau[ilev - 1], mu[imu])) * L_down + alpha(tau[ilev - 1], mu[imu]) * cplkavg(Consts::lamdas[i_rad], Consts::lamdas[i_rad+1], T[ilev - 1]);
      E_down[ilev] += 2 * M_PI * L_down * mu[imu] * dmu;
    }
    for (int ilev=Consts::nlevel-2; ilev >= 0; --ilev) {
      L_up = (1 - alpha(tau[ilev], mu[imu]))*L_up + alpha(tau[ilev], mu[imu]) * cplkavg(Consts::lamdas[i_rad], Consts::lamdas[i_rad+1], T[ilev]);
      E_up[ilev] += 2 * M_PI * L_up * mu[imu] * dmu;
    }
  }

  return;
}

void radiative_transfer(vector<double> &T, vector<double> &E_down, vector<double> &E_up, 
                        vector<double> &mu, const double &dmu) {
  
  for (int i_rad=0; i_rad<Consts::nlamda-1; i_rad++) {

    double dtau = Consts::total_taus[i_rad] / Consts::nlayer;
    vector<double> tau(Consts::nlayer, dtau);

    monochromatic_radiative_transfer(E_down, E_up, i_rad, tau, mu, T, dmu);
  }
  for (int i=0, i<Const::nlayer, i++){
  dE[i] = E_down[i] - E_down[i+1] + E_up[i+1] - E_up[i]
  }
  return;
}

int main() {

  double dp = 1000.0 / (double) Consts::nlayer;  
  double dmu = 1.0 / (double) Consts::nangle;   
  setvbuf(stdout, NULL, _IOLBF, 0);

  vector<double> plevel(Consts::nlevel); // vector of pressures between the layers
  vector<double> zlevel(Consts::nlevel); // vector of heights between the layers
  vector<double> Tlevel(Consts::nlevel); // vector of temperatures between the layers
  vector<double> Tlayer(Consts::nlayer); // vector of temperatures for each layer
  vector<double> mu(Consts::nangle); // vector of cosines of zenith angles, characterize direction of radiation
  vector<double> E_down(Consts::nlevel); // vector of downgoing thermal irradiances for each layer
  vector<double> E_up(Consts::nlevel); // vector of upgoing thermal irradiances for each layer
  vector<double> Radiances(Consts::nlayer); // initialize vector of radiances
    
  for (int i=0; i<Consts::nlevel; i++) {
    plevel[i] = dp * (double) i; // the pressure levels are spaced equally between 0 and 1000 hPa 
    zlevel[i] = p_to_z(plevel[i]); // compute height of pressure levels by barometric formula
    tau[i] = dtau;  // values of optical thickness for every layer are the same 
      
  for (int i=0; i<Consts::nlevel; ++i) {
    plevel[i] = dp * (double) i; // the pressure levels are spaced equally between 0 and 1000 hPa 
    zlevel[i] = p_to_z(plevel[i]); // compute height of pressure levels by barometric formula      

    // Initialize irradiance vectors to 0
    E_down[i] = 0.0; 
    E_up[i] = 0.0;
  }
  
  // given initial temperature profile
  Tlevel = {127.28, 187.09, 212.42, 229.22, 242.03, 252.48, 261.37, 269.13, 276.04, 282.29, 288.00};
    
      
  for (int i=0; i<Consts::nlayer; ++i) {
    Tlayer[i] = (Tlevel[i] + Tlevel[i+1]) / 2.0; // T-profile for each layer
  }
    
  for (int i=0; i<Consts::nangle; i++) {
    mu[i] = dmu / 2.0 + dmu * (double) i; // angles are spaced equally between 0 and pi/2
                                          //mu[i] is the center of the i-interval
  }
    
  /* end of initialization */


  radiative_transfer(Tlayer, E_down, E_up, mu, dmu);
    
  output(zlevel, plevel, Tlevel, E_down, E_up);

  return 0;   
}
