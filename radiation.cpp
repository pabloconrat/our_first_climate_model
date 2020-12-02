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
    static const double h;  // Planck constant [J/s]
    static const double c;  // speed of light in vacuum [m/s]
    static const double kB; // Boltzmann constant [J/K]
    static const double sigma; // Stefan–Boltzmann constant [W/m^2 K^4]
    
    static const int nlayer; // number of layers
    static const int nlevel; // number of levels
    static const int nangle; // number of angles
    static const double tau_total; // total optical thickness of whole atmosphere
    // ATTENTION: wavelengths in [cm]!
    static const vector<double> lamdas; // vector of wavelengths dividing the intervals
    static const int nlamda; // number of wavelengths

};

const double Consts::c_air = 1004.0;
const double Consts::g = 9.80665; 
const double Consts::T0 = 288.0;
const double Consts::M = 0.02896;
const double Consts::R0 = 8.3144;
const double Consts::h = 6.62607e-34;
const double Consts::c = 299792458;
const double Consts::kB = 1.380649e-23;
const double Consts::sigma = 5.670373e-8;

const int Consts::nlayer = 10; 
const int Consts::nlevel = Consts::nlayer + 1; 
const int Consts::nangle = 100; 
const double Consts::tau_total = 10; 
const vector<double> Consts::lamdas = {1, 1e10};
const int Consts::nlamda = Consts::lamdas.size();


// first version of an output function, gets called for one timestep
void output(const vector<double> &zlayer, const vector<double> &player,
            const vector<double> &T, const vector<double> &E_down, const vector<double> &E_up) {
  // print the altitude, pressure and temperature of each level
  for (int i=0; i<zlayer.size(); ++i) {
    printf("lvl: %3d z: %6.3f p: %5.1f T: %5.2f E_dn: %6.3f E_up: %6.3f\n",
           i, zlayer[i], player[i], T[i], E_down[i], E_up[i]);
  }
  return;
}

// barometric formula (includes conversion to km)
double p_to_z(const double &p_level) {
  return Consts::c_air * Consts::T0 / Consts::g * (1 - pow(p_level / 1000.0, Consts::R0 / (Consts::c_air * Consts::M))) / 1000.0;
}

// absorption coeficient 
double alpha (const double &tau, const double &mu){
  return 1.0 - exp(- tau / mu);
}

void monochromatic_radiative_transfer(vector<double> &E_down, vector<double> &E_up,
                                      const int &i_rad, const vector<double> &tau,
                                      vector<double> &mu, const vector<double> &T, const double &dmu) {



  // substitute this for loop by for_each?
  for (int imu=0; imu<mu.size(); imu++) {

    // boundary conditions
    double L_down = 0.0;
    double L_up = cplkavg(Consts::lamdas[i_rad], Consts::lamdas[i_rad+1], T[Consts::nlevel-1]);
    E_up[E_up.size()] += 2 * M_PI * L_up * mu[imu] * dmu;
    
    for (int ilev=1; ilev<E_down.size(); ++ilev) {
      L_down = (1 - alpha(tau[ilev - 1], mu[imu])) * L_down + alpha(tau[ilev - 1], mu[imu]) * cplkavg(Consts::lamdas[i_rad], Consts::lamdas[i_rad+1], T[ilev]);
      E_down[ilev] += 2 * M_PI * L_down * mu[imu] * dmu;
    }
    for (int ilev=E_up.size(); ilev >= 0; --ilev) {
      L_up = (1 - alpha(tau[ilev], mu[imu]))*L_up + alpha(tau[ilev], mu[imu]) * cplkavg(Consts::lamdas[i_rad], Consts::lamdas[i_rad+1], T[ilev]);;
      E_up[ilev] += 2 * M_PI * L_up * mu[imu] * dmu;
    }
  }

  return;
}

void radiative_transfer(vector<double> &T, vector<double> &E_down, vector<double> &E_up, 
                        const vector<double> &tau, vector<double> &mu, const double &dmu) {
  

  for (int i_rad=0; i_rad<Consts::nlamda-1; i_rad++) {
    
    
    monochromatic_radiative_transfer(E_down, E_up, i_rad, tau, mu, T, dmu);
  }
  return;
}

int main() {

  double dp = 1000.0 / (double) Consts::nlayer;  
  double dmu = 1.0 / (double) Consts::nangle;   
  double dtau = Consts::tau_total / (double) Consts::nlayer;
  setvbuf(stdout, NULL, _IOLBF, 0);

  vector<double> p(Consts::nlevel); // vector of pressures between the layers
  vector<double> z(Consts::nlevel); // vector of altitudes between the layers
  vector<double> T(Consts::nlayer); // vector of temperatures for each layer, last value is surface temperature
  vector<double> tau(Consts::nlayer); // vector of optical thickness
  vector<double> mu(Consts::nangle); // vector of cosines of zenith angles, characterize direction of radiation
  vector<double> lamda(Consts::nlamda); // vector of wavelengths
  vector<double> B(Consts::nlayer); // vector of spectral radiances for each layer
  vector<double> E_down(Consts::nlayer); // vector of downgoing thermal irradiances for each layer
  vector<double> E_up(Consts::nlayer); // vector of upgoing thermal irradiances for each layer
  vector<double> player(Consts::nlayer); // vector of pressures for each layer
  vector<double> zlayer(Consts::nlayer); // vector of heights for each layer
  vector<double> Radiances(Consts::nlayer); // Initialize vector of radiances

  vector<double> Test_T(Consts::nlevel);
  Test_T = {127.28, 187.09, 212.42, 229.22, 242.03, 252.48, 261.37, 269.13, 276.04, 282.29, 288.00};
  T = Test_T;

  for (int i=0; i<p.size(); i++) {
    p[i] = dp * (double) i; // the pressure levels are spaced equally between 0 and 1000 hPa
    z[i] = p_to_z(p[i]); // altitude levels
  }

  for (int i=0; i<player.size(); i++) {
    player[i] = (p[i]+p[i+1])/2.0; // computation of pressure between the levels
    tau[i] = dtau;  // values of optical thickness for every layer are the same  
    zlayer[i] = p_to_z(p[i]); // compute height of pressure layers by barometric formula

    // Initialize irradiance vectors to 0
    E_down[i] = 0.0; 
    E_up[i] = 0.0;
  }

  for (int i=0; i<mu.size(); i++) {
    mu[i] = dmu / 2.0 + dmu * (double) i; // angles are spaced equally between 0 and pi/2
                                          //mu[i] is the center of the i-interval
  }

  /* end of initialization */

  
  radiative_transfer(T, E_down, E_up, tau, mu, dmu);
    
  output(zlayer, player, T, E_down, E_up);

  return 0;   
}
