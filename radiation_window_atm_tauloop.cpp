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
    static const double T0; // sea level standard temperature [K]
    
    static const int nlayer; // number of layers
    static const int nlevel; // number of levels
    static const int nangle; // number of angles
    static const vector<double> lambdas; // vector of wavelengths dividing the intervals [nm]
    static const int nlambda; // number of wavelengths packages
};

const double Consts::T0 = 288.0;

const int Consts::nlayer = 10; 
const int Consts::nlevel = Consts::nlayer + 1; 
const int Consts::nangle = 20;  
const vector<double> Consts::lambdas = {1000, 8000, 12000, 1000000}; 
const int Consts::nlambda = Consts::lambdas.size() - 1;

// output function addapted for window atmosphere output
void output(const double &tau_total, const vector<double> &E_down, const vector<double> &E_up) {
  
    printf("tau: %6.1f Edn_surface: %6.3f Eup_TOA: %6.3f\n",
           tau_total, E_down[Consts::nlevel - 1], E_up[0]);
  return;
}

// absorption coeficient 
double alpha (const double &tau, const double &mu){
  return 1.0 - exp(- tau / mu);
}

// define tau(wavelengths band, altitude) for deffirent tau_total
void tau_definition(vector<vector<double>> &tau, const double &tau_total_value){
  
  double tau_total[Consts::nlambda] = {tau_total_value, 0, tau_total_value}; // total optical thickness of whole atmosphere for every interval
  for (int i=0; i<Consts::nlambda; ++i){
    for (int ilev=0; ilev<Consts::nlayer; ++ilev){
      double dtau = tau_total[i] / (double) Consts::nlayer;
      tau[i][ilev] = dtau;
    }
  } 
    
  return;
}

void monochromatic_radiative_transfer(vector<double> &E_down, vector<double> &E_up,
                                      const int &i_rad, vector<double> &tau,
                                      vector<double> &mu, const vector<double> &T, const double &dmu) {

  for (int imu=0; imu<Consts::nangle; ++imu) {

    // boundary conditions
    double L_down = 0.0;
    double L_up = cplkavg(Consts::lambdas[i_rad], Consts::lambdas[i_rad+1], Consts::T0);
    E_up[Consts::nlevel-1] += 2 * M_PI * L_up * mu[imu] * dmu;
    
    for (int ilev=1; ilev<Consts::nlevel; ++ilev) {
      L_down = (1 - alpha(tau[ilev - 1], mu[imu])) * L_down + alpha(tau[ilev - 1], mu[imu]) * cplkavg(Consts::lambdas[i_rad], Consts::lambdas[i_rad+1], T[ilev - 1]);
      E_down[ilev] += 2 * M_PI * L_down * mu[imu] * dmu;
    }
    for (int ilev=Consts::nlevel-2; ilev >= 0; --ilev) {
      L_up = (1 - alpha(tau[ilev], mu[imu]))*L_up + alpha(tau[ilev], mu[imu]) * cplkavg(Consts::lambdas[i_rad], Consts::lambdas[i_rad+1], T[ilev]);
      E_up[ilev] += 2 * M_PI * L_up * mu[imu] * dmu;
    }
  }

  return;
}

void radiative_transfer(vector<double> &T, vector<double> &E_down, vector<double> &E_up, 
                        vector<vector<double>> &tau, vector<double> &mu, const double &dmu) {
  
  for (int i_rad=0; i_rad<Consts::nlambda; i_rad++) {      
    monochromatic_radiative_transfer(E_down, E_up, i_rad, tau[i_rad], mu, T, dmu);
  }
    
  return;
}

int main() {
    
  double dmu = 1.0 / (double) Consts::nangle;
  setvbuf(stdout, NULL, _IOLBF, 0);

  vector<double> Tlevel(Consts::nlevel); // vector of temperatures between the layers
  vector<double> Tlayer(Consts::nlayer); // vector of temperatures for each layer
  vector<double> mu(Consts::nangle); // vector of cosines of zenith angles, characterize direction of radiation
  vector<double> taus(101); // vector of total optical thickness of whole atmosphere
  vector<vector<double>> tau; // vector of optical thickness for every layer are the same
  vector<double> E_down(Consts::nlevel); // vector of downgoing thermal irradiances for each layer
  vector<double> E_up(Consts::nlevel); // vector of upgoing thermal irradiances for each layer
    
  for (int i=0; i<Consts::nlevel; i++) {
    // initialize irradiance vectors to 0
    E_down[i] = 0.0; 
    E_up[i] = 0.0;
  }
  
  // given initial temperature profile
  Tlevel = {127.28, 187.09, 212.42, 229.22, 242.03, 252.48, 261.37, 269.13, 276.04, 282.29, 288.00};
    
  for (int i=0; i<Consts::nlayer; i++) {
    Tlayer[i] = (Tlevel[i] + Tlevel[i+1]) / 2.0; // T-profile for each layer
  }
    
  for (int i=0; i<Consts::nangle; i++) {
    mu[i] = dmu / 2.0 + dmu * (double) i; // angles are spaced equally between 0 and pi/2
                                          //mu[i] is the center of the i-interval
  }
  
    
  // initialization of tau as vector<vector<double>>  
  for (int i=0; i<Consts::nlambda; ++i){
    tau.push_back(vector<double>());
    for (int ilev=0; ilev<Consts::nlayer; ++ilev){
      tau[i].push_back(0);
    }
  }
    
  for (int i=0; i<=100; i++) {
    taus[i] = 0.1 * (double) i; // from 0 to 10 with step 0.1
  }
    
  /* end of initialization */
  
    
  // loop over different tau_total
  for (int i_tau=0; i_tau<=100; ++i_tau){ 
  
    tau_definition(tau, taus[i_tau]);
      
    radiative_transfer(Tlayer, E_down, E_up, tau, mu, dmu);
    
    output(taus[i_tau], E_down, E_up);
    
    // initialize irradiance vectors to 0 for new tau_total
    for (int i=0; i<Consts::nlevel; i++) {
      E_down[i] = 0.0; 
      E_up[i] = 0.0;
    }
  }
        
  return 0;   
}
