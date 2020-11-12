#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

// First version of an output function, gets called for one timestep
void output(const float &time, const vector<double> &player,
            const vector<double> &T, const vector<double> &theta) {
  // Print at what timestep the model is
  printf("Output at %2.1f hours \n", time);
  // print the pressure, temperature and potential temperature of each layer
  for (int i=0; i<player.size(); i++) {
    printf("%3d %6.1f %5.1f %5.1f\n",
           i, player[i], T[i], theta[i]);
  }
  return;
}

// heating function - so far only surface heating
void heating(vector<double> &theta, const double &dp, const double &E_abs,
             const double &dt, const double &g, const double &c_air,
             const vector<double> &t_to_theta) {
  // computation of heating
  double delta_T = E_abs * dt * g / (c_air * dp * 100.0);
  theta[theta.size()-1] += delta_T * t_to_theta[theta.size()-1];
  return;
}

/* Initialization of physical constants */
//'Python' way
struct con_struct 
{
    const double kappa = 2.0 / 7.0; // adiabatic exponent [/]
    const double c_air = 1004; // specific heat capacity [J/kg K]
    const double g = 9.80665; // gravity acceleration [m/s^2]
    const double E_abs = 235; // heating rate from surface [W/m^2]   
};

/*
// 'C++' way
class con 
{
// declaration
public:
    static const double kappa; // adiabatic exponent [/]
    static const double c_air; // specific heat capacity [J/kg K]
    static const double g; // gravity acceleration [m/s^2]
    static const double E_abs; // heating rate from surface [W/m^2]
};
// initialization
const double con::kappa = 2.0 / 7.0;
const double con::c_air = 1004;
const double con::g = 9.80665; 
const double con::E_abs = 235;
*/

int main() {
    
  struct con_struct con; // make an object of struct in 'main' for 'Python' way
    
  int nlayer = 25; // number of layers
  double dt = 360.0; // time step length [s]
  double n_steps = 1000; // number of timesteps [/]
  int output_steps = 100; // Intervall in which the model produces output [/]
  
  double dp = 1000.0 / (double) nlayer;
  double dT = 100.0 / (double) nlayer;

  int nlevel = nlayer + 1;
  vector<double> player(nlayer); // vector of pressures for each layer
  vector<double> p(nlevel); // vector of pressures between the layers
  vector<double> T(nlayer); // vector of temperatures for each layer
  vector<double> theta(nlayer); // vector of pot. temperatures for each layer
  vector<double> t_to_theta(nlayer); // vector for the conversion factors between t and theta

  /* define pressure levels */  
  for (int i=0; i<nlevel; i++) {
    p[i] = dp * (double) i; // the pressure levels are spaced equally between 0 and 1000 hPa
  }

  /* define T, p and theta for each layer */
  for (int i=0; i<nlayer; i++) {
    T[i] = 180.0 + dT * (double) i; // just a first guess for the T-profile for each layer
    player[i] = (p[i]+p[i+1])/2.0; // computation of pressure between the levels
    t_to_theta[i] = pow(1000.0 / player[i], con.kappa); // computation of conversion factors
    theta[i] = T[i] * t_to_theta[i]; // computation of theta for each layer
  }

  /* end of initialization */

  // loop over time steps
  for (int i=0; i<=n_steps; i++) {
    // Compute current time in hours from start
    float time = (float) i * dt / 360; // [hours]

    // Sort theta to simulate a stabilizing mixing
    sort(theta.begin(), theta.end(), greater<double>());

    // call output function every n=output_steps times
    if (i % output_steps == 0) {
      // calculate temperature (vectorized). This is an elementwise division
      transform(theta.begin(), theta.end(),
                t_to_theta.begin(), T.begin(), divides<double>());
      // call output function
      output(time, player, T, theta);
    }

    heating(theta, dp, con.E_abs, dt, con.g, con.c_air, t_to_theta);
  }

  return 0;
}
