#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

int main() {
  int nlayer = 10;
  double kappa = 2.0 / 7.0; // adiabatic exponent
  
  double dp = 1000.0 / (double) nlayer;
  double dT = 100.0 / (double) nlayer;
  
  int nlevel = nlayer + 1;
  vector<double> player(nlayer);
  vector<double> p(nlevel);
  vector<double> T(nlayer);
  vector<double> theta(nlayer);

  for (int i=0; i<nlevel; i++) {
    p[i] = dp * (double) i;
  }
  
  for (int i=0; i<nlayer; i++) {
    T[i] = 180.0 + dT * (double) i;
    player[i] = (p[i]+p[i+1])/2.0;
    theta[i] = T[i] * pow(1000.0 / player[i], kappa);
    printf("%3d %6.1f %5.1f %5.1f\n",
           i, player[i], T[i], theta[i]);
  }
  
  /* end of initialization */

  sort(theta.begin(), theta.end(), greater<double>());

  
  for (int i=0; i<nlayer; i++) {
    printf("%3d %6.1f %5.1f %5.1f\n",
           i, player[i], T[i], theta[i]);
  }
  
  return 0;
}