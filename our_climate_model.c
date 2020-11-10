#include <stdio.h>
#include <math.h>
#define kappa 2.0/7.0
#define g 9.80665
#define cp 1004.0

double T2theta (double T, double p){
    return T * pow (1000.0 / p, kappa);
}

double theta2T (double T, double p){
    return theta * pow (p / 1000.0, kappa);
}

int main() {
    int nlayer=20;
    //double kappa = 2.0 / 7.0; // adiabatic exponent
    double x=3;
    
    double dp = 1000.0 / (double) nlayer;
    double dT = 100.0 / (double) nlayer;
    freopen("output.csv","w",stdout);
    
    /*comment*/
    int nlevel=nlayer+1;
    double p[nlevel];
    double player[nlayer];
    double T[nlayer];
    double theta[nlayer];
    double delta_t = 15.0 * 60.0; // 15 minutes
    double Eabs = 235; // W/m2, absorbed solar radiation
    
    /* initialize pressure profile */
    for (int i=0; i<nlevel; i++) {
        p[i] = dp * (double) i;
        //printf ("%3d %6.1f\n", i, p[i]);
    }
    
    for (int i=0; i<nlevel; i++) {
        T[i] = 180 + dT * (double) i;
        player[i] = (p[i] + p[i+1]) / 2.0;
        theta[i] = T2theta(T[i],player[i]);
        //theta[i] = T[i] * pow(1000.0 / player[i], kappa);
        printf ("%d %6.1f %5.1f %5.1f\n",
               i,player[i], T[i], theta[i]);
    }
    
    /*end of initialization*/
    for (int j=0; j<10; j++){
        T[nlayer-1] = T[nlayer-1] + g/cp*Eabs/dp * delta_t;
    
    int unstable = 1;
    while (unstable!=0){          
        unstable = 0;
        printf("\n");
        for (int i=nlayer-1; i>=1; i--) {
            if (theta[i]>theta[i-1]){
                // convection = exchange pot. temperatures
                double temp = theta[i];
                theta [i] = theta [i-1];
                theta [i-1] = temp;
                unstable = 1;
            }
          }
        }

        for (int i=0; i<nlayer; i++) {
            T[i] = theta2T(theta[i], player[i]);
            printf ("%d %6.1f %5.1f %5.1f\n",
                   i,player[i], T[i], theta[i]);
        }
    }
    return 0;
}
