#include <stdio.h>
#include <math.h>
int main() {
    int nlayer=25;
    double kappa = 2.0 / 7.0; // adiabatic exponent
    
    double dp = 1000.0 / (double) nlayer;
    double dT = 100.0 / (double) nlayer;
    freopen("output.csv","w",stdout);
    /*comment*/
    int nlevel=nlayer+1;
    double p[nlevel];
    double player[nlayer];
    double T[nlayer];
    double theta[nlayer];
    
    for (int i=0; i<nlevel; i++) {
        p[i] = dp * (double) i;
        //printf ("%3d %6.1f\n", i, p[i]);
    }
    
    for (int i=0; i<nlevel; i++) {
        T[i] = 180 + dT * (double) i;
        player[i] = (p[i] + p[i+1]) / 2.0;
        theta[i] = T[i] * pow(1000.0 / player[i], kappa);
        printf ("%d %6.1f %5.1f %5.1f\n",
               i,player[i], T[i], theta[i]);
    }
    int unstable = 1;
    int counter =0;
    while (unstable!=0){
        counter ++;
        if (counter>50)
            break;
            
        unstable =0;
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

        for (int i=0; i<nlayer; i++) {
            printf ("%d %6.1f %5.1f %5.1f\n",
                   i,player[i], T[i], theta[i]);
        }
    }
    return 0;
}
