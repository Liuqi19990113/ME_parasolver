#include<iostream>
#include<cmath>
#include<iomanip>
#include"elements.h"
#include"Integral.h"
#include"rootsolve.h"

int main(int argc, char* argv[])
{

    double e_real = atof(argv[1]);
    double P_real = atof(argv[2]);
    double pi1_real = atof(argv[3]);
    double pi2_real = atof(argv[4]);
    double x_guess[4] = {atof(argv[5]),atof(argv[6]),atof(argv[7]),atof(argv[8])};


    struct rparams real_value = {e_real, P_real, pi1_real, pi2_real};
    solve(real_value, x_guess);

}