#include<cmath>
#include"Integral.h"

double integral_q(double (*f)(double, double, double, double, double, double, double), double m0, double q_max, double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    double dq = 0.2;
    double integral_q_res = 0;
    long long n = q_max/dq;
    double fa = f(m0 ,z, lambda_e, lambda_Pi, gamma_1, gamma_2, 0);
    double fb = f(m0, z, lambda_e, lambda_Pi, gamma_1, gamma_2, q_max);
    integral_q_res += fa;
    integral_q_res += fb;

    //simpson's rule
    long long j = 0;
    for(long long i = 1; i < n; i++)
    {
        double q_i = i*dq;
        double q_imhalf = dq*double(i+j)/2;
        j = i;
        double f_i = f(m0, z, lambda_e, lambda_Pi, gamma_1, gamma_2, q_i);
        double f_imhalf = f(m0, z, lambda_e, lambda_Pi, gamma_1, gamma_2, q_imhalf);
        integral_q_res += (4*f_imhalf);
        integral_q_res += (2*f_i);
    }

    integral_q_res += 4*f(m0, z, lambda_e, lambda_Pi, gamma_1, gamma_2, q_max-dq/2);
    double h = q_max/n;
    return integral_q_res*h/6;
}

double integral_z_for_e(double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    double z_max = 1;
    double dz = 0.05;
    long long m = z_max/dz;
    double m0 = 1.;
    double q_max = 1000.*m0;
    double integral_z_res = 0;
    double fa = integral_q(dedqdz_withm0, m0, q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2)
                -integral_q(dedqdz_withm0, 0., q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2)
                +analytic_dedz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2);
    //std::cout << integral_q(dedqdz_withm0, m0, q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    //std::cout << integral_q(dedqdz_withm0, 0., q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    //std::cout << analytic_dedz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    //std::cout << fa <<std::endl;
    double fb = integral_q(dedqdz_withm0, m0, q_max , z_max, lambda_e, lambda_Pi, gamma_1, gamma_2)
                -integral_q(dedqdz_withm0, 0., q_max , z_max, lambda_e, lambda_Pi, gamma_1, gamma_2)
                +analytic_dedz_withoutm0(z_max, lambda_e, lambda_Pi, gamma_1, gamma_2);
    //std::cout << fb <<std::endl;    
    integral_z_res += fa;
    integral_z_res += fb;
    //simpson's rule
    long long j = 0;
    for(long long i = 1; i < m; i++)
    {
        double z_i = i*dz;
        double z_imhalf = dz*double(i+j)/2;
        j = i;
        double f_i = integral_q(dedqdz_withm0, m0, q_max, z_i, lambda_e, lambda_Pi, gamma_1, gamma_2)
                     - integral_q(dedqdz_withm0, 0., q_max, z_i, lambda_e, lambda_Pi, gamma_1, gamma_2)
                     + analytic_dedz_withoutm0(z_i, lambda_e, lambda_Pi, gamma_1, gamma_2);
        double f_imhalf = integral_q(dedqdz_withm0, m0, q_max, z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2)
                         - integral_q(dedqdz_withm0, 0., q_max, z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2)
                         + analytic_dedz_withoutm0(z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2);
        integral_z_res += (4*f_imhalf);
        integral_z_res += (2*f_i);
    }
    integral_z_res += 4*(integral_q(dedqdz_withm0, m0, q_max, z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2)
                         - integral_q(dedqdz_withm0, 0., q_max, z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2)
                        + analytic_dedz_withoutm0(z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2));
    double h = z_max/m;
    return integral_z_res*h/6;
}

double integral_z_for_P(double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    double z_max = 1;
    double dz = 0.05;
    long long m = z_max/dz;
    double m0 = 1.;
    double q_max = 1000.*m0;
    double integral_z_res = 0;
    double fa = integral_q(dPdqdz_withm0, m0, q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2)
                -integral_q(dPdqdz_withm0, 0., q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2)
                +analytic_dPdz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2);
    //std::cout << integral_q(dPdqdz_withm0, m0, q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    //std::cout << integral_q(dPdqdz_withm0, 0., q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    //std::cout << analytic_dPdz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    //std::cout << fa <<std::endl;
    double fb = integral_q(dPdqdz_withm0, m0, q_max , z_max, lambda_e, lambda_Pi, gamma_1, gamma_2)
                -integral_q(dPdqdz_withm0, 0., q_max , z_max, lambda_e, lambda_Pi, gamma_1, gamma_2)
                +analytic_dPdz_withoutm0(z_max, lambda_e, lambda_Pi, gamma_1, gamma_2);
    integral_z_res += fa;
    integral_z_res += fb;
    //simpson's rule
    long long j = 0;
    for(long long i = 1; i < m; i++)
    {
        double z_i = i*dz;
        double z_imhalf = dz*double(i+j)/2;
        j = i;
        double f_i = integral_q(dPdqdz_withm0, m0, q_max, z_i, lambda_e, lambda_Pi, gamma_1, gamma_2)
                     - integral_q(dPdqdz_withm0, 0., q_max, z_i, lambda_e, lambda_Pi, gamma_1, gamma_2)
                     + analytic_dPdz_withoutm0(z_i, lambda_e, lambda_Pi, gamma_1, gamma_2);
        double f_imhalf = integral_q(dPdqdz_withm0, m0, q_max, z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2)
                         - integral_q(dPdqdz_withm0, 0., q_max, z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2)
                         + analytic_dPdz_withoutm0(z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2);
        integral_z_res += (4*f_imhalf);
        integral_z_res += (2*f_i);
    }
    integral_z_res += 4*(integral_q(dPdqdz_withm0, m0, q_max, z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2)
                         - integral_q(dPdqdz_withm0, 0., q_max, z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2)
                        + analytic_dPdz_withoutm0(z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2));
    double h = z_max/m;
    return integral_z_res*h/6;
}

double integral_z_for_pi1(double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    double z_max = 1;
    double dz = 0.05;
    long long m = z_max/dz;
    double m0 = 1.;
    double q_max = 1000.*m0;
    double integral_z_res = 0;
    double fa = integral_q(dpi1dqdz_withm0, m0, q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2)
                -integral_q(dpi1dqdz_withm0, 0., q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2)
                +analytic_dpi1dz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2);
    //std::cout << integral_q(dpi1dqdz_withm0, 0., q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    //std::cout << analytic_dpi1dz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    double fb = integral_q(dpi1dqdz_withm0, m0, q_max , z_max, lambda_e, lambda_Pi, gamma_1, gamma_2)
                -integral_q(dpi1dqdz_withm0, 0., q_max , z_max, lambda_e, lambda_Pi, gamma_1, gamma_2)
                +analytic_dpi1dz_withoutm0(z_max, lambda_e, lambda_Pi, gamma_1, gamma_2);
    integral_z_res += fa;
    integral_z_res += fb;
    //simpson's rule
    long long j = 0;
    for(long long i = 1; i < m; i++)
    {
        double z_i = i*dz;
        double z_imhalf = dz*double(i+j)/2;
        j = i;
        double f_i = integral_q(dpi1dqdz_withm0, m0, q_max, z_i, lambda_e, lambda_Pi, gamma_1, gamma_2)
                     - integral_q(dpi1dqdz_withm0, 0., q_max, z_i, lambda_e, lambda_Pi, gamma_1, gamma_2)
                     + analytic_dpi1dz_withoutm0(z_i, lambda_e, lambda_Pi, gamma_1, gamma_2);
        double f_imhalf = integral_q(dpi1dqdz_withm0, m0, q_max, z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2)
                         - integral_q(dpi1dqdz_withm0, 0., q_max, z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2)
                         + analytic_dpi1dz_withoutm0(z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2);
        integral_z_res += (4*f_imhalf);
        integral_z_res += (2*f_i);
    }
    integral_z_res += 4*(integral_q(dpi1dqdz_withm0, m0, q_max, z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2)
                         - integral_q(dpi1dqdz_withm0, 0., q_max, z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2)
                        + analytic_dpi1dz_withoutm0(z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2));
    double h = z_max/m;
    return integral_z_res*h/6;
}

double integral_z_for_pi2(double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    double z_max = 1;
    double dz = 0.05;
    long long m = z_max/dz;
    double m0 = 1.;
    double q_max = 1000.*m0;
    double integral_z_res = 0;
    double fa = integral_q(dpi2dqdz_withm0, m0, q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2)
                -integral_q(dpi2dqdz_withm0, 0., q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2)
                +analytic_dpi2dz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2);
    //std::cout << integral_q(dpi2dqdz_withm0, 0., q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    //std::cout << analytic_dpi2dz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    double fb = integral_q(dpi2dqdz_withm0, m0, q_max , z_max, lambda_e, lambda_Pi, gamma_1, gamma_2)
                -integral_q(dpi2dqdz_withm0, 0., q_max , z_max, lambda_e, lambda_Pi, gamma_1, gamma_2)
                +analytic_dpi2dz_withoutm0(z_max, lambda_e, lambda_Pi, gamma_1, gamma_2);
    integral_z_res += fa;
    integral_z_res += fb;
    //simpson's rule
    long long j = 0;
    for(long long i = 1; i < m; i++)
    {
        double z_i = i*dz;
        double z_imhalf = dz*double(i+j)/2;
        j = i;
        double f_i = integral_q(dpi2dqdz_withm0, m0, q_max, z_i, lambda_e, lambda_Pi, gamma_1, gamma_2)
                     - integral_q(dpi2dqdz_withm0, 0., q_max, z_i, lambda_e, lambda_Pi, gamma_1, gamma_2)
                     + analytic_dpi2dz_withoutm0(z_i, lambda_e, lambda_Pi, gamma_1, gamma_2);
        double f_imhalf = integral_q(dpi2dqdz_withm0, m0, q_max, z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2)
                         - integral_q(dpi2dqdz_withm0, 0., q_max, z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2)
                         + analytic_dpi2dz_withoutm0(z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2);
        integral_z_res += (4*f_imhalf);
        integral_z_res += (2*f_i);
    }
    integral_z_res += 4*(integral_q(dpi2dqdz_withm0, m0, q_max, z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2)
                         - integral_q(dpi2dqdz_withm0, 0., q_max, z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2)
                        + analytic_dpi2dz_withoutm0(z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2));
    double h = z_max/m;
    return integral_z_res*h/6;
}
