#include<iostream>
#include<random>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include"elements.h"
#include"Integral.h"
#include"rootsolve.h"



int func(const gsl_vector * x, void *params,
	gsl_vector * f)
{
	double e_ = ((struct rparams *) params)->e_real; 
	double P_ = ((struct rparams *) params)->P_real;
 	double pi1_ = ((struct rparams *) params)->pi1_real; 
	double pi2_ = ((struct rparams *) params)->pi2_real;
     

	const double lambda_e = gsl_vector_get(x, 0); 
    const double lambda_Pi = gsl_vector_get(x, 1);
    const double gamma_1 =  gsl_vector_get(x, 2);
    const double gamma_2 = gsl_vector_get(x, 3);

    const double e = integral_z_for_e(lambda_e, lambda_Pi, gamma_1, gamma_2) - e_;
    const double P = integral_z_for_P(lambda_e, lambda_Pi, gamma_1, gamma_2) - P_;
    const double pi1 = integral_z_for_pi1(lambda_e, lambda_Pi, gamma_1, gamma_2) - pi1_;
    const double pi2 = integral_z_for_pi2(lambda_e, lambda_Pi, gamma_1, gamma_2) - pi2_;

	gsl_vector_set(f, 0, e);
	gsl_vector_set(f, 1, P);
	gsl_vector_set(f, 2, pi1);
	gsl_vector_set(f, 3, pi2);
    return GSL_SUCCESS;
}

int print_state_f(size_t iter, gsl_multiroot_fsolver * s)
{
	printf("iter = %3u x = % .3f % .3f  % .3f % .3f"
		"  f(x) = % .3e % .3e % .3e % .3e \n",
		iter,
		gsl_vector_get(s->x, 0),
		gsl_vector_get(s->x, 1),
		gsl_vector_get(s->x, 2),
		gsl_vector_get(s->x, 3),
		gsl_vector_get(s->f, 0),
		gsl_vector_get(s->f, 1),
        gsl_vector_get(s->f, 2),
		gsl_vector_get(s->f, 3));
	return 0;
}

int solve(struct rparams p, double x_init[])
{
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
	int status;
	size_t i, iter = 0;

	const size_t n = 4;

    gsl_multiroot_function f = {&func,n, &p };

    gsl_vector *x = gsl_vector_alloc (n);

    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);
    gsl_vector_set (x, 2, x_init[2]);
    gsl_vector_set (x, 3, x_init[3]);



    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, n);
    gsl_multiroot_fsolver_set (s, &f, x);

    //print_state_f (iter, s);

    do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      //print_state_f (iter, s);

      if (status)   /* check if solver is stuck */
        break;

        status =
        gsl_multiroot_test_residual (s->f, 1e-6);
    }while (status == GSL_CONTINUE && iter < 1000);

    //printf ("status = %s\n", gsl_strerror (status));

    //print_state_f (iter, s);
    std::cout << std::fixed << std::setprecision(7) << " "
              << gsl_vector_get(s->x, 0) << "  "
		      <<gsl_vector_get(s->x, 1)  << "  "
		      <<gsl_vector_get(s->x, 2) << "  "
		      <<gsl_vector_get(s->x, 3);


    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
    return 0;
}