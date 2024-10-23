#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
/* function prototype */
/* dfunct - defines the first order differencial equations to solve */
int dfunc (double t, const double y[], double f[], void *params_ptr) {
	double *lparams = (double *) params_ptr;
	/* get parameter(s) from params_ptr */
	double alpha = lparams[0];
	double beta = lparams[1];
    double sigma = lparams[2];
	/* evaluate the right-hand-side functions at t */
	f[0] = (1/beta - alpha)*y[0] + y[2] + y[0]*y[1];
    f[1] = -beta*y[1] - y[0]*y[0];
    f[2] = -y[0] - sigma*y[2];
	return GSL_SUCCESS; /* GSL_SUCCESS defined in gsl/errno.h as 0 */
}
int main () {
	//FILE * outdata = fopen ("output_data.txt","w");
	const int dimension = 3;   /* number of differential equations */
	int status;                /* status of driver function */
	const double eps_abs = 1.e-8;  /* absolute error requested  */
	const double eps_rel = 1.e-10; /* relative error requested  */
	double alpha = 0.00001;	         /* parameter for the diff eq */
	double beta = 0.1; /* parameter for the diff eq */
    double sigma = 1;
	double myparams[3];            /* array for parameters      */
	double y[dimension];	        /* current solution vector */
	double t, t_next;             /* current and next independent variable */
	double tmin, tmax, delta_t;   /* range of t and step size for output */
	double h = 1.0e-6;            /* starting step size for ode solver */
	gsl_odeiv2_system ode_system;	/* structure with the dfunc function, etc. */
	myparams[0] = alpha;          /* problem parameters */
	myparams[1] = beta;
    myparams[2] = sigma;
	//printf("\nThis program solves a system with a single diff  equation\n\n");
	/* load values into the ode_system structure */
	ode_system.function = dfunc;       /* the right-hand-side of equation */
	ode_system.dimension = dimension;  /* number of diffeq's */
	ode_system.params = myparams;      /* parameters to pass to dfunc */
	tmin = 0.0;			/* starting t value */
	//tmax = 365.0;			/* final t value */
    tmax = 1000.0;
	delta_t = 0.001;
	y[0] = 0.1;
    y[1] = 0.23;
    y[2] = 0.31;			/* initial value of x */
	gsl_odeiv2_driver * drv =
	    gsl_odeiv2_driver_alloc_y_new (&ode_system, gsl_odeiv2_step_rkf45, h, eps_abs, eps_rel);
	//printf("Input data: \n");
	//printf(" alpha = %g; beta = %g; sigma = %g\n", alpha, beta, sigma);
	//printf(" Starting step size (h): %0.5e\n", h);
	//printf(" Time parameters: %f %f %f \n", tmin, tmax, delta_t);
	//printf(" Absolute and relative error requested: %0.6e %0.6e \n", eps_abs, eps_rel);
	//printf(" Number of equations (dimension): %d \n\n", dimension);
	//printf("    Time         x          \n");
	t = tmin;             /* initialize t */
	//printf ("%.5e %.5e %.5e %.5e \n", t, y[0], y[1], y[2]);	/* initial values */
	/* step from tmin to tmax */
	for (t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
		{
			status = gsl_odeiv2_driver_apply (drv, &t, t_next, y);
			if (status != GSL_SUCCESS) {
					printf("Error: status = %d \n", status);
					break;
				}
			//printf ("%.5e %.5e %.5e %.5e\n", t, y[0], y[1], y[2]); /* print at t=t_next */
            printf ("%f %f %f %f\n", t, y[0], y[1], y[2]);
			//fprintf(outdata,"%.5e %.5e %.5e %.5e\n", t, y[0], y[1], y[2]);
		} // end for
	gsl_odeiv2_driver_free (drv);
	//fclose (outdata);
	return 0;
}