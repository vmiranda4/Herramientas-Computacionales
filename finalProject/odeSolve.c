#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

int
func (double t, const double y[], double f[],
      void *params)
{
  (void)(t);
  double alpha = *(double *)params;
  double beta = *((double *)params + 1);
  double sigma = *((double *)params +2);
  f[0] = (1/beta - alpha)*y[0] + y[2] + y[0]*y[1];
  f[1] = -beta*y[1] - y[0]*y[0];
  f[2] = -y[0] - sigma*y[2];
  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy,
     double dfdt[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  double alpha = *(double *)params;
  double beta = *((double *)params + 1);
  double sigma = *((double *)params +2);
  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 3, 3);
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, 1/beta - alpha);
  gsl_matrix_set (m, 0, 1, y[0]);
  gsl_matrix_set (m, 0, 2, 1);
  gsl_matrix_set (m, 1, 0, -2*y[0]);
  gsl_matrix_set (m, 1, 1, -beta);
  gsl_matrix_set (m, 1, 2, 0.0);
  gsl_matrix_set (m, 2, 0, -1);
  gsl_matrix_set (m, 2, 1, 0.0);
  gsl_matrix_set (m, 2, 2, -sigma);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0;
  return GSL_SUCCESS;
}

int
main (void)
{
  const double eps_abs = 1.e-8;  /* absolute error requested  */
	const double eps_rel = 1.e-10; /* relative error requested  */
  double h = 1.0e-6;            /* starting step size for ode solver */
  //double myparams [3]={0.00001,0.1,1};
  FILE * pFile;
  pFile = fopen("parm.dat","r");
  int j;
  double x=0;
  double myparams[3];
  for(j=1; j<4; j++){
    fscanf(pFile, "%lf", &x);
    myparams[j-1]=x;
  }
  fclose(pFile);
    gsl_odeiv2_system sys = {func, jac, 3, &myparams};

  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys,  gsl_odeiv2_step_rkf45, h, eps_abs, eps_rel);
    int i;
    double t = 0.0, t1 = 1000.0;
    //double y[3] = { 0.1, 0.23, 0.31};
    FILE * vFile;
  vFile = fopen("var.dat","r");
  double y[3];
  for(j=1; j<4; j++){
    fscanf(vFile, "%lf", &x);
    y[j-1]=x;
  }
  fclose(vFile);
    FILE *output_data = fopen("coordinates.txt","w");
    fprintf (output_data,"%f %f %f %f\n", t, y[0], y[1], y[2]);	/*initial values to file*/
    for (i = 1; i <= 50000; i++)
    {
      double ti = i * t1 / 50000.0;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%d\n", status);
          break;
        }

      //printf ("%f %f %f %f\n", t, y[0], y[1], y[2]);
      fprintf (output_data,"%f %f %f %f\n", t, y[0], y[1], y[2]);
    }

  gsl_odeiv2_driver_free (d);
  return 0;
}